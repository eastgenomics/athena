"""Main entrypoint to control all running of Athena"""

from timeit import default_timer as timer

from utils import calculate
from utils import log_handle
from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils.io import read_annotated_bed, write_file
from utils import plot
from utils.report import generate_summary_text, populate_template
from utils.util_functions import format_timer, strip_html_markup, unbin


def main():
    start = timer()
    args = parse_args()

    print("Beginning generating coverage stats and coverage report")

    if args.debug:
        log_handle.setLevel("DEBUG")

    if args.annotated_bed:
        annotated_bed_file = args.annotated_bed
    else:
        annotated_bed_file = call_bedtools_intersect(
            regions=args.regions, coverage=args.coverage, overwrite=args.force
        )

    per_base_df = read_annotated_bed(annotated_bed=annotated_bed_file)
    per_base_df = unbin(coverage_data=per_base_df)

    # generate stats
    gene_df, exon_df = calculate.region_coverage(
        coverage_data=per_base_df, thresholds=args.thresholds
    )
    panel_coverage_pct = calculate.total_pct_coverage(
        coverage_data=per_base_df, threshold=args.minimum
    )

    # generate plots
    sub_threshold_plot_data = plot.sub_threshold_regions(
        coverage_data=per_base_df, threshold=args.minimum
    )
    all_region_plots = plot.all_regions(
        coverage_data=per_base_df, threshold=args.minimum
    )
    summary_plot = plot.gene_summary(
        gene_coverage=gene_df, threshold=args.minimum
    )

    if args.summary:
        summary_text = generate_summary_text(
            gene_df=gene_df,
            threshold=args.minimum,
            panel_coverage_pct=panel_coverage_pct,
            indication=args.clinical_indication,
        )
    else:
        summary_text = ""

    populated_report = populate_template(
        summary_text=summary_text,
        per_base_df=exon_df,
        gene_df=gene_df,
        region_df=exon_df,
        sub_threshold_plot_data=sub_threshold_plot_data,
        all_region_plots=all_region_plots,
        summary_plot=summary_plot,
        chromosome_plot=None,
        threshold=args.minimum,
        sample=args.output,
        build=args.build,
        panel=args.panel,
        panel_coverage_pct=panel_coverage_pct,
        panel_filters=args.panel_filters,
    )

    output_file = f"{args.output}_coverage_report.html"
    write_file(file=output_file, contents=populated_report)

    per_base_df.write_csv(
        file=f"{args.output}.coverage.bed", include_header=True, separator="\t"
    )
    gene_df.write_csv(
        file=f"{args.output}.gene_coverage.tsv",
        include_header=True,
        separator="\t",
    )
    exon_df.write_csv(
        file=f"{args.output}.region_coverage.tsv",
        include_header=True,
        separator="\t",
    )

    if args.summary_file:
        write_file(
            file=f"{args.output}_summary.txt",
            contents=strip_html_markup(summary_text),
        )

    print(
        "Completed generating report in"
        f" {format_timer(start=start, end=timer())}. Report written to"
        f" {output_file}",
    )


if __name__ == "__main__":
    main()
