"""Main entrypoint to control all running of Athena"""

import argparse
from timeit import default_timer as timer

from utils import calculate

from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils.calculate import normalise_to_sample
from utils.io import (
    read_annotated_bed,
    read_hsmetrics,
    read_normal_coverage,
    read_raw_coverage,
    read_sample_files,
    write_dataframe_to_compressed_file,
    write_file,
    write_multi_sample_coverage,
)
from utils import plot
from utils.report import generate_summary_text, populate_template
from utils.util_functions import (
    call_in_parallel,
    format_timer,
    pair_up_sample_files,
    strip_html_markup,
    unbin,
)
from version import VERSION

from utils.log import get_logger

log_handle = get_logger("athena")


def generate_report(args: argparse.Namespace) -> None:
    """
    Call all methods for calculating coverage and generating report

    Parameters
    ----------
    args : argparse.Namespace
        Command line argument Namespace object
    """
    print("Beginning generating coverage stats and coverage report")
    start = timer()

    annotated_bed_file = call_bedtools_intersect(
        regions=args.regions,
        coverage=args.coverage,
        overwrite=args.force,
        build=args.build,
    )

    per_base_df = read_annotated_bed(annotated_bed=annotated_bed_file)
    per_base_df = unbin(coverage_data=per_base_df)

    # generate stats
    gene_df, region_df = calculate.region_coverage(
        coverage_data=per_base_df, thresholds=args.thresholds
    )
    panel_coverage_pct = calculate.total_pct_coverage(
        coverage_data=per_base_df, threshold=args.minimum
    )

    if args.normal_coverage:
        hsmetrics_df = read_hsmetrics(hsmetrics_file=args.hsmetrics)
        normal_coverage_df, norm_value = read_normal_coverage(
            coverage_file=args.normal_coverage
        )

        normal_coverage_df = normalise_to_sample(
            normal_coverage=normal_coverage_df,
            hsmetrics=hsmetrics_df,
            norm_value=norm_value,
        )

        per_base_df = per_base_df.join(
            normal_coverage_df, how="left", on=["chrom", "position"]
        )

    # generate plots
    sub_threshold_plot_data = all_region_plots = summary_plot = (
        summary_text
    ) = chromosome_plots = "null"

    summary_plot = plot.gene_summary(
        gene_coverage=gene_df, threshold=args.minimum
    )

    if args.plot_sub_threshold:
        sub_threshold_plot_data = plot.sub_threshold_regions(
            coverage_data=per_base_df, threshold=args.minimum
        )

    if (
        args.limit == -1
        or gene_df.select("transcript").unique().height <= args.limit
    ):
        all_region_plots = plot.all_regions(
            coverage_data=per_base_df, threshold=args.minimum
        )

    if args.plot_chromosomes:
        raw_coverage = read_raw_coverage(coverage_file=args.coverage)
        chromosome_plots = plot.all_chromosomes(raw_coverage=raw_coverage)

    if args.summary:
        summary_text = generate_summary_text(
            gene_df=gene_df,
            threshold=args.minimum,
            panel_coverage_pct=panel_coverage_pct,
            indication=args.clinical_indication,
        )

    populated_report = populate_template(
        summary_text=summary_text,
        gene_df=gene_df,
        region_df=region_df,
        sub_threshold_plot_data=sub_threshold_plot_data,
        all_region_plots=all_region_plots,
        summary_plot=summary_plot,
        chromosome_plots=chromosome_plots,
        threshold=args.minimum,
        sample=args.output,
        build=args.build,
        panel=args.panel,
        panel_coverage_pct=panel_coverage_pct,
        panel_filters=args.panel_filters,
        version=VERSION,
    )

    output_file = f"{args.output}_coverage_report.html"
    write_file(file=output_file, contents=populated_report)

    if args.write_data:
        write_dataframe_to_compressed_file(
            dataframe=per_base_df, filename=f"{args.output}.coverage.bed.gz"
        )
        write_dataframe_to_compressed_file(
            dataframe=region_df,
            filename=f"{args.output}.region_coverage.tsv.gz",
        )
        write_dataframe_to_compressed_file(
            dataframe=gene_df, filename=f"{args.output}.gene_coverage.tsv.gz"
        )

    if args.summary_file:
        write_file(
            file=f"{args.output}_summary.txt",
            contents=strip_html_markup(summary_text),
        )

    print(
        "Completed all steps in"
        f" {format_timer(start=start, end=timer())}. Report written to"
        f" {output_file}",
    )


def generate_multi_sample_coverage(args: argparse.Namespace) -> None:
    """
    Calculates mean per base coverage and standard deviation from the
    mean for all provided samples. This is to generate a file to define
    'normal' coverage for adding context to whole gene plots in the
    output report.

    Parameters
    ----------
    args : argparse.Namespace
        Command line argument Namespace object
    """
    log_handle.info(
        "Calculating multi sample coverage from %s samples", len(args.coverage)
    )

    annotated_beds = call_in_parallel(
        call_bedtools_intersect,
        items=args.coverage,
        progress=True,
        regions=args.regions,
        build=args.build,
        overwrite=True,
    )

    sample_files = pair_up_sample_files(
        first_file_list=annotated_beds, second_file_list=args.hsmetrics
    )

    sample_dfs = call_in_parallel(read_sample_files, sample_files.values())

    normalised_coverage_df = calculate.multi_sample_mean_and_std_dev(
        sample_dfs=sample_dfs
    )

    write_multi_sample_coverage(
        filename=f"{args.output}.tsv.gz",
        coverage_df=normalised_coverage_df,
        total_samples=len(sample_files),
    )

    log_handle.info("Completed calculating multi sample coverage.")


def main():
    args = parse_args()

    if args.verbose:
        log_handle.setLevel("DEBUG")

    if args.mode == "report":
        generate_report(args=args)
    else:
        generate_multi_sample_coverage(args=args)


if __name__ == "__main__":
    main()
