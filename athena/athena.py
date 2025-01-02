"""Main entrypoint to control all running of Athena"""

import polars as pl

from utils import calculate
from utils import log_handle
from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils.io import read_annotated_bed
from utils import plot
from utils.report import populate_template
from utils.util_functions import unbin


def main():
    args = parse_args()

    if args.debug:
        log_handle.setLevel("DEBUG")

    if args.annotated_bed:
        annotated_bed_file = args.annotated_bed
    else:
        annotated_bed_file = call_bedtools_intersect(
            regions=args.regions, coverage=args.coverage
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
    all_regions_plot_data = plot.all_regions(coverage_data=per_base_df)
    summary_plot = plot.gene_summary(
        gene_coverage=gene_df, threshold=args.minimum
    )

    # with pl.Config() as cfg:
    #     cfg.set_tbl_cols(100)
    #     print(exon_df)
    #     print(gene_df)

    populate_template(
        per_base_df=exon_df,
        gene_df=gene_df,
        region_df=exon_df,
        sub_threshold_plot_data=sub_threshold_plot_data,
        all_regions_plot_data=all_regions_plot_data,
        summary_plot=summary_plot,
        chromosome_plot=None,
        threshold=args.minimum,
        sample="foo",
        build=args.build,
        panel="bar",
        panel_coverage_pct=panel_coverage_pct,
    )


if __name__ == "__main__":
    main()
