"""Main entrypoint to control all running of Athena"""

import polars as pl

from utils import calculate
from utils import log_handle
from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils.io import read_annotated_bed
from utils import plot
from utils.util_functions import unbin


def main():
    args = parse_args()

    if args.debug:
        log_handle.setLevel("DEBUG")

    annotated_bed_file = call_bedtools_intersect(
        regions=args.regions, coverage=args.coverage
    )

    per_base_df = read_annotated_bed(annotated_bed=annotated_bed_file)
    per_base_df = unbin(coverage_data=per_base_df)

    generate_low_covered_regions_plot_data(
        coverage_data=per_base_df, threshold=500
    )

    gene_df, exon_df = calculate.region_coverage(
        coverage_data=per_base_df, thresholds=args.thresholds
    )

    panel_coverage_pct = calculate.total_pct_coverage(
        coverage_data=per_base_df, threshold=args.minimum
    )

    plot.low_covered_regions(coverage_data=per_base_df, threshold=args.minimum)

    with pl.Config() as cfg:
        cfg.set_tbl_cols(100)
        print(exon_df)
        print(gene_df)


if __name__ == "__main__":
    main()
