"""
Main entrypoint to control all running of Athena
"""

from pathlib import PurePath
import sys

import polars as pl

from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils import io
from utils import calculate
from utils.util_functions import unbin


def main():
    """
    Main function to do all things Athena
    """
    args = parse_args()

    exit()

    annotated_bed_file = call_bedtools_intersect(
        regions=args.regions, coverage=args.coverage
    )

    df = io.read_annotated_bed(annotated_bed=annotated_bed_file)
    df = unbin(df)

    exon_df = calculate.min_mean_max(
        coverage_data=df, group_by_cols=("transcript", "region"), join=True
    )
    exon_df = calculate.pct_thresholds(
        coverage_data=df,
        group_by_cols=("transcript", "region"),
        thresholds=(100, 500, 1000, 1500),
    )
    exon_df = exon_df.drop(["position", "depth"]).unique(keep="first")

    gene_df = calculate.min_mean_max(
        coverage_data=df, group_by_cols=["transcript"], join=True
    )
    gene_df = calculate.pct_thresholds(
        coverage_data=df,
        group_by_cols=["transcript"],
        thresholds=(100, 500, 1000, 1500),
    )
    gene_df = gene_df.drop(
        ["position", "depth", "region", "region_start", "region_end"]
    ).unique(keep="first")

    with pl.Config() as cfg:
        cfg.set_tbl_cols(100)
        print(exon_df)
        print(gene_df)


if __name__ == "__main__":
    main()
