"""Test data for calculate.min_mean_max"""

import polars as pl


def input_calculated_columns_df():
    """Minimal example of expected min, mean and max columns"""
    return pl.DataFrame(
        [
            ["chr1", 1, "BRCA1", "1", 20],
            ["chr1", 2, "BRCA1", "1", 30],
            ["chr1", 3, "BRCA1", "1", 40],
            ["chr1", 1, "BRCA1", "2", 10],
            ["chr1", 2, "BRCA1", "2", 15],
            ["chr1", 3, "BRCA1", "2", 20],
            ["chr1", 1, "PALB2", "1", 30],
            ["chr1", 2, "PALB2", "1", 30],
            ["chr1", 1, "PALB2", "2", 0],
        ],
        schema={
            "chrom": pl.Categorical,
            "position": pl.UInt32,
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "depth": pl.UInt32,
        },
        orient="row",
    )


def expected_grouped_by_gene_df():
    """Expected DataFrame structure when grouping by just genes"""
    return pl.DataFrame(
        [["BRCA1", 10, 22.5, 40], ["PALB2", 0, 20.0, 30]],
        schema={
            "gene": pl.Categorical,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )


def expected_grouped_by_gene_and_region_df():
    """Expected DataFrame structure when grouping just gene and regions"""
    return pl.DataFrame(
        [
            ["BRCA1", "1", 20, 30, 40],
            ["BRCA1", "2", 10, 15.0, 20],
            ["PALB2", "1", 30, 30.0, 30],
            ["PALB2", "2", 0, 0.0, 0],
        ],
        schema={
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )


def expected_joined_to_input_df():
    """Expected DataFrame structure when join=True specified"""
    return pl.DataFrame(
        [
            ["chr1", 1, "BRCA1", "1", 20, 20, 30, 40],
            ["chr1", 2, "BRCA1", "1", 30, 20, 30, 40],
            ["chr1", 3, "BRCA1", "1", 40, 20, 30, 40],
            ["chr1", 1, "BRCA1", "2", 10, 10, 15.0, 20],
            ["chr1", 2, "BRCA1", "2", 15, 10, 15.0, 20],
            ["chr1", 3, "BRCA1", "2", 20, 10, 15.0, 20],
            ["chr1", 1, "PALB2", "1", 30, 30, 30.0, 30],
            ["chr1", 2, "PALB2", "1", 30, 30, 30.0, 30],
            ["chr1", 1, "PALB2", "2", 0, 0, 0.0, 0],
        ],
        schema={
            "chrom": pl.Categorical,
            "position": pl.UInt32,
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "depth": pl.UInt32,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )
