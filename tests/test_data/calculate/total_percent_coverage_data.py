"""Test data for calculate.percent_coverage"""

import polars as pl


def dataframe():
    """Minimal test input dataframe of values to calculate percentage from"""
    return pl.DataFrame(
        [["chr1", 1, 25], ["chr1", 2, 28], ["chr1", 2, 32], ["chr1", 2, 35]],
        schema={
            "chrom": pl.Categorical,
            "position": pl.Int32,
            "depth": pl.UInt32,
        },
        orient="row",
    )


def dataframe_99_99_pct():
    """DataFrame with 1/10,000 base under 20 depth => 99.99% covered"""
    return pl.DataFrame(
        [["chr1", 1, 19], *[["chr1", x, 25] for x in range(2, 10_001)]],
        schema={
            "chrom": pl.Categorical,
            "position": pl.Int32,
            "depth": pl.UInt32,
        },
        orient="row",
    )
