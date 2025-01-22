""" "Data for report.get_sub_threshold_regions tests"""

import polars as pl


def empty_dataframe() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chr": [],
            "pos": [],
            "10x": [],
            "30x": [],
        },
        schema={
            "chr": pl.Categorical,
            "pos": pl.UInt32,
            "10x": pl.Float32,
            "30x": pl.Float32,
        },
    )


def regions_dataframe() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chr": ["chr1", "chr2"],
            "pos": [123, 456],
            "10x": [100.0, 100.0],
            "30x": [100.0, 95.5],
        },
        schema={
            "chr": pl.Categorical,
            "pos": pl.UInt32,
            "10x": pl.Float32,
            "30x": pl.Float32,
        },
    )


def sub_30x_regions_dataframe() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chr": ["chr2"],
            "pos": [456],
            "10x": [100.0],
            "30x": [95.5],
        },
        schema={
            "chr": pl.Categorical,
            "pos": pl.UInt32,
            "10x": pl.Float32,
            "30x": pl.Float32,
        },
    )
