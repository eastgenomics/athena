"""Functions to handle styling of elements for the HTML output"""

import re
from typing import List, Tuple

import polars as pl

from .util_functions import natsort


def sub_threshold_regions_table(
    coverage_df: pl.DataFrame,
) -> Tuple[List[list], List[dict]]:
    """
    Styles the table of regions with < 100% coverage at the threshold.

    Parameters
    ----------
    coverage_df : pl.DataFrame
        DataFrame of regions with sub threshold coverage

    Returns
    -------
    list
        coverage data formatted as list of lists
    list
        list of dicts of column names, formatted for DataTables
    """
    threshold_columns = [
        x for x in coverage_df.columns if re.match(r"\d+x", x)
    ]

    coverage_df = coverage_df.with_columns(
        pl.col(threshold_columns + ["mean"])
        .cast(pl.Utf8)
        .str.extract(r"^(\d+\.\d{1,2})")
        .cast(pl.Float64)
    )

    coverage_df = coverage_df.rename(
        mapping={"region_start": "start", "region_end": "end"}
    )
    coverage_df = coverage_df.rename(
        mapping=lambda column: column.capitalize()
    )

    # set order for displaying
    coverage_df = coverage_df.select(
        [
            "Gene",
            "Transcript",
            "Chrom",
            "Region",
            "Start",
            "End",
            "Min",
            "Mean",
            "Max",
        ]
        + threshold_columns
    )

    coverage_df = natsort(
        dataframe=coverage_df, columns=("Transcript", "Region")
    )

    columns = [{"title": x} for x in coverage_df.columns]

    return [list(x) for x in coverage_df.rows()], columns
