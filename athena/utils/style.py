"""Functions to handle styling of elements for the HTML output"""

import re
from typing import List, Tuple

import polars as pl

from .util_functions import natsort


def dataframe_for_html(
    coverage_df: pl.DataFrame, sort_by: tuple
) -> Tuple[List[list], List[dict]]:
    """
    Styles the table of regions with < 100% coverage at the threshold.

    Parameters
    ----------
    coverage_df : pl.DataFrame
        DataFrame of regions with sub threshold coverage
    sort_by : tuple
        Columns by which to sort the dataframe

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

    # limit floats to 2 dp
    coverage_df = coverage_df.with_columns(
        pl.col(threshold_columns + ["mean"])
        .cast(pl.Utf8)
        .str.extract(r"^(\d+\.\d{1,2})")
        .cast(pl.Float64)
    )

    if "region_start" in coverage_df.columns:
        coverage_df = coverage_df.rename(
            mapping={"region_start": "start", "region_end": "end"}
        )

    coverage_df = coverage_df.rename(
        mapping=lambda column: column.capitalize()
    )

    # set order for displaying
    column_order = [
        "Gene",
        "Transcript",
        "Chrom",
        "Region",
        "Start",
        "End",
        "Min",
        "Mean",
        "Max",
    ] + threshold_columns

    coverage_df = coverage_df.select(
        [x for x in column_order if x in coverage_df.columns]
    )

    coverage_df = natsort(dataframe=coverage_df, columns=(sort_by))

    columns = [{"title": x} for x in coverage_df.columns]

    return [list(x) for x in coverage_df.rows()], columns
