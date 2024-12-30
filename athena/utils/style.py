"""Functions to handle styling of elements for the HTML output"""

import re
from typing import List

import polars as pl


def sub_threshold_regions_table(coverage_df: pl.DataFrame) -> List[list]:
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
    """
    float_columns = [
        x for x in coverage_df.columns if re.match(r"\d+x", x)
    ] + ["mean"]

    coverage_df = coverage_df.with_columns(
        pl.col(float_columns)
        .cast(pl.Utf8)
        .str.extract(r"^(\d+\.\d{1,2})")
        .cast(pl.Float64)
    )

    return [list(x) for x in coverage_df.rows()]
