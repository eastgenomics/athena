"""Functions for calculating coverage values"""

from __future__ import annotations

import polars as pl


def min_mean_max(
    coverage_data: pl.DataFrame, group_by_cols: tuple, join: bool
) -> pl.DataFrame:
    """
    Calculates the min, mean and max values for all regions in the specified
    `group_by_cols` columns.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame on which to calculate values
    group_by_cols : tuple
        columns by which to group by
    join : bool
        controls if to join the calculated values back to the input DataFrame,
        if False will return the grouped by DataFrame

    Returns
    -------
    pl.DataFrame
        DataFrame with added min, mean and max values

    Raises
    ------
    ValueError
        Raised when invalid columns provided to `group_by_cols`
    """
    if not all(col in coverage_data.columns for col in group_by_cols):
        raise ValueError(
            f"Specified group_by columns {group_by_cols} not present in"
            f" dataframe, available columns: {coverage_data.columns}"
        )

    grouped_stats = coverage_data.group_by(*group_by_cols).agg(
        pl.min("depth").alias("min").cast(pl.UInt32),
        pl.mean("depth").alias("mean").cast(pl.UInt32),
        pl.max("depth").alias("max").cast(pl.UInt32),
    )

    if join:
        grouped_stats = coverage_data.join(
            grouped_stats,
            on=group_by_cols,
            how="left",
        )

    return grouped_stats


def pct_thresholds(
    coverage_data: pl.DataFrame, group_by_cols: tuple, thresholds: tuple
) -> pl.DataFrame:
    """
    Calculates the % bases at or above each of the given thresholds, adding
    each threshold as an additional column named `{threshold}x`.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame on which to calculate values
    group_by_cols : tuple
        columns by which to group by
    thresholds : tuple
        integer thresholds at which to calculate coverage

    Returns
    -------
    pl.DataFrame
        DataFrame with added `thresholds` columns
    """
    if not all(col in coverage_data.columns for col in group_by_cols):
        raise ValueError(
            f"Specified group_by columns {group_by_cols} not present in"
            f" dataframe, available columns: {coverage_data.columns}"
        )

    if not all(type(x) == int for x in thresholds):
        raise TypeError(
            f"Given threshold values not all integers: {thresholds}"
        )

    for threshold in thresholds:
        coverage_data = coverage_data.join(
            coverage_data.group_by(*group_by_cols).agg(
                ((pl.col("depth") >= threshold).sum() / pl.len())
                .alias(f"{threshold}x")
                .cast(pl.Float32)
                * 100
            ),
            on=group_by_cols,
            how="left",
        )

    return coverage_data
