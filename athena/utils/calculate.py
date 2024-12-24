"""Functions for calculating coverage values"""

from __future__ import annotations
from timeit import default_timer as timer

import polars as pl

from .util_functions import format_timer
from .utils import log_handle


def region_coverage(
    coverage_data: pl.DataFrame, thresholds: list
) -> tuple(pl.DataFrame, pl.DataFrame):
    """
    Calculates the coverage at both gene (transcript) and exon / intron
    level, returning both as separate DataFrames.

    Simple wrapper function to call both min_mean_max() and pct_thresholds()
    for both transcript and transcript-region levels.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
    thresholds : list
        list of thresholds to calculate coverage at

    Returns
    -------
    pl.DataFrame
        DataFrame of per gene / transcript coverage
    pl.DataFrame
        DataFrame of per exon / intron coverage
    """
    exon_df = min_mean_max(
        coverage_data=coverage_data,
        group_by_cols=("transcript", "region"),
        join=True,
    )
    exon_df = pct_thresholds(
        coverage_data=exon_df,
        group_by_cols=("transcript", "region"),
        thresholds=thresholds,
    )
    exon_df = exon_df.drop(["position", "depth"]).unique(keep="first")

    gene_df = min_mean_max(
        coverage_data=coverage_data, group_by_cols=["transcript"], join=True
    )
    gene_df = pct_thresholds(
        coverage_data=gene_df,
        group_by_cols=["transcript"],
        thresholds=thresholds,
    )
    gene_df = gene_df.drop(
        ["position", "depth", "region", "region_start", "region_end"]
    ).unique(keep="first")

    return gene_df, exon_df


def total_pct_coverage(coverage_data: pl.DataFrame, threshold: int) -> float:
    """
    Calculates the total percent coverage of all unique bases above the
    given threshold (i.e. the total panel coverage at threshold). The
    value is returned truncated to 2 dp to prevent misleading rounding
    errors (i.e. round(99.999, 2) -> 100.0).

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
    threshold : int
        threshold at which to calculate percent coverage

    Returns
    -------
    float
        total percent coverage
    """
    total_pct = (
        coverage_data.select(
            pl.col("chrom"), pl.col("position"), pl.col("depth")
        )
        .unique()
        .select(((pl.col("depth") >= threshold).sum() / pl.len()))
        .item()
    )

    return int(total_pct * 100) / 100.0


def min_mean_max(
    coverage_data: pl.DataFrame, group_by_cols: tuple, join: bool
) -> pl.DataFrame:
    """
    Calculates the min, mean and max values for all regions in the specified
    `group_by_cols` columns.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
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
    log_handle.debug(
        "Calculating min, mean and max for %s rows with column(s) %s",
        coverage_data.height,
        ", ".join(group_by_cols),
    )
    start = timer()

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

    log_handle.debug(
        "Calculated in %s", format_timer(start=start, end=timer())
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
    log_handle.debug(
        "Calculating percent thresholds for %s rows with column(s) %s against"
        " thresholds %s",
        coverage_data.height,
        ", ".join(group_by_cols),
        ", ".join(map(str, thresholds)),
    )
    start = timer()

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

    log_handle.debug(
        "Calculated in %s", format_timer(start=start, end=timer())
    )

    return coverage_data
