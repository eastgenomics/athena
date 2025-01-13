"""Functions for calculating coverage values"""

from __future__ import annotations
from timeit import default_timer as timer
from typing import List, Tuple

import polars as pl

from utils import log_handle
from .constants import NORM_VALUE
from .util_functions import format_timer


def region_coverage(
    coverage_data: pl.DataFrame, thresholds: list
) -> Tuple[pl.DataFrame, pl.DataFrame]:
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
    given threshold (i.e. the total panel coverage at threshold).

    The value is returned truncated to 2 dp to prevent misleading rounding
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
    ) * 100

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
        f"{coverage_data.height:,}",
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
        pl.mean("depth").alias("mean").cast(pl.Float32),
        pl.max("depth").alias("max").cast(pl.UInt32),
    )

    if join:
        grouped_stats = coverage_data.join(
            grouped_stats,
            on=group_by_cols,
            how="left",
        )

    log_handle.debug(
        "Calculated min, mean and max in %s",
        format_timer(start=start, end=timer()),
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
        f"{coverage_data.height:,}",
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
        "Calculated percent thresholds in %s",
        format_timer(start=start, end=timer()),
    )

    return coverage_data


def calculate_normalisation_factor(
    hsmetrics_df: pl.DataFrame, norm_value: int
) -> int:
    """
    Calculates the factor for which to normalise against. This will use
    values from the hsmetrics file and the provided normalisation value.

    Parameters
    ----------
    hsmetrics_df : pl.DataFrame
        DataFrame of hsmetrics values
    norm_value : int
        Normalisation value to use

    Returns
    -------
    int
        Normalisation factor
    """
    sample_bases = hsmetrics_df.select(
        pl.col("ON_TARGET_BASES").cast(pl.Int32)
        * pl.col("PCT_USABLE_BASES_ON_TARGET").cast(pl.Float64)
    ).item()

    return norm_value / sample_bases


def multi_sample_mean_and_std_dev(
    sample_dfs: List[Tuple[pl.DataFrame, pl.DataFrame]],
) -> pl.DataFrame:
    """
    Calculates the normalised mean and std deviation across all positions.

    Normalisation is calculated as 1,000,000 over the fraction of usable,
    de-deduplicated on target bases.

    Parameters
    ----------
    sample_dfs : List[pl.DataFrame, pl.DataFrame]
        List of per sample coverage and hsmetrics dataframes

    Returns
    -------
    pl.DataFrame
        DataFrame of mean and std deviation per position
    """
    start = timer()
    log_handle.debug(
        "Calculating mean and std deviation across %s samples with"
        " normalisation value of %s",
        len(sample_dfs),
        NORM_VALUE,
    )

    combined_coverage_df = sample_dfs[0][0].select(
        pl.col("chrom"), pl.col("position")
    )
    sample_columns = []

    for idx, dfs in enumerate(sample_dfs):
        coverage_df, hsmetrics_df = dfs
        sample_columns.append(str(idx))

        norm_factor = calculate_normalisation_factor(
            hsmetrics_df=hsmetrics_df, norm_value=NORM_VALUE
        )

        coverage_df = coverage_df.with_columns(
            (pl.col("depth") * norm_factor)
        ).rename({"depth": str(idx)})

        combined_coverage_df = combined_coverage_df.join(
            coverage_df, on=["chrom", "position"], how="left"
        )

    combined_coverage_df = (
        combined_coverage_df.with_columns(
            pl.concat_list(sample_columns).alias("all")
        )
        .drop(sample_columns)
        .with_columns(
            pl.col("all").list.mean().alias("mean"),
            pl.col("all").list.std().alias("std"),
        )
        .drop("all")
    )

    log_handle.debug(
        "Completed calculating mean and std dev in %s",
        format_timer(start=start, end=timer()),
    )

    return combined_coverage_df


def normalise_to_sample(
    normal_coverage: pl.DataFrame, hsmetrics: pl.DataFrame, norm_value: int
) -> pl.DataFrame:
    """
    Normalises the normal coverage values to the given sample.

    This will normalise the normal coverage against the amount of
    sequencing for the given sample, adjusting the normal for the amount
    of given sequencing.

    This will add the normalised mean and +/- 3 std deviations as
    separate columns to the returned dataframe.

    Parameters
    ----------
    normal_coverage : pl.DataFrame
        Per base dataframe of normal coverage
    hsmetrics : pl.DataFrame
        DataFrame of hsmetrics for sample
    norm_value : int
        Normalisation value to use, required to be same value used for
        generating the normal data

    Returns
    -------
    pl.DataFrame
        Per base dataframe of normal coverage, normalised to sample
    """
    norm_factor = calculate_normalisation_factor(
        hsmetrics_df=hsmetrics, norm_value=norm_value
    )

    normal_coverage = (
        normal_coverage.with_columns(
            (pl.col("mean") * norm_factor).alias("normal_mean"),
            (pl.col("std") * norm_factor).alias("normal_std"),
        )
        .drop("mean", "std")
        .with_columns(
            pl.col("normal_mean"),
            (pl.col("normal_mean") - (pl.col("normal_std")) * 3).alias(
                "mean_-_std"
            ),
            (pl.col("normal_mean") + (pl.col("normal_std")) * 3).alias(
                "mean_+_std"
            ),
        )
    )

    return normal_coverage
