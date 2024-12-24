"""Plotting related functions"""

from __future__ import annotations
from timeit import default_timer as timer

import polars as pl

from .util_functions import format_timer
from utils import log_handle


def low_covered_regions(coverage_data: pl.DataFrame, threshold: int) -> str:
    """
    Generates the HTML formatted data of all exons with at least one base
    beneath given threshold depth for displaying in the low covered
    regions plots in the report.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
    threshold : int
        threshold beneath which regions are considered low coverage

    Returns
    -------
    str
        HTML formatted string representation of plot data
    """
    log_handle.debug("Generating data for low coverage regions plots")
    start = timer()

    # get the rows where depth for any position in the region under threshold
    low_coverage = (
        coverage_data.group_by("transcript", "region")
        .agg(
            pl.col("depth")
            .filter(pl.col("depth") < threshold)
            .alias("sub_threshold")
        )
        .filter(pl.col("sub_threshold") != [])
        .join(coverage_data, on=("transcript", "region"), how="left")
    )

    # format as a HTML string with transcript, positions and depth
    low_coverage = low_coverage.group_by("transcript").agg(
        pl.col("position").str.join(","),
        pl.col("depth").str.join(","),
    )
    low_coverage = low_coverage.select(
        [
            pl.format(
                "<div class='sub_plot'>{},{},{}</div>",
                "transcript",
                "position",
                "depth",
            )
        ]
    )
    low_coverage = ",".join([x[0] for x in low_coverage.rows()])

    log_handle.debug(
        "Generated plot data in %s", format_timer(start=start, end=timer())
    )

    return low_coverage


def all_regions(coverage_data: pl.DataFrame) -> list(dict):
    """
    Generates the data for plotting all regions in the report.

    Data are returned as all depths for each transcript, with plotting
    happening on the fly using Plotly in the report. This is formatted as:

    {
        'NM_000123.4': [
            (1, [36, 39, 47, 51, 43, 52, ...]),
            (2, [33, 32, 38, 44, 40, 41, ...]),
            ...
        ],
        'NM_000567.8': ...
    }

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data

    Returns
    -------
    list
        list of dicts of data per transcript
    """
    log_handle.debug('Generating plot data for all regions')

    plot_data = (
        coverage_data.group_by("transcript", "region")
        .agg("depth")
        .partition_by("transcript")
    )

    return [dict(x.rows_by_key(key="transcript")) for x in plot_data]
