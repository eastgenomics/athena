"""Plotting related functions"""

from timeit import default_timer as timer

import polars as pl

from .util_functions import format_timer
from utils import log_handle


def generate_low_covered_regions_plot_data(
    coverage_data: pl.DataFrame, threshold: int
) -> str:
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
