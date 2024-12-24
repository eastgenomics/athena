"""General utility functions"""

from timeit import default_timer as timer

import polars as pl

from utils import log_handle


def unbin(coverage_data: pl.DataFrame) -> pl.DataFrame:
    """
    Unbin binned coverage data to per base records in DataFrame, dropping
    those outside of region boundaries where the bin originally spanned
    the boundary. This will return a DataFrame with one row per position
    in the given regions.

    ┌───────┬──────────────┬────────────┬───┬────────────┬────────────┬───────┐
    │ chrom ┆ region_start ┆ region_end ┆ … ┆ depth_bin  ┆ depth_bin  ┆ depth │
    │ ---   ┆ ---          ┆ ---        ┆   ┆ _start     ┆ _end       ┆ ---   │
    │ ---   ┆ ---          ┆ ---        ┆   ┆ ---        ┆ ---        ┆ ---   │
    │ cat   ┆ u32          ┆ u32        ┆   ┆ u32        ┆ u32        ┆ u32   │
    ╞═══════╪══════════════╪════════════╪═══╪═══════════ ╪══════════ ═╪═══════╡
    │ 11    ┆ 108098346    ┆ 108098428  ┆ … ┆ 108098346  ┆ 108098347  ┆ 783   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ … ┆ 108098347  ┆ 108098348  ┆ 771   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ … ┆ 108098348  ┆ 108098349  ┆ 772   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ … ┆ 108098349  ┆ 108098350  ┆ 777   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ … ┆ 108098350  ┆ 108098351  ┆ 792   │
    │ …     ┆ …            ┆ …          ┆ … ┆ …          ┆ …          ┆ …     │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ … ┆ 29130709   ┆ 29130710   ┆ 1001  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ … ┆ 29130710   ┆ 29130711   ┆ 999   │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ … ┆ 29130711   ┆ 29130712   ┆ 1003  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ … ┆ 29130712   ┆ 29130713   ┆ 990   │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ … ┆ 29130713   ┆ 29130714   ┆ 993   │
    └───────┴──────────────┴────────────┴───┴────────────┴────────────┴───────┘

                                ↓

    ┌───────┬──────────────┬────────────┬─────┬───────┬───────────┐
    │ chrom ┆ region_start ┆ region_end ┆ ... | depth ┆ position  │
    │ ---   ┆ ---          ┆ ---        ┆     | ---   ┆ ---       │
    │ cat   ┆ u32          ┆ u32        ┆     | u32   ┆ u32       │
    ╞═══════╪══════════════╪════════════╪═════|═══════╪═══════════╡
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ... | 783   ┆ 108098346 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ... | 771   ┆ 108098347 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ... | 772   ┆ 108098348 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ... | 777   ┆ 108098349 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ... | 792   ┆ 108098350 │
    │ …     ┆ …            ┆ …          ┆ ... | …     ┆ …         │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ ... | 1001  ┆ 29130709  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ ... | 999   ┆ 29130710  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ ... | 1003  ┆ 29130711  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ ... | 990   ┆ 29130712  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ ... | 993   ┆ 29130713  │
    └───────┴──────────────┴────────────┴─────┴───────┴───────────┘

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame on which to calculate values

    Returns
    -------
    pd.DataFrame
        unbinned coverage data
    """
    log_handle.debug("Unbinning data from %s rows", coverage_data.height)
    start = timer()

    coverage_data = (
        coverage_data.with_columns(
            position=pl.int_ranges(
                start="depth_bin_start", end="depth_bin_end", dtype=pl.UInt32
            )
        )
        .drop(["depth_bin_start", "depth_bin_end"])
        .explode("position")
        .filter(
            (pl.col("region_start") <= pl.col("position"))
            & (pl.col("position") < pl.col("region_end"))
        )
    )

    log_handle.debug(
        "Completed unbinning in %s, data now has %s rows",
        format_timer(start=start, end=timer()),
        coverage_data.height,
    )

    return coverage_data


def format_timer(start: float, end: float) -> str:
    """
    Format timer to minutes and seconds

    Parameters
    ----------
    start : float
        start of timer
    end : float
        end of timer

    Returns
    -------
    str
        formatted total time string
    """
    return (
        f"{int(float(f'{end - start}') // 60)}m "
        f"{round(float(f'{end - start}') % 60, 2)}s"
    )
