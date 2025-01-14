"""General utility functions"""

import concurrent.futures
from multiprocessing import get_context
from os import cpu_count
from pathlib import Path
import re
from timeit import default_timer as timer
from typing import Callable, Dict, Iterable, List, Tuple

import polars as pl

from utils import log_handle
from utils.constants import DATAFRAME_TYPES


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
    log_handle.debug(
        "Unbinning data from %s rows", f"{coverage_data.height:,}"
    )
    start = timer()

    coverage_data = (
        coverage_data.with_columns(
            position=pl.int_ranges(
                start="depth_bin_start", end="depth_bin_end", dtype=pl.UInt32
            )
        )
        .drop(["depth_bin_start", "depth_bin_end"])
        .explode("position")
    )

    if (
        "region_start" in coverage_data.columns
        and "region_end" in coverage_data.columns
    ):
        coverage_data = coverage_data.filter(
            (pl.col("region_start") <= pl.col("position"))
            & (pl.col("position") < pl.col("region_end"))
        )

    log_handle.debug(
        "Completed unbinning in %s, data now has %s rows",
        format_timer(start=start, end=timer()),
        f"{coverage_data.height:,}",
    )

    return coverage_data


def call_in_parallel(
    func: Callable,
    items: Iterable,
    progress: bool = False,
    cores: int = cpu_count(),
    **kwargs,
) -> list:
    """
    Calls the given function in parallel using
    concurrent.futures.ProcessPoolExecutor on the given set of items.

    Additional arguments specified to kwargs are directly passed to the
    specified function.

    Parameters
    ----------
    func : callable
        function to call on each item
    items : list
        iterable to call function on
    progress : bool
        controls if to print progress to debug log channel
    cores : int
        no. CPU cores to split across (defaults to all available)
    Returns
    -------
    list
        list of responses
    """
    start = timer()

    log_handle.debug(
        "Calling function %s.%s for %s item(s) using %s CPU cores",
        func.__module__,
        func.__name__,
        len(items),
        cores,
    )
    results = []

    pool_executor = concurrent.futures.ProcessPoolExecutor(
        max_workers=cores, mp_context=get_context("spawn")
    )

    concurrent_jobs = {
        pool_executor.submit(func, item, **kwargs): item for item in items
    }

    n_completed = 0

    for future in concurrent.futures.as_completed(concurrent_jobs):
        # access returned output as each is returned in any order
        try:
            results.append(future.result())
            n_completed += 1

            if progress:
                log_handle.debug(
                    "Completed %s/%s processes", n_completed, len(items)
                )

        except Exception as exc:
            # catch any errors that might get raised
            print(
                "\nError calling function for input data"
                f" {concurrent_jobs[future]}: {exc}"
            )
            raise exc

    pool_executor.shutdown(wait=True)

    log_handle.debug(
        "Completed parallel calling of %s.%s in %s",
        func.__module__,
        func.__name__,
        format_timer(start=start, end=timer()),
    )

    return results


def get_column_dtypes(columns: list) -> dict:
    """
    Filter the predefined DataFrame types against given columns.

    Parameters
    ----------
    columns : list
        List of columns to get dtypes for

    Returns
    -------
    dict
        Mapping of columns to defined dtype
    """
    return {column: DATAFRAME_TYPES[column] for column in columns}


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


def pair_up_sample_files(
    hsmetrics_files: List[str], coverage_files: List[str]
) -> Dict[str, Tuple[str, str]]:
    """
    Pair up files from both lists to their file prefix.

    Ensures we have exactly one of each file for each files prefix (i.e.
    one of each file per sample).

    Parameters
    ----------
    hsmetrics_files : list
        List of hsmetrics files
    coverage_files : list
        list of per base coverage files

    Returns
    -------
    Dict[str, Tuple[str, str]]
        mapping of sample prefix to coverage file and hsmetrics file

    Raises
    ------
    ValueError
        Raised when one or more samples do not have exactly 2 files
    """
    sample_hsmetrics = {
        remove_file_extensions(file): file for file in hsmetrics_files
    }
    sample_coverage = {
        remove_file_extensions(file): file for file in coverage_files
    }
    sample_files = {
        sample: (
            sample_coverage.get(sample),
            sample_hsmetrics.get(sample),
        )
        for sample in sample_hsmetrics.keys()
    }

    missing_files = {k: v for k, v in sample_files.items() if len(v) != 2}

    if missing_files:
        raise ValueError(
            f"One or more samples with mismatched files: {missing_files}"
        )

    return sample_files


def remove_file_extensions(file: str) -> str:
    """
    Strips all file extensions from a given filename

    Parameters
    ----------
    file : str
        filename with extensions

    Returns
    -------
    str
        filename without extions
    """
    return file.replace("".join(Path(file).suffixes), "")


def strip_html_markup(html_text: str) -> str:
    """
    Strips HTML marked up text back to plain text

    Parameters
    ----------
    html_text : str
        HTML marked up text

    Returns
    -------
    str
        Input text with no markup
    """
    stripped_text = re.sub(r"<br><\/br>", "\n", html_text)
    stripped_text = re.sub(re.compile("<.*?>", re.DOTALL), "", stripped_text)

    stripped_text = "\n".join(
        [s.strip() for s in stripped_text.split("\n") if s.strip()]
    )

    return stripped_text


def natsort(dataframe: pl.DataFrame, columns: tuple) -> pl.DataFrame:
    """
    Apply natural sorting using the given columns.

    There is no current implementation of natural sorting within Polars
    (i.e. like how the natsort package can sort). This function is based
    from the following: https://github.com/pola-rs/polars/issues/17604

    Parameters
    ----------
    dataframe : pl.DataFrame
        DataFrame to sort
    columns : tuple
        columns upon which to sort

    Returns
    -------
    pl.DataFrame
        natural sorted DataFrame
    """
    return dataframe.sort(
        pl.col(*columns)
        .cast(pl.String)
        .str.extract_all(r"\D+\d*|\d+")
        .list.eval(
            pl.struct(
                string=pl.element().str.replace(r"\d+", ""),
                number=pl.element()
                .str.replace(r"\D+", "")
                .cast(pl.Int64, strict=False),
            )
        )
        .list.to_struct("max_width")
    )
