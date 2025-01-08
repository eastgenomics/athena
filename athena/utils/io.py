"""General io related functions"""

from base64 import b64encode
from pathlib import Path
from timeit import default_timer as timer
from typing import Tuple

import polars as pl

from utils import log_handle
from .constants import DATAFRAME_TYPES
from .util_functions import unbin, format_timer


def read_file(file: Path) -> str:
    """
    Generic method to read contents of file.

    Parameters
    ----------
    file : pathlib.Path
        path to file to read from

    Returns
    -------
    str
        contents of specified file
    """
    with open(Path(file), encoding="utf-8", mode="r") as fh:
        return fh.read()


def read_image(file: Path) -> str:
    """
    Reads an image file (i.e. png) to base64 encoded string

    Parameters
    ----------
    file : Path
        path to image file to read

    Returns
    -------
    str
        base64 string of image
    """
    with open(file, "rb") as f:
        return b64encode(f.read()).decode("utf-8")


def read_annotated_bed(
    annotated_bed: Path, call_unbin: bool = False
) -> pl.DataFrame:
    """
    Read in annotated bed file with per base coverage information for
    the target regions output from `bedtools intersect`.

    Parameters
    ----------
    annotated_bed : pathlib.Path
        filename of annotated bed file
    call_unbin : bool
        Controls if to call util_functions.unbin

    Returns
    -------
    pl.DataFrame
        DataFrame of annotated bed file

    Raises
    ------
    FileNotFoundError
        Raised when given `annotated_bed` does not exist
    """
    log_handle.debug("Reading annotated bed file from %s", annotated_bed)
    start = timer()

    if not Path(annotated_bed).exists():
        raise FileNotFoundError(
            f"expected file does not exist: {annotated_bed}"
        )

    columns = [
        "chrom",
        "region_start",
        "region_end",
        "gene",
        "transcript",
        "region",
        "depth_bin_start",
        "depth_bin_end",
        "depth",
    ]

    column_types = {column: DATAFRAME_TYPES[column] for column in columns}

    coverage_data = pl.read_csv(
        source=annotated_bed,
        new_columns=columns,
        separator="\t",
        has_header=False,
        schema=column_types,
    )

    log_handle.debug(
        "Read %s rows and %s columns from bed file in %s",
        f"{coverage_data.height:,}",
        coverage_data.width,
        format_timer(start=start, end=timer()),
    )

    if call_unbin:
        coverage_data = unbin(coverage_data=coverage_data)

    return coverage_data


def read_hsmetrics(hsmetrics_file: Path) -> pl.DataFrame:
    """
    Read in contents of given hsmetrics file.

    Parameters
    ----------
    hsmetrics_file : Path
        hsmetrics file to read from
    Returns
    -------
    pl.DataFrame
        DataFrame of hsmetrics_file contents

    Raises
    ------
    AssertionError
        Raised if '### METRICS CLASS' not present in file
    """
    hsmetrics_contents = read_file(file=hsmetrics_file).splitlines()

    metrics = []

    for idx, line in enumerate(hsmetrics_contents):
        if line.startswith("## METRICS CLASS"):
            metrics.extend(hsmetrics_contents[idx + 1 : idx + 3])
            break

    assert metrics, "METRICS CLASS could not be parsed from hsmetrics file"

    return pl.DataFrame(
        [metrics[1].split("\t")], schema=metrics[0].split("\t"), orient="row"
    )


def read_sample_files(
    sample_files: Tuple[str, str],
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Convenience wrapper to call both read_annotated_bed and read_hsmetrics
    for a given sample.

    Parameters
    ----------
    sample_files : tuple
        Tuple of annotated bed and hsmetrics file to read in

    Returns
    -------
    pl.DataFrame
        DataFrame of coverage data
    pl.DataFrame
        DataFrame of hsmetrics data
    """
    return read_annotated_bed(
        annotated_bed=sample_files[0], call_unbin=True
    ), read_hsmetrics(hsmetrics_file=sample_files[1])


def write_file(file: Path, contents: str) -> None:
    """
    Generic method to write lines to file.

    Parameters
    ----------
    file : Path
        filepath to write to
    contents : str
        lines to write to file
    """
    with open(file, mode="w") as fh:
        fh.write(contents)
