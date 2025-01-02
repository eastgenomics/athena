"""General io related functions"""

from base64 import b64encode
from timeit import default_timer as timer
from pathlib import Path

import polars as pl

from .constants import DATAFRAME_TYPES
from .util_functions import format_timer
from utils import log_handle


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
    return b64encode(open(file, "rb").read()).decode("utf-8")


def read_annotated_bed(annotated_bed: Path) -> pl.DataFrame:
    """
    Read in annotated bed file with per base coverage information for
    the target regions output from `bedtools intersect`.

    Parameters
    ----------
    annotated_bed : pathlib.Path
        filename of annotated bed file

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
        "%s rows and %s columns read from bed file in %s",
        coverage_data.height,
        coverage_data.width,
        format_timer(start=start, end=timer()),
    )

    return coverage_data


def read_raw_coverage_data(file: Path) -> pl.DataFrame:
    """
    Reads the full raw coverage output bed file from samtools / mosdepth.

    Parameters
    ----------
    file : Path
        Path to file to read from

    Returns
    -------
    pl.DataFrame
        DataFrame of raw data
    """
    pass


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
    with open(file) as fh:
        fh.write(contents)
