"""General io related functions"""

from base64 import b64encode
from datetime import datetime
import gzip
from pathlib import Path
from timeit import default_timer as timer
from typing import Tuple

import polars as pl

from .constants import NORM_VALUE
from .log import get_logger
from .util_functions import format_timer, get_column_dtypes, unbin

log_handle = get_logger("athena")

pl.enable_string_cache()


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

    Raises
    ------
    FileNotFoundError
        Raised if `file` is not a valid path
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
    EOFError
        Raised when given annotated bed is empty
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

    coverage_data = pl.read_csv(
        source=annotated_bed,
        new_columns=columns,
        separator="\t",
        has_header=False,
        schema=get_column_dtypes(columns=columns),
    )

    log_handle.debug(
        "Read %s rows and %s columns from bed file in %s",
        f"{coverage_data.height:,}",
        coverage_data.width,
        format_timer(start=start, end=timer()),
    )

    if coverage_data.height == 0:
        raise EOFError(
            f"Reading from file {annotated_bed} has returned an empty"
            " DataFrame"
        )

    if call_unbin:
        coverage_data = unbin(coverage_data=coverage_data)

    return coverage_data


def read_hsmetrics(hsmetrics_file: Path) -> pl.DataFrame:
    """
    Read in contents of given hsmetrics file to DataFrame.

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

    metrics = pl.DataFrame(
        [metrics[1].split("\t")], schema=metrics[0].split("\t"), orient="row"
    )

    return metrics


def read_normal_coverage(coverage_file: Path) -> Tuple[pl.DataFrame, int]:
    """
    Reads in the normal coverage file calculated from multiple samples.

    File is expected to contain a header line with the value used for
    normalisation stored as NORM_VALUE, along with 4 columns of chromosome,
    position, mean and 1 std deviation.

    Parameters
    ----------
    coverage_file : Path
        Path to normal coverage file

    Returns
    -------
    pl.DataFrame
        DataFrame of normal coverage
    int
        Normalisation value used for generating the normal data

    Raises
    ------
    ValueError
        Raised when fails to parse line beginning with #NORM_VALUE from
        provided coverage_file
    """
    log_handle.debug("Reading normal coverage from %s", coverage_file)

    norm_value = generated_at = generated_from = None

    with gzip.open(coverage_file, mode="rb") as fh:
        while True:
            line = fh.readline().decode("utf8").strip("\n")

            if not line.startswith("#"):
                break
            elif line.startswith("#NORM_VALUE"):
                norm_value = int(line.split("=")[1])
            elif line.startswith("#GENERATED_AT"):
                generated_at = line.split("=")[1]
            elif line.startswith("#GENERATED_FROM"):
                generated_from = line.split("=")[1]

    if not norm_value:
        raise ValueError(
            "Failed to parse #NORM_VALUE line from provided normal coverage"
            f" file {coverage_file}"
        )

    coverage_df = pl.read_csv(
        source=coverage_file,
        separator="\t",
        comment_prefix="#",
        schema=get_column_dtypes(
            columns=["chrom", "position", "obs_mean", "obs_min", "obs_max"]
        ),
    )

    log_handle.debug(
        "Read normal coverage data generated from %s samples at %s with %s"
        " positions",
        generated_from,
        generated_at,
        coverage_df.height,
    )

    return coverage_df, norm_value


def read_raw_coverage(coverage_file: Path) -> pl.DataFrame:
    """
    Reads the raw coverage (i.e. mosdepth output) into a DataFrame.

    This expects the data to be a tab separated file that is binned and
    having 4 columns consisting of chrom, bin start, bin end and depth.

    Parameters
    ----------
    coverage_file : Path
        Raw coverage file to read in

    Returns
    -------
    pl.DataFrame
        DataFrame of raw coverage
    """
    log_handle.debug("Reading raw coverage")
    coverage_df = pl.read_csv(
        source=coverage_file,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "depth_bin_start", "depth_bin_end", "depth"],
        schema={
            "chrom": pl.Categorical,
            "depth_bin_start": pl.UInt32,
            "depth_bin_end": pl.UInt32,
            "depth": pl.UInt32,
        },
    )

    log_handle.debug("finished")

    return coverage_df


def read_sample_files(
    sample_files: Tuple[str, str],
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    """
    Convenience wrapper to call both read_annotated_bed and read_hsmetrics
    for a given sample.

    Used for calculating the multi sample normal coverage.

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
    annotated_bed = read_annotated_bed(
        annotated_bed=sample_files[0], call_unbin=True
    )

    # only keep required columns to reduce memory usage
    annotated_bed = annotated_bed.select("chrom", "position", "depth")

    hsmetrics = read_hsmetrics(hsmetrics_file=sample_files[1])

    return annotated_bed, hsmetrics


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
    log_handle.debug("Writing %s lines to %s", contents.count("\n"), file)

    with open(file, mode="w") as fh:
        fh.write(contents)


def write_dataframe_to_compressed_file(
    dataframe: pl.DataFrame, filename: str
) -> None:
    """
    Writes given dataframe as compressed tsv file

    Parameters
    ----------
    dataframe : pl.DataFrame
        DataFrame to write to file
    filename : str
        Filename to write to
    """
    log_handle.debug("Writing %s lines to %s", dataframe.height, filename)

    with gzip.open(filename, mode="wb") as fh:
        dataframe.write_csv(
            file=fh,
            include_header=True,
            separator="\t",
        )


def write_normal_coverage_file(
    filename: str, coverage_df: pl.DataFrame, total_samples: int
) -> None:
    """
    Writes the multi sample dataframe of per base positions with mean
    and std deviation

    Parameters
    ----------
    filename : str
        Filename to write to
    coverage_df : pl.DataFrame
        DataFrame of coverage values to write
    total_samples : int
        Total number of samples normal generated from
    """
    log_handle.info("Writing multi sample coverage data to %s", filename)

    with gzip.open(filename, mode="wb") as fh:
        fh.write(f"#NORM_VALUE={NORM_VALUE}\n".encode())
        fh.write(
            f"#GENERATED_AT={datetime.now().strftime('%H:%M %Y-%m-%d')}\n"
            .encode()
        )
        fh.write(f"#GENERATED_FROM={total_samples} samples\n".encode())
        coverage_df.write_csv(file=fh, separator="\t", include_header=True)
