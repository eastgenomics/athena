"""General io related functions"""

from pathlib import Path

import polars as pl

from .constants import DATAFRAME_TYPES


def read_annotated_bed(annotated_bed):
    """
    Read in annotated bed file with per base coverage information for
    the target regions output from `bedtools intersect`

    Parameters
    ----------
    annotated_bed : str
        filename of annotated bed file

    Returns
    -------
    pd.DataFrame
        DataFrame of annotated bed file

    Raises
    ------
    FileNotFoundError
        Raised when given `annotated_bed` does not exist
    """
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

    column_types = {k: v for k, v in DATAFRAME_TYPES.items() if k in columns}

    coverage_data = pl.read_csv(
        source=annotated_bed,
        new_columns=columns,
        separator="\t",
        has_header=False,
        schema=column_types,
    )

    return coverage_data
