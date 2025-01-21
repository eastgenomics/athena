"""Test data for test_io.TestReadAnnotatedBed"""

from pathlib import Path

import polars as pl
import pytest

from athena.utils import constants


@pytest.fixture
def input_coverage_bed_file(tmp_path):
    """Example annotated coverage bed file from bedtools intersect"""
    target_output = Path(tmp_path).joinpath("coverage.bed")

    with open(target_output, "w+") as fh:
        fh.write(
            "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t2556664\t2556666\t604\n"
            "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t2556666\t2556668\t605\n"
            "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t2556668\t2556669\t607\n"
        )

    return target_output


def expected_binned_coverage_bed_df():
    """Expected data frame when unbin=False"""
    required_cols = [
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

    # ensure we use the same predefined dtypes
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}

    return pl.DataFrame(
        [
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                2556664,
                2556666,
                604,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                2556666,
                2556668,
                605,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                2556668,
                2556669,
                607,
            ],
        ],
        schema=schema,
        orient="row",
    )


def expected_unbinned_coverage_bed_df():
    """Expected data frame when unbin=True"""

    required_cols = [
        "chrom",
        "region_start",
        "region_end",
        "gene",
        "transcript",
        "region",
        "depth",
        "position",
    ]

    # ensure we use the same predefined dtypes
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}

    return pl.DataFrame(
        [
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                604,
                2556664,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                604,
                2556665,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                605,
                2556666,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                605,
                2556667,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                607,
                2556668,
            ],
        ],
        schema=schema,
        orient="row",
    )
