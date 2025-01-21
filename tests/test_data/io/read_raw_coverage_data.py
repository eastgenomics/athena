"""Test data for test_io.TestReadRawCoverage"""

import os

import polars as pl
import pytest


@pytest.fixture
def input_raw_coverage_file(tmp_path):
    """Example annotated coverage bed file from bedtools intersect"""
    target_output = os.path.join(tmp_path, "coverage.bed")

    with open(target_output, "w+") as fh:
        fh.write(
            "chr1\t2556664\t2556666\t604\n"
            "chr1\t2556666\t2556669\t605\n"
            "chr1\t2556669\t2556674\t607\n"
        )

    return target_output


def expected_raw_coverage_df():
    return pl.DataFrame(
        {
            "chrom": [
                "chr1",
                "chr1",
                "chr1",
            ],
            "depth_bin_start": [2556664, 2556666, 2556669],
            "depth_bin_end": [2556666, 2556669, 2556674],
            "depth": [604, 605, 607],
        },
        schema={
            "chrom": pl.Categorical,
            "depth_bin_start": pl.UInt32,
            "depth_bin_end": pl.UInt32,
            "depth": pl.UInt32,
        },
    )
