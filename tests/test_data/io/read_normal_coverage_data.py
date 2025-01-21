"""Test data for test_io.TestReadNormalCoverage"""

from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

from athena.utils.io import write_normal_coverage_file


@pytest.fixture
@patch("athena.utils.io.NORM_VALUE", 123456)
def input_normal_coverage_file(tmp_path):
    """Example normal coverage file contents"""
    target_output = Path(tmp_path).joinpath("normal_coverage.tsv.gz")

    example_output_values = pl.DataFrame(
        {
            "chrom": ["1", "1", "1", "1"],
            "position": [10000, 10001, 10002, 10003],
            "mean": [14.123, 16.262, 12.222, 13.333],
            "std": [1.112, 1.545, 1.234, 1.443],
        }
    )

    write_normal_coverage_file(
        filename=target_output,
        coverage_df=example_output_values,
        total_samples=32,
    )

    return target_output


def expected_normal_coverage_file_contents():
    return (
        pl.DataFrame(
            {
                "chrom": ["1", "1", "1", "1"],
                "position": [10000, 10001, 10002, 10003],
                "mean": [14.123, 16.262, 12.222, 13.333],
                "std": [1.112, 1.545, 1.234, 1.443],
            }
        ),
        123456,
    )
