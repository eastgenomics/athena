from unittest.mock import patch

import polars as pl
import pytest

from athena.utils import calculate
from tests.test_data import calculate_test_data


class TestRegionCoverage:
    pass


class TestTotalPctCoverage:
    @pytest.mark.parametrize(
        "dataframe, threshold, expected_pct",
        [
            (calculate_test_data.total_pct_coverage_df, 20, 100.00),
            (calculate_test_data.total_pct_coverage_df, 30, 50.00),
            (calculate_test_data.total_pct_coverage_df, 40, 0.00),
        ],
    )
    def test_percentage_correct_against_threshold(
        self, dataframe, threshold, expected_pct
    ):
        pct_coverage = calculate.total_pct_coverage(
            coverage_data=dataframe, threshold=threshold
        )

        assert pct_coverage == expected_pct

    def test_empty_dataframe_returns_zero(self):
        empty_dataframe = pl.DataFrame(
            schema={
                "chrom": pl.Categorical,
                "position": pl.Int32,
                "depth": pl.UInt32,
            },
        )

        calculated_pct = calculate.total_pct_coverage(
            coverage_data=empty_dataframe, threshold=20
        )

        assert calculated_pct == 0.00

    def test_coverage_at_99_99_does_not_round_to_100_pct(self):
        calculated_pct = calculate.total_pct_coverage(
            coverage_data=calculate_test_data.total_pct_coverage_df_99_99,
            threshold=20,
        )

        assert calculated_pct == 99.99


class TestMinMeanMax:
    pass


class TestPctThresholds:
    pass


class TestCalculateNormalisationFactor:
    pass


class TestMultiSampleMeanAndStdDev:
    pass


class TestNormaliseToSample:
    pass
