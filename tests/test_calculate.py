"""Tests for utils.calculate"""

import polars as pl
import polars.testing as pl_testing
import pytest

from athena.utils import calculate
from tests.test_data import (
    calculate_min_mean_max_data,
    calculate_total_percent_coverage_data,
)


class TestRegionCoverage:
    pass


class TestTotalPctCoverage:
    """
    Data and fixture(s) for the following tests are stored in
    tests/test_data/calculate_total_percent_coverage_data.py
    """

    @pytest.mark.parametrize(
        "dataframe, threshold, expected_pct",
        [
            (calculate_total_percent_coverage_data.dataframe(), 20, 100.00),
            (calculate_total_percent_coverage_data.dataframe(), 30, 50.00),
            (calculate_total_percent_coverage_data.dataframe(), 40, 0.00),
        ],
    )
    def test_percentage_correct_against_threshold(
        self, dataframe, threshold, expected_pct
    ):
        pct_coverage = calculate.total_pct_coverage(
            coverage_data=dataframe, threshold=threshold
        )

        assert pct_coverage == expected_pct

    def test_empty_dataframe_returns_zero_percent(self):
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
            coverage_data=calculate_total_percent_coverage_data.dataframe_99_99_pct(),
            threshold=20,
        )

        assert calculated_pct == 99.99


class TestMinMeanMax:
    """
    Data and fixture(s) for the following tests are stored in
    tests/test_data/calculate_total_min_mean_max.py
    """

    def test_value_error_raised_when_invalid_groupby_columns_provided(self):
        with pytest.raises(ValueError):
            calculate.min_mean_max(
                coverage_data=pl.DataFrame(
                    {"chrom": [], "pos": [], "depth": []}
                ),
                group_by_cols="gene",
                join=False,
            )

    def test_empty_dataframe_returns_empty_dataframe_with_additional_columns(
        self,
    ):
        calculated_values = calculate.min_mean_max(
            coverage_data=pl.DataFrame({"gene": pl.Categorical, "depth": []}),
            group_by_cols=("gene",),
            join=False,
        )

        expected_df = pl.DataFrame(
            {
                "gene": [],
                "min": [],
                "mean": [],
                "max": [],
            }
        )

        pl_testing.assert_frame_equal(
            calculated_values, expected_df, check_dtypes=False
        )

    def test_min_mean_max_correct_when_groupd_by_only_gene(self):
        calculated_values = calculate.min_mean_max(
            coverage_data=calculate_min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene",),
            join=False,
        )

        expected_values = (
            calculate_min_mean_max_data.expected_grouped_by_gene_df()
        )

        pl_testing.assert_frame_equal(
            calculated_values, expected_values, check_row_order=False
        )

    def test_min_mean_max_correct_when_grouped_by_gene_and_region(self):
        calculated_values = calculate.min_mean_max(
            coverage_data=calculate_min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene", "region"),
            join=False,
        )

        expected_values = (
            calculate_min_mean_max_data.expected_grouped_by_gene_and_region_df()
        )

        pl_testing.assert_frame_equal(
            calculated_values,
            expected_values,
            check_row_order=False,
        )

    def test_columns_correct_when_joined_to_input_dataframe(self):
        calculated_values = calculate.min_mean_max(
            coverage_data=calculate_min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene", "region"),
            join=True,
        )

        expected_values = (
            calculate_min_mean_max_data.expected_joined_to_input_df()
        )

        pl_testing.assert_frame_equal(
            calculated_values,
            expected_values,
            check_row_order=False,
        )


class TestPctThresholds:
    pass


class TestCalculateNormalisationFactor:
    pass


class TestMultiSampleMeanAndStdDev:
    pass


class TestNormaliseToSample:
    pass
