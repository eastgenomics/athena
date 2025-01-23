"""Tests for utils.calculate"""

from unittest.mock import patch

import polars as pl
import polars.testing as pl_testing
import pytest

from athena.utils import calculate

from tests.test_data.calculate import (
    min_mean_max_data,
    pct_thresholds_data,
    total_percent_coverage_data,
)


class TestRegionCoverage:
    pass


class TestTotalPctCoverage:
    """
    Data and fixture(s) for the following tests are stored in
    tests/test_data/calculate/total_percent_coverage_data.py
    """

    @pytest.mark.parametrize(
        "dataframe, threshold, expected_pct",
        [
            (total_percent_coverage_data.dataframe(), 20, 100.00),
            (total_percent_coverage_data.dataframe(), 30, 50.00),
            (total_percent_coverage_data.dataframe(), 40, 0.00),
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
            coverage_data=total_percent_coverage_data.dataframe_99_99_pct(),
            threshold=20,
        )

        assert calculated_pct == 99.99


class TestMinMeanMax:
    """
    Data and fixture(s) for the following tests are stored in
    tests/test_data/calculate/min_mean_max_data.py
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
            coverage_data=min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene",),
            join=False,
        )

        expected_values = min_mean_max_data.expected_grouped_by_gene_df()

        pl_testing.assert_frame_equal(
            calculated_values, expected_values, check_row_order=False
        )

    def test_min_mean_max_correct_when_grouped_by_gene_and_region(self):
        calculated_values = calculate.min_mean_max(
            coverage_data=min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene", "region"),
            join=False,
        )

        expected_values = (
            min_mean_max_data.expected_grouped_by_gene_and_region_df()
        )

        pl_testing.assert_frame_equal(
            calculated_values,
            expected_values,
            check_row_order=False,
        )

    def test_columns_correct_when_joined_to_input_dataframe(self):
        calculated_values = calculate.min_mean_max(
            coverage_data=min_mean_max_data.input_calculated_columns_df(),
            group_by_cols=("gene", "region"),
            join=True,
        )

        expected_values = min_mean_max_data.expected_joined_to_input_df()

        pl_testing.assert_frame_equal(
            calculated_values,
            expected_values,
            check_row_order=False,
        )


class TestPctThresholds:
    def test_value_error_raised_when_group_by_columns_not_in_dataframe(self):
        with pytest.raises(ValueError):
            calculate.pct_thresholds(
                coverage_data=pct_thresholds_data.input_df(),
                group_by_cols=("invalid_column",),
                thresholds=(10, 20, 30),
            )

    def test_type_error_raised_specified_threshold_not_an_int(self):
        with pytest.raises(TypeError):
            calculate.pct_thresholds(
                coverage_data=pct_thresholds_data.input_df(),
                group_by_cols=("gene",),
                thresholds=("10", "20"),
            )

    def test_threshold_percentages_correctly_calculated_by_gene(self):
        returned_df = calculate.pct_thresholds(
            coverage_data=pct_thresholds_data.input_df(),
            group_by_cols=("gene",),
            thresholds=(10, 20, 30),
        )

        pl_testing.assert_frame_equal(
            returned_df, pct_thresholds_data.expected_pct_coverage_per_gene()
        )

    def test_threshold_percentages_correctly_calculated_by_gene_and_region(
        self,
    ):
        returned_df = calculate.pct_thresholds(
            coverage_data=pct_thresholds_data.input_df(),
            group_by_cols=("gene", "region"),
            thresholds=(10, 20, 30),
        )

        pl_testing.assert_frame_equal(
            returned_df, pct_thresholds_data.expected_pct_coverage_per_region()
        )


class TestCalculateNormalisationFactor:

    def test_value_error_raised_when_required_columns_not_in_hsmetrics_dataframe(
        self,
    ):
        hsmetrics_df = pl.DataFrame(
            {
                "BAIT_SET": "targets",
                "BAIT_TERRITORY": "685969",
                "ON_TARGET_BASES": "870783730",
            }
        )

        with pytest.raises(ValueError):
            calculate.calculate_normalisation_factor(hsmetrics_df=hsmetrics_df)

    def test_norm_fatcor_correct_when_default_norm_value_used_if_not_specified(
        self,
    ):
        hsmetrics_df = pl.DataFrame(
            {
                "BAIT_SET": "targets",
                "BAIT_TERRITORY": "23131242345",
                "PCT_USABLE_BASES_ON_TARGET": "0.2",
                "ON_TARGET_BASES": "50000",
            }
        )

        with patch("athena.utils.calculate.NORM_VALUE", 400):
            calculated_value = calculate.calculate_normalisation_factor(
                hsmetrics_df=hsmetrics_df
            )

        assert calculated_value == 25.0

    def test_norm_factor_correct_when_passed_norm_value(self):
        hsmetrics_df = pl.DataFrame(
            {
                "BAIT_SET": "targets",
                "BAIT_TERRITORY": "23131242345",
                "PCT_USABLE_BASES_ON_TARGET": "0.2",
                "ON_TARGET_BASES": "50000",
            }
        )

        calculated_value = calculate.calculate_normalisation_factor(
            hsmetrics_df=hsmetrics_df, norm_value=1000
        )

        assert calculated_value == 50.0


class TestMultiSampleMeanAndStdDev:
    pass


class TestNormaliseToSample:
    pass
