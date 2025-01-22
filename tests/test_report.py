"""Tests for functions related to report generation"""

import os
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import polars as pl
import polars.testing as pl_testing
import pytest

from athena.utils import report
from tests import TEST_DATA_DIR
from tests.test_data.report import (
    get_sub_threshold_regions_data,
    get_total_fully_covered_genes_data,
)


class TestGenerateSummaryText:
    pass


class TestGeneratePanelFiltrs:
    def test_value_error_raised_when_invalid_formatted_filters_provided(self):
        with pytest.raises(ValueError):
            report.generate_panel_filters(["foo:bar", "invalid"])

    def test_empty_list_returns_empty_string(self):
        returned_filters = report.generate_panel_filters(filters=[])

        assert returned_filters == ""

    def test_valid_filter_strings_returned_in_expected_format(self):
        returned_filters = report.generate_panel_filters(
            ["panel_1:BRCA1,BRCA2", "panel_2:EGFR,CFTR"]
        )

        expected_filters = (
            '<option value="BRCA1,BRCA2">panel_1</option>'
            '<option value="EGFR,CFTR">panel_2</option>'
        )

        assert returned_filters == expected_filters


class TestGetSubThresholdRegions:
    def test_value_error_raised_when_threshold_not_in_dataframe_columns(self):
        with pytest.raises(ValueError):
            report.get_sub_threshold_regions(
                region_df=get_sub_threshold_regions_data.regions_dataframe(),
                threshold=50,
            )

    def test_empty_dataframe_returned_when_all_regions_covered_to_100_pct(
        self,
    ):
        returned_df = report.get_sub_threshold_regions(
            region_df=get_sub_threshold_regions_data.regions_dataframe(),
            threshold=10,
        )

        pl_testing.assert_frame_equal(
            returned_df, get_sub_threshold_regions_data.empty_dataframe()
        )

    def test_correct_regions_returned_when_under_100_pct_at_given_threshold(
        self,
    ):
        returned_df = report.get_sub_threshold_regions(
            region_df=get_sub_threshold_regions_data.regions_dataframe(),
            threshold=30,
        )

        pl_testing.assert_frame_equal(
            returned_df,
            get_sub_threshold_regions_data.sub_30x_regions_dataframe(),
        )


class TestGetTotalUniqueRegions:
    def test_empty_dataframe_returns_zero_unique_genes_and_transcripts(self):
        empty_df = pl.DataFrame(
            {"gene": [], "transcript": []},
            schema={"gene": pl.Categorical, "transcript": pl.Categorical},
        )

        unique_genes, unique_transcripts = report.get_total_unique_regions(
            empty_df
        )

        assert unique_genes == 0 and unique_transcripts == 0

    def test_correct_unique_genes_and_transcripts_returned(self):
        regions_df = pl.DataFrame(
            {
                "gene": ["gene_1", "gene_1", "gene_2"],
                "transcript": ["transcript_1", "transcript_2", "transcript_3"],
            }
        )
        unique_genes, unique_transcripts = report.get_total_unique_regions(
            regions_df
        )

        with TestCase().subTest("unique genes"):
            assert unique_genes == 2

        with TestCase().subTest("unique transcripts"):
            assert unique_transcripts == 3


class TestGetTotalFullyCoveredRegions:
    def test_value_error_raised_when_threshold_not_in_dataframe_columns(self):
        with pytest.raises(ValueError):
            report.get_total_fully_covered_genes(
                gene_df=get_total_fully_covered_genes_data.regions_df(),
                threshold=50,
            )

    def test_empty_dataframe_returns_zero(self):
        fully_covered_genes = report.get_total_fully_covered_genes(
            gene_df=get_total_fully_covered_genes_data.empty_df(), threshold=10
        )

        assert fully_covered_genes == 0

    def test_correct_total_fully_covered_genes_returned_for_threshold(self):
        fully_covered_genes = report.get_total_fully_covered_genes(
            gene_df=get_total_fully_covered_genes_data.regions_df(),
            threshold=30,
        )

        assert fully_covered_genes == 1


class TestPopulateTemplate:
    pass
