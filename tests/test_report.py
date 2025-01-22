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
from tests.test_data.io import (
    read_annotated_bed_data,
    read_hsmetrics_data,
    read_normal_coverage_data,
    read_raw_coverage_data,
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
    pass


class TestGetTotalUniqueRegions:
    pass


class TestGetTotalFullyCoveredRegions:
    pass


class TestPopulateTemplate:
    pass
