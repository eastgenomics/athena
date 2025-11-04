import os
import sys
from unittest import TestCase, loader
import pathlib

sys.path.append(os.path.abspath(os.path.join(os.path.realpath(__file__), "../../")))

from bin import utils
from bin.load_data import loadData
import pandas as pd
import pandas.testing as pd_testing


class TestLoadDataFile(TestCase):
    """
    Tests for utils.load_data_file
    """

    def test_file_not_found_raises_exception(self):
        """
        Test that a FileNotFoundError is raised if the file does not exist
        """
        loader = loadData()
        with self.assertRaises(FileNotFoundError):
            loader.read_panel_bed("non_existent_file.txt")

    def test_loading_bed_with_header(self):
        """
        Test that a BED file with a header is loaded correctly
        """
        test_file_path = pathlib.Path("athena/tests/test_data/example_bed.bed")

        columns = ["chrom", "start", "end", "transcript"]

        expected_df = pd.DataFrame(
            [
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
            ],
            columns=columns,
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})

        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)

        pd_testing.assert_frame_equal(loaded_df, expected_df)

    def test_loading_bed_with_no_header(self):
        """
        Test that a BED file with no header is loaded correctly
        """
        test_file_path = pathlib.Path(
            "athena/tests/test_data/example_bed_no_header.bed"
        )

        columns = ["chrom", "start", "end", "transcript"]

        expected_df = pd.DataFrame(
            [
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
            ],
            columns=columns,
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})

        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)

        pd_testing.assert_frame_equal(loaded_df, expected_df)

    def test_loading_bed_with_chr_prefix(self):
        """
        Test that a BED file with 'chr' prefix in chromosome names is loaded
        correctly and 'chr' prefix is removed
        """
        test_file_path = pathlib.Path(
            "athena/tests/test_data/example_bed_chr_prefix.bed"
        )

        columns = ["chrom", "start", "end", "transcript"]
        expected_df = pd.DataFrame(
            [
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["1", 17345376, 17345453, "NM_003000.2"],
            ],
            columns=columns,
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        pd_testing.assert_frame_equal(loaded_df, expected_df)
