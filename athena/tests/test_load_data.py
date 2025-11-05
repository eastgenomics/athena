import os
import pathlib
import sys
from unittest import TestCase

sys.path.append(os.path.abspath(os.path.join(os.path.realpath(__file__), "../../")))

import pandas as pd
import pandas.testing as pd_testing
from bin.load_data import loadData


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

    def test_loading_empty_bed_file(self):
        """
        Test that an empty BED file is loaded correctly as an empty DataFrame
        """
        test_file_path = pathlib.Path("athena/tests/test_data/example_bed_empty.bed")
        columns = ["chrom", "start", "end", "transcript"]
        expected_df = pd.DataFrame(columns=columns).astype(
            {"chrom": str, "start": "Int64", "end": "Int64", "transcript": str}
        )
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        pd_testing.assert_frame_equal(loaded_df, expected_df)

    def test_loading_bed_with_mixed_chr_prefix(self):
        """
        Test that a BED file with mixed 'chr' prefix (some with, some without)
        has all 'chr' prefixes removed
        """
        test_file_path = pathlib.Path(
            "athena/tests/test_data/example_bed_mixed_chr.bed"
        )
        columns = ["chrom", "start", "end", "transcript"]
        expected_df = pd.DataFrame(
            [
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["2", 17345376, 17345453, "NM_003001.2"],
                ["X", 17345376, 17345453, "NM_003002.2"],
                ["Y", 17345376, 17345453, "NM_003003.2"],
            ],
            columns=columns,
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        pd_testing.assert_frame_equal(loaded_df, expected_df)

    def test_loading_bed_with_different_chromosomes(self):
        """
        Test that a BED file with different chromosome types (numeric, X, Y, MT)
        is loaded correctly
        """
        test_file_path = pathlib.Path(
            "athena/tests/test_data/example_bed_diff_chroms.bed"
        )
        columns = ["chrom", "start", "end", "transcript"]
        expected_df = pd.DataFrame(
            [
                ["1", 17345376, 17345453, "NM_003000.2"],
                ["10", 17345376, 17345453, "NM_003001.2"],
                ["X", 17345376, 17345453, "NM_003002.2"],
                ["Y", 17345376, 17345453, "NM_003003.2"],
                ["MT", 17345376, 17345453, "NM_003004.2"],
            ],
            columns=columns,
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        pd_testing.assert_frame_equal(loaded_df, expected_df)

    def test_loading_bed_dtype_preservation(self):
        """
        Test that data types are correctly preserved when loading BED file
        """
        test_file_path = pathlib.Path("athena/tests/test_data/example_bed.bed")
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        assert loaded_df["chrom"].dtype == "object"  # str
        assert loaded_df["start"].dtype == "Int64"
        assert loaded_df["end"].dtype == "Int64"
        assert loaded_df["transcript"].dtype == "object"  # str

    def test_loading_bed_with_single_row(self):
        """
        Test that a BED file with a single row is loaded correctly
        """
        test_file_path = pathlib.Path(
            "athena/tests/test_data/example_bed_single_row.bed"
        )
        columns = ["chrom", "start", "end", "transcript"]
        expected_df = pd.DataFrame(
            [["1", 17345376, 17345453, "NM_003000.2"]], columns=columns
        ).astype({"chrom": str, "start": "Int64", "end": "Int64", "transcript": str})
        loader = loadData()
        loaded_df = loader.read_panel_bed(test_file_path)
        pd_testing.assert_frame_equal(loaded_df, expected_df)
