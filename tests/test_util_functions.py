from unittest.mock import patch

import pytest

from athena.utils import io, util_functions
from tests import TEST_DATA_DIR


class TestUnbin:
    """
    Tests for util_functions.unbin()

    Function takes the binned coverage data output from mosdepth that is
    then added to our bedfile by annoatate.call_bedtools_intersect, and
    splits the bins of coverage so that each row represents a single
    position.
    """

    # minimal binned data example
    binned_data = io.read_annotated_bed(
        f"{TEST_DATA_DIR}/files/binned_data.bed"
    )

    unbinned_data = util_functions.unbin(coverage_data=binned_data)

    def test_ubinned_positions_have_correct_range(self):
        """
        Test that positions are correctly expanded out from bins

        Bin starts at 2488096 and exon starts at 2488098 => positions
        should start at 2488098 and end at 2488105
        """
        expected_positions = list(range(2488098, 2488106))

        assert (
            self.unbinned_data["position"].to_list() == expected_positions
        ), "Incorrect positions unbinned"

    def test_coverage_values_as_expected_from_bins(self):
        """
        Test that coverage per position is correct
        """
        expected_coverage = [233, 236, 237, 235, 235, 235, 238, 238]

        assert (
            self.unbinned_data["depth"].to_list() == expected_coverage
        ), "Coverage in unbinned data incorrect"

    def test_expected_columns_after_unbinning(self):
        """
        Test that we get the expected columns returned
        """
        expected_columns = [
            "chrom",
            "region_start",
            "region_end",
            "gene",
            "transcript",
            "region",
            "depth",
            "position",
        ]

        assert sorted(self.unbinned_data.columns) == sorted(
            expected_columns
        ), "Incorrect columns returned in unbinned data"


class TestCompressAndEncode:
    @pytest.mark.parametrize(
        "input,expected",
        [
            ("foo", "eNpLy88HAAKCAUU="),
            ("foo bar", "eNpLy89XSEosAgAKcAKa"),
            ([["foo"], ["bar"]], "eNqLjlZPy89Xj9VRiFZPSixSj40FADS7BYo="),
        ],
    )
    def test_base64_compressed_string_returned(self, input, expected):
        compressed_data = util_functions.compress_and_encode(data=input)

        assert compressed_data == expected


class TestGetColumnTypes:
    def test_correct_columns_and_types_returned(self):
        all_defined_dtypes = {
            "chrom": str,
            "pos": int,
            "start": int,
            "end": int,
            "cov": float,
        }

        with patch(
            "athena.utils.util_functions.DATAFRAME_TYPES", all_defined_dtypes
        ):
            selected_types = util_functions.get_column_dtypes(
                ["chrom", "pos", "cov"]
            )

            expected_types = {"chrom": str, "pos": int, "cov": float}

            assert selected_types == expected_types

    def test_key_error_raised_on_passing_non_defined_column(self):
        all_defined_dtypes = {
            "chrom": str,
            "pos": int,
            "start": int,
            "end": int,
            "cov": float,
        }

        with patch(
            "athena.utils.util_functions.DATAFRAME_TYPES", all_defined_dtypes
        ) and pytest.raises(KeyError):
            util_functions.get_column_dtypes(["chrom", "foo"])


class TestFormatTimer:
    @pytest.mark.parametrize(
        "start,end,expected_time",
        [
            (1.0, 10.0, "0m 9.0s"),
            (1.0, 1.55, "0m 0.55s"),
            (1.0, 62.0, "1m 1.0s"),
            (1.0, 666.666, "11m 5.67s"),
        ],
    )
    def test_time_delta_formatted_correctly_as_string(
        self, start, end, expected_time
    ):
        pretty_time = util_functions.format_timer(start=start, end=end)

        assert pretty_time == expected_time


class TestPairUpSampleFiles:
    def test_file_pairs_correctly_returned_when_all_are_paired(self):
        first_sample_files = ["sample_1.bed.gz", "sample_2.bed.gz"]
        second_sample_files = ["sample_1.per-bed.gz", "sample_2.per-bed.gz"]

        expected_paired_files = {
            "sample_1": ("sample_1.bed.gz", "sample_1.per-bed.gz"),
            "sample_2": ("sample_2.bed.gz", "sample_2.per-bed.gz"),
        }

        actual_paired_files = util_functions.pair_up_sample_files(
            first_file_list=first_sample_files,
            second_file_list=second_sample_files,
        )

        assert expected_paired_files == actual_paired_files

    def test_value_error_raised_when_sample_missing_file(self):
        first_sample_files = ["sample_1.bed.gz", "sample_2.bed.gz"]
        second_sample_files = ["sample_1.per-bed.gz"]

        with pytest.raises(ValueError):
            util_functions.pair_up_sample_files(
                first_file_list=first_sample_files,
                second_file_list=second_sample_files,
            )

    def test_value_error_raised_when_sample_has_more_than_two_files(self):
        first_sample_files = ["sample_1.bed.gz", "sample_2.bed.gz"]
        second_sample_files = [
            "sample_1.per-bed.gz",
            "sample_2.per-bed.gz",
            "sample_2.bonus.bed",
        ]

        with pytest.raises(ValueError):
            util_functions.pair_up_sample_files(
                first_file_list=first_sample_files,
                second_file_list=second_sample_files,
            )


class TestRemoveFileExtension:
    @pytest.mark.parametrize(
        "file,expected_str",
        [
            ("sample_1.bed.gz", "sample_1"),
            ("sample_1.foo.bar.baz.gz", "sample_1"),
            ("sample_1_bed.gz", "sample_1_bed"),
        ],
    )
    def test_all_suffixes_correctly_removed(self, file, expected_str):
        unsuffixed_file = util_functions.remove_file_extensions(file=file)

        assert unsuffixed_file == expected_str
