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
    binned_data = io.read_annotated_bed(f"{TEST_DATA_DIR}/binned_data.bed")

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
