import os
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from bin import  load, utils


TEST_DATA_DIR = (
    os.path.join(os.path.dirname(__file__), 'test_data')
)


class TestUnbin():
    """
    Tests for utils.unbin()

    Function takes the binned coverage data output from mosdepth that is
    then added to our bedfile by annotate.annotateBed(), and splits the
    bins of coverage so that each row represents a single position
    """
    # minimal binned data example
    binned_data = load.LoadData().read_annotated_bed(
        f'{TEST_DATA_DIR}/binned_data.bed')

    unbinned_data = utils.unbin(binned_data=binned_data)


    def test_position(self):
        """
        Test that positions are correctly expanded out from bins

        Bin starts at 2488096 and exon starts at 2488098 => positions
        should start at 2488098 and end at 2488105
        """
        expected_positions = list(range(2488098, 2488106))

        assert self.unbinned_data['position'].to_list() == expected_positions, (
            'Incorrect positions unbinned'
        )

    def test_coverage(self):
        """
        Test that coverage per position is correct
        """
        expected_coverage = [
            233, 236, 237, 235, 235, 235, 238, 238
        ]

        assert self.unbinned_data['cov'].to_list() == expected_coverage, (
            'Covergae in unbinned data incorrect'
        )

    def test_expected_columns(self):
        """
        Test that we get the expected columns returned
        """
        expected_columns = [
            'chrom', 'exon_start', 'exon_end', 'gene',
            'transcript', 'exon', 'cov', 'position'
        ]

        assert self.unbinned_data.columns.to_list() == expected_columns, (
            'Incorrect columns returned in unbinned data'
        )
