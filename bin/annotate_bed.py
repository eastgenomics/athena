"""
Script to annotate a panel bed file with per base coverage data.

Jethro Rainford
20/06/2021
"""

import argparse
import os
import pybedtools as bedtools

from load_data import loadData


class annotateBed():

    def add_transcript_info(self, panel_bed, transcript_info_df):
        """
        Use pybedtools to annotate panel bed file with coverage data
        Args:
            - panel_bed (df): panel bed file regions df
            - transcript_info_df (df): transcript info file df
        Returns:
            - bed_w_transcript (df): panel bed file with transcript information
        """
        print("Calling bedtools to add transcript info")

        # turn dfs into BedTools objects
        bed = bedtools.BedTool.from_dataframe(panel_bed)
        transcript_info = bedtools.BedTool.from_dataframe(transcript_info_df)

        # intersecting panel bed file with transcript/gene/exon information
        # requires regions overlap 100% due to potential for overlapping regions
        bed_w_transcript = bed.intersect(
            transcript_info, wb=True, f=1.0, r=True
        )

        # convert pybedtools object to df
        bed_w_transcript = bed_w_transcript.to_dataframe(names=[
            "p_chrom", "p_start", "p_end", "p_transcript",
            "t_chrom", "t_start", "t_end", "t_gene", "t_transcript", "t_exon"
        ])

        # panel bed file defines transcript to use, filter transcript file for
        # just those transcripts
        bed_w_transcript = bed_w_transcript[
            bed_w_transcript["p_transcript"] == bed_w_transcript["t_transcript"]
        ]

        # drop duplicate columns
        bed_w_transcript = bed_w_transcript.drop(columns=[
            'p_chrom', 'p_start', 'p_end', 'p_transcript'
        ])

        return bed_w_transcript


    def add_coverage(self, bed_w_transcript, coverage_df):
        """
        Use pybedtools to add coverage bin data to selected panel regions
        Args:
            - bed_w_transcript (df): panel bed file with transcript information
            - coverage_df (df): coverage bin data df
        Returns:
            - bed_w_coverage (df): panel bed with transcript and coverage info
        """
        print("Calling bedtools to add coverage info")

        # turn dfs into BedTools objects
        bed_w_transcript = bedtools.BedTool.from_dataframe(bed_w_transcript)
        coverage_df = bedtools.BedTool.from_dataframe(coverage_df)

        bed_w_coverage = bed_w_transcript.intersect(
            coverage_df, wa=True, wb=True
        )

        bed_w_coverage = bed_w_coverage.to_dataframe(names=[
            "t_chrom", "t_start", "t_end", "t_gene", "t_transcript", "t_exon",
            "c_chrom", "cov_start", "cov_end", "cov"
        ])

        # drop duplicate chromosome col and rename
        bed_w_coverage.drop(columns=["c_chrom"], inplace=True)

        bed_w_coverage.columns = [
            "chrom", "exon_start", "exon_end", "gene", "tx", "exon",
            "cov_start", "cov_end", "cov"
        ]

        return bed_w_coverage


def write_file(bed_w_coverage):
    """
    Write annottated bed to file
    Args: bed_w_coverage (df): bed file with transcript and coverage info
    Outputs: annotated_bed.tsv
    """
    bin_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(bin_dir, "../output/")
    outfile = os.path.join(out_dir, "annotated_bed.tsv")

    bed_w_coverage.to_csv(outfile, sep="\t", index=False)


def parse_args():
    """
    Parse cmd line arguments

    Args: None

    Returns:
        - args (arguments): args passed from cmd line
    """
    parser = argparse.ArgumentParser(
        description='Annotate panel bed file with transcript & coverage data.'
    )
    parser.add_argument(
        '--panel_bed', '-p',
        help='panel bed file'
    )
    parser.add_argument(
        '--transcript_file', '-t',
        help='file with gene and exon information'
    )
    parser.add_argument(
        '--coverage_file', '-c',
        help='per base coverage data file'
    )

    args = parser.parse_args()

    return args


def main():
    annotate = annotateBed()
    load = loadData()

    args = parse_args()

    # read in files
    panel_bed = load.read_panel_bed(args.panel_bed)
    transcript_info_df = load.read_transcript_info(args.transcript_file)
    pb_coverage_df = load.read_coverage_data(args.coverage_file)

    # add transcript info
    bed_w_transcript = annotate.add_transcript_info(
        panel_bed, transcript_info_df
    )

    # add coverage
    bed_w_coverage = annotate.add_coverage(bed_w_transcript, pb_coverage_df)

    # sense check generated file isn't empty
    assert len(bed_w_coverage) > 0, (
        'An error has occured: annotated bed file is empty. This is likely ',
        'due to an error in regions defined in bed file (i.e. different ',
        'transcripts to those in the transcripts file). Start debugging by ',
        'intersecting files manually...'
    )

    write_file(bed_w_coverage)


if __name__ == "__main__":

    main()
