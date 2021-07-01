"""
Script to annotate a panel bed file with transcript information and per
base coverage data.

Requires: bedtools

Jethro Rainford
20/06/2021
"""

import argparse
import os
import pandas as pd
from pathlib import Path
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
        print("calling bedtools to add transcript info")

        # get total number of transcripts before to ensure none are dropped
        panel_transcripts = panel_bed.transcript.unique().tolist()

        # turn dfs into BedTools objects
        bed = bedtools.BedTool.from_dataframe(panel_bed)
        transcript_info = bedtools.BedTool.from_dataframe(transcript_info_df)

        # intersecting panel bed file with transcript/gene/exon information
        # requires 100% overlap on panel -> transcript coordinates
        bed_w_transcript = bed.intersect(
            transcript_info, wa=True, wb=True, F=1.0
        )

        # convert pybedtools object to df
        bed_w_transcript = bed_w_transcript.to_dataframe(names=[
            "p_chrom", "p_start", "p_end", "p_transcript",
            "t_chrom", "t_start", "t_end", "t_gene", "t_transcript", "t_exon"
        ])

        # check for empty file
        assert len(bed_w_transcript.index) > 0, """Empty file returned from
            intersecting panel bed and transcript file. Check if flanks are
            being used as 100% coordinate overlap is currently required."""

        # panel bed file defines transcript to use, filter transcript file for
        # just those transcripts
        bed_w_transcript = bed_w_transcript[
            bed_w_transcript["p_transcript"] == bed_w_transcript["t_transcript"]
        ]

        # drop duplicate columns
        bed_w_transcript = bed_w_transcript.drop(columns=[
            't_chrom', 't_start', 't_end', 'p_transcript'
        ])

        intersect_transcripts = bed_w_transcript.t_transcript.unique().tolist()

        # ensure no transcripts dropped from panel due to missing from
        # transcripts file
        assert len(panel_transcripts) == len(intersect_transcripts), (
            f"Transcript(s) dropped from panel during intersecting with "
            f"transcript file. Total before {len(panel_transcripts)}. Total "
            f"after {len(intersect_transcripts)}. Dropped transcripts: "
            f"{set(panel_transcripts) - set(intersect_transcripts)}"
        )

        return bed_w_transcript


    def add_coverage(self, bed_w_transcript, coverage_df, chunks=False):
        """
        Use pybedtools to add coverage bin data to selected panel regions
        Args:
            - bed_w_transcript (df): panel bed file with transcript information
            - coverage_df (df / list): coverage bin data df / list of dfs if
                chunks value passed
        Returns:
            - bed_w_coverage (df): panel bed with transcript and coverage info
        """
        print("calling bedtools to add coverage info")

        # turn dfs into BedTools objects
        bed_w_transcript = bedtools.BedTool.from_dataframe(bed_w_transcript)

        col_names = [
            "t_chrom", "t_start", "t_end", "t_gene", "t_transcript", "t_exon",
            "c_chrom", "cov_start", "cov_end", "cov"
        ]

        if not chunks:
            # per-base coverage all in one df
            coverage_df = bedtools.BedTool.from_dataframe(coverage_df)

            bed_w_coverage = bed_w_transcript.intersect(
                coverage_df, wa=True, wb=True
            )
            bed_w_coverage = bed_w_coverage.to_dataframe(names=col_names)
        else:
            # coverage data in chunks, loop over each df and intersect
            bed_w_coverage = pd.DataFrame(columns=col_names)

            for num, df in enumerate(coverage_df):
                print(f"intersecting {num + 1}/{len(coverage_df)} coverage chunks")
                # read each to bedtools object, intersect and add back to df
                chunk_df = bedtools.BedTool.from_dataframe(df)

                bed_w_coverage_chunk = bed_w_transcript.intersect(
                    chunk_df, wa=True, wb=True
                )

                bed_w_coverage_chunk = bed_w_coverage_chunk.to_dataframe(
                    names=col_names
                )

                bed_w_coverage = pd.concat(
                    [bed_w_coverage, bed_w_coverage_chunk],
                    ignore_index=True
                )

        # check again for empty output of bedtools, can happen due to memory
        # maxing out and doesn't seem to raise an exception...
        assert len(bed_w_coverage) > 0, """Error intersecting with coverage
            data, empty file generated. Is this the correct coverage data for
            the panel used? bedtools may also have reached memory limit and
            died, try re-running with --chunk_size 1000000"""

        # drop duplicate chromosome col and rename
        bed_w_coverage.drop(columns=["c_chrom"], inplace=True)

        bed_w_coverage.columns = [
            "chrom", "exon_start", "exon_end", "gene", "tx", "exon",
            "cov_start", "cov_end", "cov"
        ]

        return bed_w_coverage


def write_file(bed_w_coverage, outfile):
    """
    Write annotated bed to file
    Args:
        - bed_w_coverage (df): bed file with transcript and coverage info
        - output_prefix (str): prefix for naming output file
    Outputs: annotated_bed.tsv
    """
    # tiny function but want this separate for writing a wrapper script later
    bed_w_coverage.to_csv(outfile, sep="\t", header=False, index=False)
    print(f"annotated bed file written to {outfile}")


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
    parser.add_argument(
        '--chunk_size', '-s', type=int,
        help='number lines to read per-base coverage file in one go'
    )
    parser.add_argument(
        '--output_name', '-n',
        help='name preifx for output file, if none will use coverage file'
    )

    args = parser.parse_args()

    return args


def main():
    annotate = annotateBed()
    load = loadData()  # class of functions for reading in data

    args = parse_args()

    if not args.output_name:
        # output name not defined, use sample identifier from coverage file
        args.output_name = Path(args.coverage_file).name.split('_')[0]

    # set dir for writing to
    bin_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(bin_dir, "../output/")
    outfile_name = f"{args.output_name}_annotated.bed"
    outfile = os.path.join(out_dir, outfile_name)

    # read in files
    panel_bed_df = load.read_panel_bed(args.panel_bed)
    transcript_info_df = load.read_transcript_info(args.transcript_file)
    pb_coverage_df = load.read_coverage_data(
        args.coverage_file, args.chunk_size
    )

    # add transcript info
    bed_w_transcript = annotate.add_transcript_info(
        panel_bed_df, transcript_info_df
    )

    # add coverage
    if args.chunk_size:
        # per-base coverage split to multiple dfs to limit memory usage
        bed_w_coverage = annotate.add_coverage(
            bed_w_transcript, pb_coverage_df, chunks=True
        )
    else:
        bed_w_coverage = annotate.add_coverage(
            bed_w_transcript, pb_coverage_df, chunks=False
        )

    # sense check generated file isn't empty, should be caught earlier
    assert len(bed_w_coverage.index) > 0, (
        'An error has occured: annotated bed file is empty. This is likely ',
        'due to an error in regions defined in bed file (i.e. different ',
        'transcripts to those in the transcripts file). Start debugging by ',
        'intersecting files manually...'
    )

    write_file(bed_w_coverage, outfile)


if __name__ == "__main__":

    main()
