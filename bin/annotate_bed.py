"""
Script to annotate a panel bed file with per base coverage data.

Jethro Rainford
20/06/2021
"""

import argparse
import math
import multiprocessing
import numpy as np
import os
import pandas as pd
from pathlib import Path
import pybedtools as bedtools
import re
import sys

import load_data


class annotateBed():

    def add_transcript_info(panel_bed, transcript_info_df):
        """
        Use pybedtools to annotate panel bed file with coverage data
        """
        # turn dfs into BedTools objects
        bed = bedtools.BedTool.from_dataframe(panel_bed)
        transcript_info = bedtools.BedTool.from_dataframe(transcript_info_df)

        bed_w_transcript = bed.intersect(
            transcript_info_df, wa=True, wb=True
        )


def main():
    # load files into dfs
    annotate = annotateBed()
    load = load_data()

    args = parse_args

    panel_bed = load.read_panel_bed(args.panel_bed)
    transcript_info_df = load.read_transcript_info(args.transcript_file)
    pb_coverage_df = load.read_coverage_data(args.coverage_file)



def parse_args():
    """
    Parse cmd line arguments

    Args: None

    Returns:
        - args (arguments): args passed from cmd line
    """
    parser = argparse.ArgumentParser(
        description='Annotate panel bed file with transcript information and coverage data.'
    )
    parser.add_argument(
        '--panel_bed',
        help='panel bed file'
    )
    parser.add_argument(
        '--transcript_file',
        help='file with gene and exon information'
    )
    parser.add_argument(
        '--coverage_file',
        help='per base coverage data file'
    )

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    main()