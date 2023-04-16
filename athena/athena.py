"""
Main script to generate coverage statistics and HTML report for given sample

Jethro Rainford
03/07/2021
"""
import argparse
import numpy as np
import os
import pandas as pd
from pathlib import Path

from bin import load, stats


class SubCommands():

    @staticmethod
    def annotate_bed():
        """
        Calls functions to annotate panel bed file with transcript info and coverage

        Parameters
        ----------

        Returns
        -------

        Outputs
        -------
        """
        pass

    
    @staticmethod
    def calculate_sample_stats(data, thresholds):
        """
        Calls functions to generate per exon and per gene stats

        Args:

        Returns:

        Outputs:
        """
        exon_stats = stats.stats().calculate_exon_stats_parallel(
            coverage_data=data,
            thresholds=thresholds
        )
        gene_stats = stats.stats().calculate_gene_stats(
            coverage_data=data,
            exon_data=exon_stats,
            thresholds=thresholds
        )

        return exon_stats, gene_stats


def call_sub_command(args):
    """
    Calls one or more of the subcommands dependent on cmd line args

    Parameters
    ----------
    args : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    sub = SubCommands()
    if args.sub == 'calculate_sample_stats':
        per_base_coverage = load.loadData().read_raw_coverage(
            args.annotated_bed
        )
        exon_stats, gene_stats = sub.calculate_sample_stats(
            data=per_base_coverage,
            thresholds=args.thresholds
        )

        if args.hsmetrics:
            metrics = load.loadData().read_hsmetrics(args.hsmetrics)
            metrics.to_csv(
                f"{args.output}_exon_stats.tsv",
                sep='\t', index=False, header=None
            )
            mode = 'a'
        else:
            mode = 'w'              

        exon_stats.to_csv(
            f"{args.output}_exon_stats.tsv",
            sep='\t', index=False, mode=mode
        )
        gene_stats.to_csv(
            f"{args.output}_gene_stats.tsv", sep='\t', index=False)



def parse_args():
    """
    Parse cmd line arguments

    Args:

    Returns:
    """
    parser = argparse.ArgumentParser(
        description='Generate coverage stats and HTML report.'
    )

    # generic inputs
    parser.add_argument(
        '--output', '-o', required=True,
        help='prefix for naming output files'
    )
    
    # sub commands for running seperate parts
    subparsers = parser.add_subparsers(
        title='sub_command', dest='sub',
        help='sub-commands for development and testing'
    )

    # parser to generate coverage stats
    stats_parser = subparsers.add_parser(
        'calculate_sample_stats',
        help='generate per gene and per exon coverage stats'
    )
    stats_parser.add_argument(
        '--annotated_bed',
        help='per base coverage bed file from annotated_bed'
    )
    stats_parser.add_argument(
        '--thresholds', nargs='*',
        default=[10, 20, 30, 50, 100],
        help='threshold values to calculate coverage for as comma\
            seperated integers (default: 10, 20, 30, 50, 100).'
    )
    stats_parser.add_argument(
        '--hsmetrics', nargs='?',
        help='Optional hsmetrics file, needed for generating run stats.'
    )
    stats_parser.add_argument(
        '--build', nargs='?',
        help='Optional text file with build number used for alignment.'
    )
    stats_parser.add_argument(
        '--outfile', nargs='?', help='Output file name prefix, if not\
        given the input file name will be used as the name prefix.',
        type=str
    )
    stats_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )

    args = parser.parse_args()

    return args


def main():
    """
    Main function to do all things Athena
    """
    args = parse_args()

    if args.sub:
        call_sub_command(args)


if __name__ == "__main__":
    main()
