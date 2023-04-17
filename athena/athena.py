"""
Main script to generate coverage statistics and HTML report for given sample

Jethro Rainford
03/07/2021
"""
import argparse
import pandas as pd

from .bin import annotate, load, stats


class SubCommands():

    @staticmethod
    def annotate_bed(panel_bed, exon_data, coverage_data, chunks) -> pd.DataFrame:
        """
        Calls functions to annotate panel bed file with transcript info and coverage

        Parameters
        ----------
        panel_bed : pd.DataFrame
            dataframe of target panel bed to annotate
        exon_data : pd.DataFrame
            dataframe of transcript and exon information
        coverage_data : pd.DataFrame
            dataframe of per base coverage information
        chunks : bool
            control if to annotate bed file in chunks to reduce peak
            memory usage

        Returns
        -------
        pd.DataFrame
            dataframe of annotated target / panel bed with exon and coverage data
        """
        annotated_bed = annotate.annotateBed().add_transcript_info(
            panel_bed=panel_bed,
            transcript_info_df=exon_data
        )

        annotated_bed = annotate.annotateBed().add_coverage(
            bed_w_transcript=annotated_bed,
            coverage_df=coverage_data,
            chunks=chunks
        )

        return annotated_bed


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
    args : argparse.NameSpace
        argparse NameSpace object

    Returns
    -------

    """
    sub = SubCommands()

    if args.sub == 'annotate_bed_file':
        # sub command to annotated target bed file with transcript/exon and
        # per-base coverage data
        per_base_coverage = load.loadData().read_coverage_data(
            coverage_file=args.per_base_coverage
        )
        target_bed = load.loadData().read_panel_bed(
            bed_file=args.target_bed
        )
        exon_data = load.loadData().read_transcript_info(
            transcript_file=args.exon_data
        )

        annotated_bed = sub.annotate_bed(
            panel_bed=target_bed,
            exon_data=exon_data,
            coverage_data=per_base_coverage,
            chunks=args.chunks
        )

        annotated_bed.to_csv(
            f"{args.output}_annotated_bed.tsv.gz",
            sep='\t', index=False
        )

    elif args.sub == 'calculate_sample_stats':
        # sub command to generate single sample exon and gene stats
        # from a pre-annotated bed file
        annotated_bed = load.loadData().read_raw_coverage(
            args.annotated_bed
        )
        exon_stats, gene_stats = sub.calculate_sample_stats(
            data=annotated_bed,
            thresholds=args.thresholds
        )

        if args.hsmetrics:
            metrics = load.loadData().read_hsmetrics(args.hsmetrics)
            metrics.to_csv(
                f"{args.output}_exon_stats.tsv",
                sep='\t', index=False, header=None
            )
            mode = 'a'  # set mode for writing df of exon data after hsmetrics
        else:
            mode = 'w'

        exon_stats.to_csv(
            f"{args.output}_exon_stats.tsv",
            sep='\t', index=False, mode=mode
        )
        gene_stats.to_csv(
            f"{args.output}_gene_stats.tsv", sep='\t', index=False)

    elif args.sub == 'calculate_run_stats':
        pass

    elif args.sub == 'generate_report':
        pass

    else:
        pass


def parse_args():
    """
    Parse cmd line arguments

    Returns
    -------
    argparse.NameSpace
        parsed command line arguments
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

    # parser to annotate bed file
    annotate_parser = subparsers.add_parser(
        'annotate_bed_file',
        help=(
            'annotate target bed file with transcript-exon information '
            'and per base coverage values'
        )
    )
    annotate_parser.add_argument(
        '--target_bed',
        help='bed file of target / panel to annotate'
    )
    annotate_parser.add_argument(
        '--exon data',
        help='tsv file of transcript-exon information'
    )
    annotate_parser.add_argument(
        '--per_base_coverage',
        help='per base coverage data (output from mosdepth)'
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
        '--hsmetrics', nargs='?', required=False,
        help=(
            'Optional hsmetrics file, needed for generating run stats. '
            'If given metrics will be written to first lines of exon stats '
            'file.'
        )
    )
    stats_parser.add_argument(
        '--build', nargs='?', required=False,
        help='Optional text file with build number used for alignment.'
    )

    return parser.parse_args()


def main():
    """
    Main function to do all things Athena
    """
    args = parse_args()

    if args.sub:
        call_sub_command(args)
    


if __name__ == "__main__":
    main()
