import argparse


def parse_args():
    """
    Parse cmd line arguments

    Returns
    -------
    argparse.Namespace
        parsed command line arguments
    """

    parser = argparse.ArgumentParser(
        description='Generate coverage stats and HTML report.'
    )

    # parse each set of args for sub commands and generic args
    parser = generic_arguments(parser)
    parser = subParsers(parser).parser

    return parser.parse_args()


def generic_arguments(parser):
    """
    Generic arguments not specific to any running mode
    
    Parameters
    ----------
    parser : argparse.Namespace
        Namespace object from argparse
    
    Returns
    -------
    argparse.Namespace
        Namespace object from argparse with generic args added
    """
    parser.add_argument(
        '--output', '-o', required=False,
        help=(
            'prefix for naming output files, if not given will default to '
            'using prefix of input sample file.'
        )
    )

    return parser


class subParsers():
    """
    Sub parsers for running individual sub commands
    """
    def __init__(self, parser) -> None:
        """
        Parameters
        ----------
        parser : argparse.Namespace
        Namespace object from argparse
        """
        self.parser = parser
        self.subparsers = self.parser.add_subparsers(
            title='sub_command', dest='sub',
            help='sub-commands for development and testing'
        )
        self.add_annotated_bed_file()
        self.add_calculate_sample_stats()
        self.add_calculate_run_stats()


    def add_annotated_bed_file(self):
        """
        Sub command for annotating panel bed file
        """
        annotate_parser = self.subparsers.add_parser(
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
            '--exon_data',
            help='tsv file of transcript-exon information'
        )
        annotate_parser.add_argument(
            '--per_base_coverage',
            help='per base coverage data (output from mosdepth)'
        )
        annotate_parser.add_argument(
            '--build', choices=[37, 38], type=int,
            help='Reference build of sample data'
        )


    def add_calculate_sample_stats(self):
        """
        Sub command for generating single sample exon and gene stats
        """
        stats_parser = self.subparsers.add_parser(
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
            help=(
                'threshold values to calculate coverage for as comma\
                seperated integers (default: 10, 20, 30, 50, 100).'
            )
        )
        stats_parser.add_argument(
            '--hsmetrics', nargs='?', required=False,
            help=(
                'Optional hsmetrics file, needed for generating run stats. '
                'If given metrics will be written to first lines of exon stats '
                'file.'
            )
        )


    def add_calculate_run_stats(self):
        """
        Sub command for generating normalised run level exon and gene stats
        """
        run_parser = self.subparsers.add_parser(
            'calculate_run_stats',
            help='generate run level gene and exon coverage stats'
        )
        run_parser.add_argument(
            '--exon_stats', nargs='+',
            help=(
                'exon stats for all samples from a run, these must have the '
                'hsmetrics given when generating per sample stats'
            )
        )
        run_parser.add_argument(
            '--run_prefix', required=False,
            help=(
                'Prefix for naming run level coverage stats files (e.g. the '
                'run ID), if not given will generate a random ID for naming'
            )
        )
