"""
Main script to control all running of Athena
"""
from pathlib import PurePath

import pandas as pd

from bin import annotate, arguments, load, stats


class SubCommands():

    @staticmethod
    def annotate_bed(panel_bed, exon_data, coverage_data, build) -> pd.DataFrame:
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

        Returns
        -------
        pd.DataFrame
            dataframe of annotated target / panel bed with exon and coverage data
        """
        annotated_bed = annotate.annotateBed().add_transcript_info(
            panel_bed=panel_bed,
            transcript_info=exon_data
        )

        annotated_bed = annotate.annotateBed().add_coverage(
            bed_w_transcript=annotated_bed,
            per_base_coverage=coverage_data,
            build=build
        )

        return annotated_bed


    @staticmethod
    def calculate_sample_stats(data, thresholds):
        """
        Calls functions to generate per exon and per gene stats

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of annotated bed file from SubCommands.annotated_bed()

        Returns:
        pd.DataFrame
            dataframe of per exon coverage stats
        pd.DataFrame
            dataframe of per gene coverage stats
        """
        exon_stats = stats.Sample().calculate_exon_stats_parallel(
            coverage_data=data,
            thresholds=thresholds
        )
        gene_stats = stats.Sample().calculate_gene_stats(
            coverage_data=data,
            exon_data=exon_stats,
            thresholds=thresholds
        )

        return exon_stats, gene_stats


def call_sub_command(args):
    """
    Calls given subcommand dependent on cmd line args

    Parameters
    ----------
    args : argparse.NameSpace
        argparse NameSpace object

    Returns
    -------

    """
    sub = SubCommands()

    if args.sub == 'annotate_bed_file':
        # sub command to annotate target bed file with transcript/exon and
        # per-base coverage data
        annotated_bed = sub.annotate_bed(
            panel_bed=args.target_bed,
            exon_data=args.exon_data,
            coverage_data=args.per_base_coverage,
            build=args.build
        )

        if not args.output:
            # set output prefix to be prefix of per base coverage file
            args.output = args.per_base_coverage.replace(
                ''.join(PurePath(args.per_base_coverage).suffixes), '')

        annotated_bed.to_csv(
            f"{args.output}_annotated_bed.tsv.gz",
            sep='\t', index=False
        )

        print(
            "\nFinished annotating bed file, output written to "
            f"{args.output}_annotated_bed.tsv.gz"
        )

    elif args.sub == 'calculate_sample_stats':
        # sub command to generate single sample exon and gene stats
        # from a pre-annotated bed file
        annotated_bed = load.loadData().read_annotated_bed(
            annotated_bed=args.annotated_bed
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

        if not args.output:
            # set output prefix to be prefix of annotated bed
            args.output = args.annotated_bed.replace(
                ''.join(PurePath(args.annotated_bed).suffixes), ''
            ).replace('_annotated', '')

        exon_stats.to_csv(
            f"{args.output}_exon_stats.tsv",
            sep='\t', index=False, mode=mode
        )
        gene_stats.to_csv(
            f"{args.output}_gene_stats.tsv", sep='\t', index=False)

        print(
            "\nFinished generating sample coverage stats, output written to:"
            f"\n\t{args.output}_exon_stats.tsv\n\t{args.output}_gene_stats.tsv"
        )

    elif args.sub == 'calculate_run_stats':
        pass

    elif args.sub == 'generate_report':
        pass

    else:
        pass


def main():
    """
    Main function to do all things Athena
    """
    args = arguments.parse_args()

    # calling a single sub command (i.e. annotate, sample stats, report...)
    if args.sub:
        call_sub_command(args)


if __name__ == "__main__":
    main()
