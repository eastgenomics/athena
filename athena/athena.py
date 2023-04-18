"""
Main script to control all running of Athena
"""
from pathlib import PurePath
from typing import Union

import pandas as pd

from bin import annotate, arguments, load, stats


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

        Parameters
        ----------

        Returns
        -------

        """
        exon_stats = stats.sample().calculate_exon_stats_parallel(
            coverage_data=data,
            thresholds=thresholds
        )
        gene_stats = stats.sample().calculate_gene_stats(
            coverage_data=data,
            exon_data=exon_stats,
            thresholds=thresholds
        )

        return exon_stats, gene_stats


    @staticmethod
    def calculate_run_stats(exon_stats) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        _summary_

        Parameters
        ----------
        exon_stats : list(pd.DataFrame)
            list of dataframes of per sample exon stats

        Returns
        -------
        Union[pd.DataFrame, pd.DataFrame]
            _description_
        """
        print(exon_stats)

        


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

        if not args.output:
            # set output prefix to be prefix of per base coverage file
            args.output = args.per_base_coverage.replace(
                ''.join(PurePath(args.per_base_coverage).suffixes), '')

        annotated_bed.to_csv(
            f"{args.output}_annotated_bed.tsv.gz",
            sep='\t', index=False
        )

        print(
            "Finished annotating bed file, output written to "
            f"{args.output}_annotated_bed.tsv.gz"
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
        # sub command to generate run level stats from a set of previously
        # calculated per sample exon stats files
        all_exon_stats = [
            load.loadData().read_exon_stats(file) for file in args.exon_stats
        ]
        SubCommands().calculate_run_stats(exon_stats=all_exon_stats)

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
