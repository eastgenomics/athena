"""
Main script to control all running of Athena
"""
from datetime import datetime
from pathlib import PurePath
from time import time
from typing import Union
from uuid import uuid1

import pandas as pd

from bin import annotate, arguments, load, stats


class SubCommands():
    """
    Functions to individual sub parts of Athena (i.e. generate sample stats,
    run stats etc). Will either be called as single functions to do parts for
    development or multiple at a time.

    """

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


    @staticmethod
    def calculate_run_stats(all_exon_stats) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        _summary_

        Parameters
        ----------
        exon_stats : list(tuple(pd.DataFrame, pd.DataFrame))
            list of tuples of dataframes of per sample exon stats and
            hsmetrics files, one tuple per sample

        Returns
        -------
        pd.DataFrame
            dataframe of run level per exon coverage stats
        pd.DataFrame
            dataframe of run level per gene coverage stats
        """
        exon_stats = stats.Run().calculate_exon_stats(
            all_exon_stats=all_exon_stats)

        gene_stats = stats.Run().calculate_gene_stats(
            run_exon_stats=exon_stats
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
        start = time()
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
            f"\nFinished annotating bed file in {round((time() - start), 2)}s.\n"
            f"Output written to {args.output}_annotated_bed.tsv.gz"
        )

    elif args.sub == 'calculate_sample_stats':
        # sub command to generate single sample exon and gene stats
        # from a pre-annotated bed file
        annotated_bed = load.LoadData().read_annotated_bed(
            annotated_bed=args.annotated_bed
        )
        exon_stats, gene_stats = sub.calculate_sample_stats(
            data=annotated_bed,
            thresholds=args.thresholds
        )

        if not args.output:
            # set output prefix to be prefix of annotated bed
            args.output = args.annotated_bed.replace(
                ''.join(PurePath(args.annotated_bed).suffixes), ''
            ).replace('_annotated', '')

        if args.hsmetrics:
            metrics = load.LoadData().read_hsmetrics(args.hsmetrics)
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

        print(
            "\nFinished generating sample coverage stats, output written to:"
            f"\n\t{args.output}_exon_stats.tsv\n\t{args.output}_gene_stats.tsv"
        )

    elif args.sub == 'calculate_run_stats':
        # sub command to generate run level stats from a set of previously
        # calculated per sample exon stats files
        all_exon_stats = [
            load.LoadData().read_exon_stats(file) for file in args.exon_stats
        ]

        run_exon_stats, run_gene_stats = SubCommands().calculate_run_stats(
            all_exon_stats=all_exon_stats)
        
        # generate string of all samples used to generate run stats from
        # the write into header of output file
        now = datetime.now().strftime('%H:%M - %m/%d/%Y')
        samples = ','.join([
            PurePath(x).name.replace('_exon_stats.tsv', '')
            for x in args.exon_stats])
        
        samples = (
            f'#This file was generated at {now} from the samples:\n#{samples}\n'
        )

        exon_file = f"{args.run_prefix}_run_exon_stats.tsv"
        gene_file = f"{args.run_prefix}_run_gene_stats.tsv"
        
        with open(exon_file, 'w') as f1, open(gene_file, 'w') as f2:
            f1.write(samples)
            f2.write(samples)
        
        run_exon_stats.to_csv(exon_file, sep='\t', index=False, mode='a')
        run_gene_stats.to_csv(gene_file, sep='\t', index=False, mode='a')

        print(
            "\nFinished generating run coverage stats, output written to:\n\t"
            f"{args.run_prefix}_exon_stats.tsv\n\t{args.run_prefix}_gene_stats.tsv"
        )

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
