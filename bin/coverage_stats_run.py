"""
Script to generate run level coverage statistics for a run of samples.
Takes as input an array of {sample_name}_exon_stats.tsv files output
from the coverage_stats_single.py, this is expected to be a sequencing
run of samples of the same capture to compare coverage across.

Run level coverage stats are calculated across all given samples, and
both an exon level and gene level output file are generated.

Jethro Rainford
jethro.rainford@addenbrookes.nhs.uk
201128
"""

import argparse
from math import sqrt
import multiprocessing
import os
from pathlib import Path
import re
import sys

import pandas as pd
import numpy as np


class runCoverage():

    def import_data(self, args):
        """
        Read in single stats coverage files to generate run stats from
        Args:
            - args (args): args passed from cmd line
        Returns:
            - sample_data (list): list of tuples with df and flagstat
                for each sample
        """
        # list to add dfs and flagstats dicts to, this will be a list of
        # tuples with (stats_df, flagstats) for each sample
        sample_data = []

        # read in exon stats files
        for file in args.files:
            with open(file, 'r') as f:
                flagstat = {}
                for line in f:
                    # read through header of stats file to get flagstats
                    if line.startswith('#'):
                        # flagstat stored in stats file as #key:value
                        stat = line.strip('#').replace("\n", "").split(':')
                        flagstat[stat[0]] = stat[1]
                    else:
                        # line after # is header, for adding data to df
                        header = line.split()
                        break
                # ensure flagstats were found
                assert bool(flagstat), (
                    'No flagstat found in exon stats header, this is required',
                    ' for generating run level stats. --flagstat is required ',
                    'when generating single stats. Exiting now.'
                )

                # add sample stats to df
                data = pd.read_csv(
                    f, sep="\t", names=header, low_memory=False, comment='#'
                )

            sample_data.append((data, flagstat))

        # check all dfs have same columns, check all against first
        assert all(
            [set(sample_data[0][0].columns) == set(df[0].columns) for df
                in sample_data]
        ), (
            'Mix of columns in input files, please use only exon'
            'stats files with the same threshold columns. Exiting.'
        )

        return sample_data


    def normalise_mean(self, sample_data, normal_reads):
        """
        Normalise mean of exons for sample against normalisation value
        and usable reads from flagstats.

        Args:
            - sample_data (tuple): tuple with sample stats df & flagstat
            - normal_reads (int): normalisation value (default: 1M)
        """
        if sample_data is None:
            # empty list passed
            return None

        sample_stats = sample_data[0][0]
        flagstat = sample_data[0][1]

        sample_stats["norm_mean"] = sample_stats["mean"].apply(
            lambda x: x * normal_reads / int(flagstat['usable_reads'])
        )

        return sample_stats


    def aggregate_exons(self, raw_stats):
        """
        Aggregates coverage stats for given exon coverage files and
        calculates standard deviation.

        Args:
            - stats_dfs (list): list of stats dfs
        Returns:
            - exon_stats (df): df of averaged run stats for given samples
        """
        # get list of genes and exons to calculate stats from
        exons = list(set(raw_stats[['gene', 'exon']].apply(tuple, axis=1)))
        exons.sort(key=lambda element: (element[0], element[1]))

        # empty df for run stats with same header
        exon_stats = raw_stats.iloc[0:0]
        exon_stats.insert(loc=8, column="std_dev", value="")

        for exon in exons:
            sample_exons = raw_stats.loc[
                (raw_stats["gene"] == exon[0]) & (raw_stats["exon"] == exon[1])
            ]
            sample_exons.index = range(len(sample_exons.index))

            row = sample_exons.iloc[0]
            num_samples = len(sample_exons.index)

            # get list of means and calculate standard deviation
            means = sample_exons["norm_mean"].tolist()
            std_dev = np.std(means)

            stats = {
                "chrom": row["chrom"],
                "exon_start": row["exon_start"],
                "exon_end": row["exon_end"],
                "gene": row["gene"],
                "tx": row["tx"],
                "exon": row["exon"],
                "exon_len": row["exon_len"],
                "min": sample_exons["min"].sum() / num_samples,
                "mean": (sample_exons["mean"].sum() / num_samples).round(2),
                "max": row["max"],
                "std_dev": std_dev,
                "norm_mean": row["norm_mean"]
            }

            # calculate mean for threshold columns
            threshold_cols = list(sample_exons.filter(regex='[0-9]+x', axis=1))

            for t in threshold_cols:
                stats[t] = (sample_exons[t].sum() / num_samples).round(2)

            exon_stats = exon_stats.append(stats, ignore_index=True)

        exon_stats.rename(columns={'mean': 'runMean'}, inplace=True)

        return exon_stats


    def gene_summary(self, exon_stats):
        """
        Gives averages per gene
        Args:
            - exon_stats (df): df of averaged exon run stats for given samples
        Returns:
            - gene_stats (df): df of averaged gene stats for run of samples
        """
        # get list of genes in data
        genes = exon_stats.gene.unique()

        # empty df for run stats with same header
        gene_stats = exon_stats.iloc[0:0]
        gene_stats = gene_stats.drop(
            ["exon_start", "exon_end", "exon", "exon_len"], axis=1
        )

        for gene in genes:
            exons = exon_stats.loc[exon_stats["gene"] == gene]
            exons.index = range(len(exons.index))
            row = exons.loc[0]

            # get fraction of gene of each exon
            exons["exon_frac"] = exons["exon_len"] / sum(exons["exon_len"])

            min = round(exons["min"].min(), 2)
            runMean = round(sum(
                [x * y for x, y in zip(exons["runMean"], exons["exon_frac"])]
            ), 2)
            norm_mean = round(sum(
                [x * y for x, y in zip(exons["norm_mean"], exons["exon_frac"])]
            ), 2)
            max = round(exons["max"].max(), 2)

            # calculate variance per exon to get std dev of gene
            exons["variance"] = exons["norm_mean"] * exons["std_dev"]
            gene_std_dev = sqrt(sum(exons["variance"]) / len(exons.index))
            gene_std_dev = round(gene_std_dev, 2)

            stats = {
                "gene": row["gene"],
                "tx": row["tx"],
                "chrom": row["chrom"],
                "min": min,
                "runMean": runMean,
                "std_dev": gene_std_dev,
                "max": max,
                "norm_mean": norm_mean
            }

            # get columns of threshold values, calculate avg of each
            thrshlds = exons.filter(regex="[0-9]+").columns.tolist()

            for thrshld in thrshlds:
                stats[thrshld] = round(
                    sum(exons[thrshld]) / len(exons.index), 2
                )

            # add gene totals to df
            gene_stats = gene_stats.append(stats, ignore_index=True)

        return gene_stats


    def write_outfile(self, exon_stats, gene_stats, outfile):
        """
        Write run level stats to tsv file
        Args:
            - exon_stats (df): run level aggregated stats of exon coverage
            - gene_stats (df): run level aggregated stats of gene coverage
            - outfile (str): output naming prefix for files

        Returns: None
        """
        # write report
        out_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "../output/"
        )
        outfile = os.path.join(out_dir, outfile)

        exon_stats.to_csv(
            outfile + "_run_exon_stats.tsv", sep="\t", index=False
        )
        gene_stats.to_csv(
            outfile + "_run_gene_stats.tsv", sep="\t", index=False
        )


    def parse_args(self):
        """
        Parse cmd line arguments
        Args: None
        Returns:
            - args (arguments): args passed from cmd line
        """
        parser = argparse.ArgumentParser(
            description='Generate run level coverage stats for multiple \
            samples.'
        )
        parser.add_argument(
            '--files', nargs='+',
            help='Exon stats files to generate run stats from.'
        )
        parser.add_argument(
            '--outfile', required=True,
            help='Output file name prefix', type=str
        )
        parser.add_argument(
            '--flank', required=False, default=5,
            help='flank included in bed file used to generate single sample\
                stats. Default: 5'
        )
        parser.add_argument(
            '--norm', required=False, default=1000000, type=int,
            help='Value to normalise against. Default: 1,000,000'
        )
        parser.add_argument(
            '--cores', nargs='?', default=None, type=int,
            help='Number of cores to utilise, for larger numbers of genes this\
            will drastically reduce run time. If not given will use maximum\
            available'
        )

        args = parser.parse_args()

        return args


def main():
    """
    Main to generate run level coverage stats
    """
    # turns off chained assignment warning - not req. as
    # intentionally writing back to df
    pd.options.mode.chained_assignment = None

    # initialise, pass args and import data
    run = runCoverage()
    args = run.parse_args()
    sample_data = run.import_data(args)

    # get total cores available for multiprocessing
    num_cores = multiprocessing.cpu_count()

    if args.cores is not None:
        # cores to use passed
        if int(args.cores) > num_cores:
            print(
                "Number cores given: {}, but only {} are available.",
                "Only using total cores available."
            )
        else:
            num_cores = args.cores

    # split samples equally for normalisation across seperate processes
    sample_data = np.array_split(np.array(sample_data), num_cores)

    # normalise means of all samples, returns list of dfs
    with multiprocessing.Pool(num_cores) as pool:
        print("Normalising all sample means")
        # if less samples than cores will result in empty lists, set to
        # None to handle in normalise_mean()
        sample_data = [x if len(x) > 0 else None for x in sample_data]

        raw_stats = pd.concat(pool.starmap(
            run.normalise_mean, map(lambda e: (e, args.norm), sample_data)
        ), ignore_index=True)

    raw_stats = raw_stats.sort_values(
        ["gene", "exon"], ascending=[True, True]
    )
    raw_stats.index = range(len(raw_stats.index))

    # split genes into equal arrays by num cores
    genes = raw_stats.gene.unique().tolist()
    genes_array = np.array_split(np.array(genes), num_cores)

    # split raw data by each gene-exon array to pass to each process
    split_exons = np.asanyarray(
        [raw_stats[raw_stats["gene"].isin(x)] for x in genes_array],
        dtype=object)

    with multiprocessing.Pool(num_cores) as pool:
        print("Aggregating exons and calculating std dev")
        exon_stats = pd.concat(pool.map(
            run.aggregate_exons, split_exons), ignore_index=True)

    # split genes into equal arrays by num cores
    genes = sorted(exon_stats.gene.unique().tolist())
    gene_array = np.array_split(np.array(genes), num_cores)

    # split exon stats df into seperate dfs by genes in each list
    split_exons = np.asanyarray(
        [exon_stats[exon_stats["gene"].isin(x)] for x in gene_array],
        dtype=object
    )

    with multiprocessing.Pool(num_cores) as pool:
        print("Calculating gene averages")
        gene_stats = pd.concat(pool.map(
            run.gene_summary, split_exons), ignore_index=True)

    run.write_outfile(exon_stats, gene_stats, args.outfile)


if __name__ == "__main__":

    main()
