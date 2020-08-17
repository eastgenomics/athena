"""
Script to generate run level coverage statistics for a run of samples.

Takes as input an array of {sample_name}_exon_stats.tsv files output 
from the coverage_stats_single.py, this is expected to be a sequencing 
run of samples to compare coverage across.

Run level coverage stats are calculated across all given samples, and 
both an exon level and gene level output file are generated.

Jethro Rainford 200814
jethro.rainford@addenbrookes.nhs.uk
"""

import argparse
import numpy as np
import os
import pandas as pd
import sys

from math import sqrt


class runCoverage():

    def import_data(self, args):
        """
        Read in single stats coverage files to generate run stats from

        Args:
            - args (args): args passed from cmd line

        Returns:
            - stats_dfs (list): list of all stats dfs 
        """
        stat_dfs = []

        for file in args.files:
            with open(file):
                data = pd.read_csv(file, sep="\t", header=0, low_memory=False)
                stat_dfs.append(data)
        
        # check all dfs are same type of file (i.e exon or gene stats)
        if not all(
            [set(stat_dfs[0].columns) == set(df.columns) for df in stat_dfs]):
            print('Mix of columns in input files, please use only exon stats\
                    files with the same threshold columns. Exiting.')
            sys.exit()

        return stat_dfs


    def aggregate_exons(self, stat_dfs):
        """
        Aggregates coverage stats for given exon coverage files

        Args:
            - stats_dfs (list): list of all stats dfs 
        
        Returns:
            - exon_stats (df): df of averaged run stats for given samples
        """

        # combine all dfs, sort by gene and exon
        raw_stats = pd.concat(stat_dfs)
        raw_stats = raw_stats.sort_values(
                            ["gene", "exon"], ascending=[True, True])
        raw_stats.index = range(len(raw_stats.index))
    
        # get list of genes and exons to calculate stats from
        exons = list(set(raw_stats[['gene', 'exon']].apply(tuple, axis=1)))
        exons.sort(key=lambda element: (element[0], element[1]))

        # empty df for run stats with same header
        exon_stats = raw_stats.iloc[0:0]
        exon_stats.insert(loc=8, column="std_dev", value="")

        for exon in exons:

            sample_exons = raw_stats.loc[(raw_stats["gene"] == exon[0]) &
                        (raw_stats["exon"] == exon[1])]
            sample_exons.index = range(len(sample_exons.index))

            row = sample_exons.iloc[0]
            num_samples = len(sample_exons.index)

            # get list of means and calculate standard deviation
            means = sample_exons["mean"].tolist()
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
                "std_dev": std_dev
                }

            # calculate mean for threshold columns
            threshold_cols = list(sample_exons.filter(regex='[0-9]+x', axis=1))

            for t in threshold_cols:
                stats[t] = (sample_exons[t].sum() / num_samples).round(2)

            exon_stats = exon_stats.append(stats, ignore_index = True)

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
                                    ["exon_start", "exon_end", 
                                    "exon", "exon_len"], axis = 1
                                    )

        for gene in genes:
            exons = exon_stats.loc[exon_stats["gene"] == gene]
            exons.index = range(len(exons.index))
            row = exons.loc[0]

            # get fraction of gene of each exon
            exons["exon_frac"] = exons["exon_len"] / sum(exons["exon_len"])

            min = round(exons["min"].min(), 2)
            mean = round(sum([x * y for x,y in zip(exons["mean"],
                                                    exons["exon_frac"])]), 2)
            max = round(exons["max"].max(), 2)

            # calculate variance per exon to get std dev of gene
            exons["variance"] = exons["mean"] * exons["std_dev"]
            gene_std_dev = sqrt(sum(exons["variance"]) / len(exons.index))
            gene_std_dev = round(gene_std_dev, 2)

            stats = {
                    "gene": row["gene"],
                    "tx": row["tx"],
                    "chrom": row["chrom"],
                    "min": min,
                    "mean": mean,
                    "std_dev": gene_std_dev,
                    "max": max
                }

            # get columns of threshold values, calculate avg of each
            thrshlds = exons.filter(regex = "[0-9]+").columns.tolist()

            for thrshld in thrshlds:
                stats[thrshld] = round(
                                    sum(exons[thrshld]) / len(exons.index), 2
                                    )
            
            # add gene totals to df
            gene_stats = gene_stats.append(stats, ignore_index = True)
        
        return gene_stats


    def write_outfile(self, exon_stats, gene_stats, outfile):
        """
        Write run level stats to tsv file

        Args:
            - stats_dfs (list): list of all stats dfs 
        
        Returns: None
        """
        # write report
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        out_dir = os.path.join(bin_dir, "../output/")
        outfile = os.path.join(out_dir, outfile)

        exon_stats.to_csv(outfile + "_exon_stats.tsv", sep="\t", index=False)
        gene_stats.to_csv(outfile + "_gene_stats.tsv", sep="\t", index=False)


    def parse_args(self):
        """
        Parse cmd line arguments

        Args: None

        Returns:
            - args (arguments): args passed from cmd line
        """
        parser = argparse.ArgumentParser(
                    description='Generate run level coverage stats for\
                                multiple samples.'
                    )               
        parser.add_argument('--files', nargs='+',
                    help='exon stats files on which to generate report from'
                    )
        parser.add_argument('--outfile', required=True,
                help='Output file name prefix', type=str
                    )
                    
        args = parser.parse_args()
        
        return args
    

    def main(self):
        """
        Main to generate run level coverage stats
        """

        # turns off chained assignment warning - not req. as intentionally
        # writing back to df
        pd.options.mode.chained_assignment = None

        args = run.parse_args()

        stat_dfs = run.import_data(args)

        exon_stats = run.aggregate_exons(stat_dfs)

        gene_stats = run.gene_summary(exon_stats)

        run.write_outfile(exon_stats, gene_stats, args.outfile)


if __name__ == "__main__":

    run = runCoverage()

    run.main()
