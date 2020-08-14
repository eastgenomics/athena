"""
Script to generate run level coverage statistics for a run of samples.

Takes as input an array of files output from the 
coverage_stats_single.py, this is expected to be a sequencing run of
samples to compare coverage across.

Jethro Rainford 200814
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
        if not all([set(stat_dfs[0].columns) == set(df.columns) for df in stat_dfs]):
            print('Mix of stats files passed, please use all exon OR gene stats. Exiting.')
            sys.exit()

        return stat_dfs
    
    
    def standard_dev(self, means):
        """
        Calculate standard deviation from given sample means

        Args:
            - means (list): list of mean values

        Returns:
            - std_dev (float): std dev from given means
        """

        sqr_sum = 0.0

        for mean in means:
            sqr_sum += float(mean)**2

        mean_sum = sum(means)
        std_dev = sqrt(
            (len(means) * sqr_sum - mean_sum * mean_sum) /
            (len(means) * (len(means) -1 ))
            )

        std_dev = round(std_dev, 2)

        return std_dev


    def aggregate_exons(self, stat_dfs):
        """
        Aggregates coverage stats for given exon coverage files

        Args:
            - stats_dfs (list): list of all stats dfs 
        
        Returns:
            - run_stats (df): df of averaged run stats for given samples
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
        run_stats = raw_stats.iloc[0:0]
        runs_stats = run_stats.insert(loc=8, column="std_dev", value="")

        for exon in exons:

            exon_stats = raw_stats.loc[(raw_stats["gene"] == exon[0]) & (raw_stats["exon"] == exon[1])]
            exon_stats.index = range(len(exon_stats.index))

            row = exon_stats.iloc[0]
            num_samples = len(exon_stats.index)

            # get list of means and calculate standard deviation
            means = exon_stats["mean"].tolist()
            std_dev = self.standard_dev(means)

            stats = {
                    "chrom": row["chrom"],
                    "exon_start": row["exon_start"],
                    "exon_end": row["exon_end"],
                    "gene": row["gene"], 
                    "tx": row["tx"],
                    "exon": row["exon"],
                    "exon_len": row["exon_len"],
                    "min": exon_stats["min"].sum() / num_samples,
                    "mean": (exon_stats["mean"].sum() / num_samples).round(2),
                    "std_dev": std_dev,
                    "max": exon_stats["max"].sum() / num_samples,
                    "10x": (exon_stats["10x"].sum() / num_samples).round(2),
                    "20x": (exon_stats["20x"].sum() / num_samples).round(2),
                    "30x": (exon_stats["30x"].sum() / num_samples).round(2),
                    "50x": (exon_stats["50x"].sum() / num_samples).round(2),
                    "100x": (exon_stats["100x"].sum() / num_samples).round(2)
                    }

            run_stats = run_stats.append(stats, ignore_index = True)
            
        return run_stats


    def write_outfile(self, run_stats, outfile):
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

        run_stats.to_csv(outfile+"_run_stats.tsv", sep="\t", index=False)


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
                    help='annotated bed file on which to generate report from'
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

        args = run.parse_args()

        stat_dfs = run.import_data(args)

        run_stats = run.aggregate_exons(stat_dfs)

        run.write_outfile(run_stats, args.outfile)


if __name__ == "__main__":

    run = runCoverage()

    run.main()