"""
Script to generate coverage statistics for a single sample.
Takes bed file annotated with gene, exon and mosdepth coverage,
example input file format in data/example_input_coverage.txt

Jethro Rainford 200721
"""

import argparse
import os
import sys
import pandas as pd

class singleCoverage():

    def import_data(self, args):
        """
        Import bed file with regions annotated with mosdepth coverage.
        Must have the expected columns in headers below and in the correct order.

        Args:
            - args (args): passed arguments from arg parse

        Returns:
            - data (df): df of imported coverage file

        """
        with open(args.file) as file:
            headers = [
                "chrom", "exon_start", "exon_end",
                "gene", "tx", "exon", "cov_start",
                "cov_end", "cov"
            ]
            data = pd.read_csv(file, sep="\t", header=None, names=headers)

        return data


    def cov_stats(self, data):
        """
        Calculate coverage stats for sample
        
        Args:
            - data (df): dataframe of sample regions and coverage
        
        Returns:
            - cov_stats (df): df of coverage stats
        """
        # get list of genes in data
        genes = data.gene.unique()

        header = [
            "chrom", "exon_start", "exon_end", "gene", "tx",
            "exon", "min", "mean", "max",
            "10x", "20x", "30x", "50x", "100x"
        ]
        cov_stats = pd.DataFrame(columns=header)

        for gene in genes:
            
            # get coverage data for current gene
            gene_cov = data.loc[data["gene"] == gene]
            
            # get list of exons for gene
            exons = list(set(gene_cov["exon"].tolist()))

            for exon in exons:
                # calculate per exon coverage metrics
                
                # get coverage data for current exon
                exon_cov = gene_cov.loc[gene_cov["exon"] == exon]
                exon_cov.index = range(len(exon_cov.index))

                start = exon_cov.iloc[0]
                end = exon_cov.iloc[-1]

                # info for adding exon stats to output df
                row = exon_cov.iloc[0]

                if start["exon_start"] != start["cov_start"]:
                    # if cov_start is diff to tx start due to mosdepth binning, use tx start
                    # avoids wrongly estimating coverage by using wrong tx length
                    exon_cov.loc[0, "cov_start"] = int(start["exon_start"])

                if end["exon_end"] != end["cov_end"]:
                    # same as start
                    exon_cov.loc[
                        exon_cov.index[-1],
                        "cov_end"
                    ] = int(end["exon_end"])

                # calculate summed coverage per bin
                exon_cov["cov_bin_len"] = exon_cov["cov_end"] - exon_cov["cov_start"]
                exon_cov["cov_sum"] = exon_cov["cov_bin_len"] * exon_cov["cov"]

                # calculate mean coverage from tx length and sum of coverage
                tx_len = int(end["exon_end"]) - int(start["exon_start"])
                mean_cov = round(exon_cov["cov_sum"].sum() / tx_len, 2)

                min_cov = exon_cov["cov"].min()
                max_cov = exon_cov["cov"].max()

                # print("min: ", min_cov)
                # print("mean: ", round(mean_cov, 2))
                # print("max: ", max_cov)
                
                # get raw no. bases at each threshold
                bases_10x = exon_cov[exon_cov["cov"] > 10]["cov_bin_len"].sum()
                bases_20x = exon_cov[exon_cov["cov"] > 20]["cov_bin_len"].sum()
                bases_30x = exon_cov[exon_cov["cov"] > 30]["cov_bin_len"].sum()
                bases_50x = exon_cov[exon_cov["cov"] > 50]["cov_bin_len"].sum()
                bases_100x = exon_cov[exon_cov["cov"] > 100]["cov_bin_len"].sum()

                # calculate % bases at each threshold to 2 dp.
                pct_10x = round(bases_10x / tx_len * 100, 2)
                pct_20x = round(bases_20x / tx_len * 100, 2)
                pct_30x = round(bases_30x / tx_len * 100, 2)
                pct_50x = round(bases_50x / tx_len * 100, 2)
                pct_100x = round(bases_100x / tx_len * 100, 2)

                # print("percent at 10x: ", pct_10x)
                # print("percent at 20x: ", pct_20x)
                # print("percent at 30x: ", pct_30x)
                # print("percent at 50x: ", round(pct_50x, 2))
                # print("percent at 100x: ", round(pct_100x, 2))

                stats = {
                    "chrom": row["chrom"], "exon_start": row["exon_start"],
                    "exon_end": row["exon_end"], "gene": gene, "tx": row["tx"],
                    "exon": row["exon"], "min": min_cov, "mean": mean_cov, "max": max_cov, 
                    "10x": pct_10x, "20x": pct_20x, "30x": pct_30x, "50x": pct_50x,
                    "100x": pct_100x
                }

                cov_stats = cov_stats.append(stats, ignore_index=True)

        return cov_stats


    def summary_stats(self, cov_stats):
        """
        Calculate per gene summary values

        Args:
            - cov_stats (df): df of per exon coverage stats

        Returns:
            - cov_summary (df): df of per gene coverage stats
        """
        
        header = [
            "gene", "tx", "min", "mean", "max",
            "10x", "20x", "30x", "50x", "100x"
            ]
        
        # empty df for summary stats
        cov_summary = pd.DataFrame(columns=header)

        # calculate each exon len to get accurate stats
        cov_stats["exon_len"] = cov_stats["exon_end"] - cov_stats["exon_start"]

        # make list of genes
        genes = sorted(list(set(cov_stats["gene"].tolist())))

        for gene in genes:
        
            gene_cov = cov_stats.loc[cov_stats["gene"] == gene]
            gene_cov.index = range(len(gene_cov.index))

            # info for adding gene info to output df
            row = gene_cov.iloc[0]

            # calculate fraction of gene each exon covers
            # used to calculate each exon proportion of total gene metrics
            gene_cov["exon_frac"] = gene_cov["exon_len"] / sum(gene_cov["exon_len"])

            # calculate gene coverage values
            min = gene_cov["min"].min()
            mean = sum([x * y for x,y in zip(gene_cov["mean"], gene_cov["exon_frac"])])
            max = gene_cov["max"].max()

            x10 = sum([x * y for x,y in zip(gene_cov["10x"], gene_cov["exon_frac"])])
            x20 = sum([x * y for x,y in zip(gene_cov["20x"], gene_cov["exon_frac"])])
            x30 = sum([x * y for x,y in zip(gene_cov["30x"], gene_cov["exon_frac"])])
            x50 = sum([x * y for x,y in zip(gene_cov["50x"], gene_cov["exon_frac"])])
            x100 = sum([x * y for x,y in zip(gene_cov["100x"], gene_cov["exon_frac"])])

            stats = {
                    "gene": gene, "tx": row["tx"], "min": min, "mean": mean, 
                    "max": max, "10x": x10, "20x": x20, "30x": x30, 
                    "50x": x50, "100x": x100
                    }

            cov_summary = cov_summary.append(stats, ignore_index=True)

        # round calculated vals to 2 dp
        round_cols = ['mean', '10x', '20x', '30x', '50x', '100x']
        cov_summary[round_cols] = cov_summary[round_cols].round(2)

        return cov_summary


    def write_outfiles(self, cov_stats, cov_summary, outfile):
        """
        If --outfile arg given, writes coverage stats to file.
        
        Args:
            - cov_stats (df): df of generated coverage stats
            - args (args): includes name for output file
        
        Returns: None
        
        Outputs:
            - $outfile_exon_stats.tsv (file): tsv file of exon stats
            - $outfile_gene_stats.tsv (file): tsv file of gene stats
        """

        # write report
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        out_dir = os.path.join(bin_dir, "../output/")
        outfile = os.path.join(out_dir, outfile)

        cov_stats.to_csv(outfile+"_exon_stats.tsv", sep="\t", index=False)
        cov_summary.to_csv(outfile+"_gene_stats.tsv", sep="\t", index=False)


    def parse_args(self):
        """
        Parse cmd line arguments

        Args: None

        Returns:
            - args (arguments): args passed from cmd line
        """
        parser = argparse.ArgumentParser(
                    description='Generate coverage stats for a single sample.'
                    )               
        parser.add_argument(
        '--file', help='annotated bed file on which to generate report from')
        parser.add_argument(
        '--outfile', help='Output file name prefix', type=str)
        
        args = parser.parse_args()

        return args


    def main(self):
        """
        Main to generate single coverage statistics output files
        """
        # parse arguments
        args = single.parse_args()

        # import data
        data = single.import_data(args)

        # functions to generate coverage stats
        cov_stats = single.cov_stats(data)
        cov_summary = single.summary_stats(cov_stats)

        # write tables to output files
        if args.outfile:
            single.write_outfiles(cov_stats, cov_summary, args.outfile)


if __name__ == "__main__":

    single = singleCoverage()

    single.main()
        