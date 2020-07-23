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


def import_data(args):
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


def cov_stats(data, output):
    """
    Calculate coverage stats for sample
    
    Args:
        - data (df): dataframe of sample regions and coverage
    
    Returns:
        - cov_stats (df): df of coverage stats
        - sub_20x (df): df of sub optimal regions
    """
    # get list of genes in data
    genes = data.gene.unique()

    header = [
        "chrom", "exon_start", "exon_end", "gene", "tx",
        "exon", "min", "mean", "max",
        "10x", "20x", "30x", "50x", "100x"
    ]
    cov_stats = pd.DataFrame(columns=header)
    sub_20x = pd.DataFrame(columns=header)

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

            if int(pct_20x) < 100:
                sub_20x = sub_20x.append(stats, ignore_index=True)

    return cov_stats, sub_20x


    def write_outfile(cov_stats, args):
        """
        If --outfile arg given, writes coverage stats to file.
        
        Args:
            - cov_stats (df): df of generated coverage stats
            - args (args): includes name for output file
        
        Returns: None
        
        Outputs:
            - outfile (file): tab seperated output file of coverage stats
        """
        cov_stats.to_csv(outfile, sep="\t", index=False)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate coverage stats for a single sample.'
    )
    parser.add_argument('--file', help='annotated bed file on which to generate report from')
    parser.add_argument('--outfile', help='Output file name')
    parser.add_argument('--plots', help='', nargs='?')
    args = parser.parse_args()

    # functions to generate coverage stats
    data = import_data(args)
    cov_stats, sub_20x = cov_stats(data, args.output)
    
    if args.outfile:
        write_outfile(cov_stats, args)
        

