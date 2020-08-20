"""
Script to generate coverage statistics for a single sample.
Takes bed file annotated with gene, exon and mosdepth coverage,
example input file format in data/example_input_coverage.txt

Jethro Rainford 200721
"""

import argparse
import ast
import os
import sys
import pandas as pd


class singleCoverage():

    def import_data(self, args):
        """
        Import bed file with regions annotated with mosdepth coverage.
        Must have the expected columns in headers below and in the
        correct order.

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
            data = pd.read_csv(file, sep="\t", header=None, names=headers,
                               low_memory=False)

        if isinstance(args.thresholds, str):
            # thresholds passed as string
            # clean list of thresholds
            if len(args.thresholds) == 1:
                # list given with commas and no spaces
                thresholds = args.thresholds[0].split(",")

            # if list has commas
            thresholds = [i.strip(",") for i in args.thresholds]
        elif "," in args.thresholds:
            # given as list at cmd line, format correctly
            t = args.thresholds
            thresholds = [x.strip(",[]") for x in t]
        else:
            # using default list
            thresholds = args.thresholds

        return data, thresholds


    def cov_stats(self, data, thresholds):
        """
        Calculate coverage stats for sample

        Args:
            - data (df): dataframe of sample regions and coverage

        Returns:
            - cov_stats (df): df of coverage stats
        """
        print("Generating per base exon stats")

        header = [
            "chrom", "exon_start", "exon_end", "gene", "tx",
            "exon", "min", "mean", "max"
        ]

        # add thresholds to header list
        threshold_header = [str(i) + "x" for i in thresholds]
        header.extend(threshold_header)

        # get list of genes in data
        genes = data.gene.unique()

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

                # sort by coordinate in case of being out of order
                exon_cov = exon_cov.sort_values(by=["cov_start"])

                start = exon_cov.iloc[0]
                end = exon_cov.iloc[-1]

                # info for adding exon stats to output df
                row = exon_cov.iloc[0]

                if start["exon_start"] != start["cov_start"]:
                    # if cov_start is diff to tx start due to mosdepth
                    # binning, use tx start avoids wrongly estimating
                    # coverage by using wrong tx length
                    exon_cov.iloc[0, exon_cov.columns.get_loc("cov_start")] = int(start["exon_start"])

                if end["exon_end"] != end["cov_end"]:
                    # same as start
                    exon_cov.loc[
                        exon_cov.index[-1], "cov_end"] = int(end["exon_end"])
                
                # calculate summed coverage per bin
                exon_cov["cov_bin_len"] = exon_cov["cov_end"] -\
                    exon_cov["cov_start"]
                exon_cov["cov_sum"] = exon_cov["cov_bin_len"] * exon_cov["cov"]

                # calculate mean coverage from tx length and sum of coverage
                tx_len = int(end["exon_end"]) - int(start["exon_start"])
                mean_cov = round(exon_cov["cov_sum"].sum() / tx_len, 2)

                min_cov = exon_cov["cov"].min()
                max_cov = exon_cov["cov"].max()

                # get raw no. bases at each threshold
                raw_bases = {}
                for thrshld, header in zip(thresholds, threshold_header):
                    raw_bases[header] = exon_cov[
                        exon_cov["cov"] > int(thrshld)
                    ]["cov_bin_len"].sum()

                # calculate % bases at each threshold  from raw to 2 dp.
                pct_bases = {}
                for key, value in raw_bases.items():
                    pct_bases[key] = round(value / tx_len * 100, 2)

                stats = {
                    "chrom": row["chrom"], "exon_start": row["exon_start"],
                    "exon_end": row["exon_end"], "gene": gene, "tx": row["tx"],
                    "exon": row["exon"], "min": min_cov, "mean": mean_cov,
                    "max": max_cov
                }

                stats.update(pct_bases)

                cov_stats = cov_stats.append(stats, ignore_index=True)

        return cov_stats


    def summary_stats(self, cov_stats, thresholds):
        """
        Calculate per gene summary values

        Args:
            - cov_stats (df): df of per exon coverage stats

        Returns:
            - cov_summary (df): df of per gene coverage stats
        """
        print("Generating gene level summary stats")

        threshold_header = [str(i) + "x" for i in thresholds]

        # empty df for summary stats, use headers from stats table
        cov_summary = cov_stats.iloc[0:0]
        cov_summary = cov_summary.drop(
            ["chrom", "exon_start", "exon_end"], axis=1
        )

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
            gene_cov["exon_frac"] =\
                gene_cov["exon_len"] / sum(gene_cov["exon_len"])

            # calculate gene coverage values
            min = gene_cov["min"].min()
            mean = sum(
                [x * y for x, y in zip(
                    gene_cov["mean"], gene_cov["exon_frac"]
                )]
            )
            max = gene_cov["max"].max()

            # average coverage % at given thresholds
            thresholds = {}
            for t in threshold_header:
                thresholds[t] = sum(
                    [x * y for x, y in zip(gene_cov[t], gene_cov["exon_frac"])]
                )

            stats = {
                "gene": gene, "tx": row["tx"], "min": min, "mean": mean,
                "max": max
            }

            stats.update(thresholds)

            cov_summary = cov_summary.append(stats, ignore_index=True)

        # round calculated vals to 2 dp
        round_cols = ['mean'] + threshold_header
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

        cov_stats.to_csv(outfile + "_exon_stats.tsv", sep="\t", index=False)
        cov_summary.to_csv(outfile + "_gene_stats.tsv", sep="\t", index=False)


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
            '--file',
            help='annotated bed file on which to generate report from'
        )
        parser.add_argument(
            '--outfile', nargs='?', help='Output file name prefix, if not\
            given the input file name will be used as the name prefix.',
            type=str
        )
        parser.add_argument(
            '--thresholds', nargs='*',
            default=[10, 20, 30, 50, 100],
            help='List of threshold values seperated integers.'
        )

        args = parser.parse_args()

        if not args.outfile:
            # output file name not given
            args.outfile = args.file.rsplit(".")[0]

        return args


    def main(self):
        """
        Main to generate single coverage statistics output files
        """
        # turns off chained assignment warning - not req. as
        # intentionally writing back to df
        pd.options.mode.chained_assignment = None

        # parse arguments
        args = single.parse_args()

        # import data
        data, thresholds = single.import_data(args)

        # functions to generate coverage stats
        cov_stats = single.cov_stats(data, thresholds)
        cov_summary = single.summary_stats(cov_stats, thresholds)

        # write tables to output files
        if args.outfile:
            single.write_outfiles(cov_stats, cov_summary, args.outfile)


if __name__ == "__main__":

    single = singleCoverage()

    single.main()
