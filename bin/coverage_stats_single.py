"""
Script to generate coverage statistics for a single sample.
Takes bed file annotated with gene, exon and mosdepth coverage,
example input file format in data/example_input_coverage.txt

Jethro Rainford 200721
"""

import argparse
from functools import partial
import os
import re
import sys
import math
import multiprocessing
import numpy as np
import pandas as pd

from pathlib import Path


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

            dtypes = {
                "chrom": str, "exon_start": int, "exon_end": int,
                "gene": str, "tx": str, "exon": str, "cov_start": int,
                "cov_end": int, "cov": int
            }

            data = pd.read_csv(
                file, sep="\t", header=None, names=headers, dtype=dtypes
            )
            # strip chr from chrom in cases of diff. formatted bed
            data["chrom"] = data["chrom"].apply(
                lambda x: str(x).replace("chr", "")
            )

        # clean list of thresholds
        if len(args.thresholds) == 1:
            # list given with commas and no spaces
            args.thresholds = args.thresholds[0].strip("[]").split(",")

        if all(isinstance(x, str) for x in args.thresholds):
            # if passed will be strings, default will be int
            # strip other items from list
            thresholds = [
                int(re.sub(r'\W', '', i)) for i in args.thresholds if i
            ]
            # remove duplicates if given
            thresholds = sorted(list(set(thresholds)))
        else:
            # using default list
            thresholds = args.thresholds

        if args.build:
            # build no. file given
            with open(args.build) as b:
                build = b.readline().strip("\n")
                if "37" in build:
                    build = "GRCh37 ({})".format(build)
                if "38" in build:
                    build = "GRCh38 ({})".format(build)
        else:
            build = ""

        flagstat = {}

        if args.flagstat:
            with open(args.flagstat) as f:
                # flagstat file passed
                # read each flagstat file and add metrics to dict
                lines = [line.rstrip('\n') for line in f]
                for line in lines:
                    if re.search('total', line):
                        flagstat['total_reads'] = re.sub(
                            r'^(\d+) .*reads.*', r'\1', line
                        )
                    elif re.search(r'duplicates', line):
                        flagstat['dups_reads'] = re.sub(
                            r'^(\d+) .*duplicates', r'\1', line
                        )
                    elif re.search(r'mapped \(', line):
                        flagstat['mapped_reads'] = re.sub(
                            r'^(\d+) .*mapped.*', r'\1', line
                        )
                    elif re.search(r'properly paired', line):
                        flagstat['properly paired'] = re.sub(
                            r'^(\d+) .*properly paired .*', r'\1', line
                        )
                    elif re.search(r'singletons', line):
                        flagstat['singletons'] = re.sub(
                            r'^(\d+) .*singletons .*', r'\1', line
                        )
                flagstat['usable_reads'] = int(
                    flagstat['mapped_reads']
                ) - int(flagstat['dups_reads'])

        return data, thresholds, flagstat, build


    def cov_stats(self, data, thresholds):
        """
        Calculate coverage stats for sample

        Args:
            - data (df): dataframe of sample regions and coverage

        Returns:
            - cov_stats (df): df of coverage stats
        """
        header = [
            "chrom", "exon_start", "exon_end", "gene", "tx",
            "exon", "min", "mean", "max"
        ]

        # add thresholds to header list
        threshold_header = [str(i) + "x" for i in thresholds]
        header.extend(threshold_header)

        # get list of genes in data
        genes = sorted(data.gene.unique().tolist())

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

                # get unique list of exon start & end, should always be
                # just one, if not exon has been split => error
                coords = exon_cov[["exon_start", "exon_end"]].to_records(
                    index=False
                )
                coords = list(np.unique(coords))

                # if more than one region for exon in bed, exit as will
                # be incorrectly calculated
                assert len(coords) == 1, (
                    "More than one region is present in the bed file for "
                    "exon {} of {}: {}.\n".format(exon, gene, coords),
                    "Currently each exon MUST be in one pair of start / "
                    "end coordinates else coverage values will be "
                    "incorrect for those regions. Exiting now."
                )

                start = exon_cov.iloc[0]
                end = exon_cov.iloc[-1]

                # info for adding exon stats to output df
                row = exon_cov.iloc[0]

                if start["exon_start"] != start["cov_start"]:
                    # if cov_start is diff to tx start due to mosdepth
                    # binning, use tx start avoids wrongly estimating
                    # coverage by using wrong tx length
                    exon_cov.iloc[0, exon_cov.columns.get_loc(
                        "cov_start")] = int(start["exon_start"])

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

                assert tx_len > 0, "transcript length for {} ({}) appears to\
                    be 0. Exon start: {}. Exon end: {}".format(
                    start["gene"], start["tx"],
                    start["exon_start"], end["exon_end"]
                )

                mean_cov = round(exon_cov["cov_sum"].sum() / tx_len, 2)

                min_cov = exon_cov["cov"].min()
                max_cov = exon_cov["cov"].max()

                # get raw no. bases at each threshold
                raw_bases = {}
                for thrshld, header in zip(thresholds, threshold_header):
                    raw_bases[header] = exon_cov[
                        exon_cov["cov"] >= int(thrshld)
                    ]["cov_bin_len"].sum()

                # calculate % bases at each threshold
                pct_bases = {}
                for key, value in raw_bases.items():
                    raw_value = value / tx_len * 100
                    pct_bases[key] = raw_value

                stats = {
                    "chrom": row["chrom"], "exon_start": row["exon_start"],
                    "exon_end": row["exon_end"], "gene": gene, "tx": row["tx"],
                    "exon": row["exon"], "min": min_cov, "mean": mean_cov,
                    "max": max_cov
                }

                stats.update(pct_bases)

                cov_stats = cov_stats.append(stats, ignore_index=True)

        # calculate each exon len to get accurate stats
        cov_stats["exon_len"] = cov_stats["exon_end"] - cov_stats["exon_start"]

        return cov_stats


    def summary_stats(self, cov_stats, thresholds):
        """
        Calculate per gene summary values

        Args:
            - cov_stats (df): df of per exon coverage stats

        Returns:
            - cov_summary (df): df of per gene coverage stats
        """
        threshold_header = [str(i) + "x" for i in thresholds]

        # empty df for summary stats, uses header from stats table
        cov_summary = cov_stats.iloc[0:0]
        cov_summary = cov_summary.drop(
            ["chrom", "exon_start", "exon_end", "exon_len"], axis=1
        )

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
            mean = round(sum(
                [x * y for x, y in zip(gene_cov["mean"], gene_cov["exon_frac"])]
            ), 12)
            max = gene_cov["max"].max()

            # average coverage % at given thresholds, round to 12 dp to
            # account for accuracy of float values when summing and
            # length of human genome
            thresholds = {}
            for t in threshold_header:
                thresholds[t] = float(round(sum(
                    [x * y for x, y in zip(gene_cov[t], gene_cov["exon_frac"])]
                ), 12))

            stats = {
                "gene": gene, "tx": row["tx"], "min": min, "mean": mean,
                "max": max
            }

            stats.update(thresholds)

            cov_summary = cov_summary.append(stats, ignore_index=True)

        return cov_summary


    def write_outfiles(self, cov_stats, cov_summary, outfile, flagstat, build):
        """
        Writes both exon and gene level coverage stats to file.

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
        outfile = os.path.join(out_dir, Path(outfile).stem)

        exon_stats = outfile + "_exon_stats.tsv"
        gene_stats = outfile + "_gene_stats.tsv"

        # check if file already exists, overwrite empty if true
        if os.path.exists(exon_stats):
            open(exon_stats, 'w').close()

        if os.path.exists(gene_stats):
            open(gene_stats, 'w').close()

        if flagstat:
            # if flagstat file given
            flags = ""
            for key, value in flagstat.items():
                flags += "#" + key + " : " + str(value) + "\n"

            with open(outfile + "_exon_stats.tsv", 'w+') as file:
                file.write(flags)

            with open(outfile + "_gene_stats.tsv", 'w+') as file:
                file.write(flags)

        if build:
            # if build file given
            with open(exon_stats, 'a+') as file:
                file.write("# build: " + build + "\n")

            with open(gene_stats, 'a+') as file:
                file.write("# build: " + build + "\n")

        # write stats files
        with open(exon_stats, 'a+') as file:
            cov_stats.to_csv(exon_stats, sep="\t", mode='a', index=False)

        with open(gene_stats, 'a+') as file:
            cov_summary.to_csv(gene_stats, sep="\t", mode='a', index=False)


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
            help='annotated bed file on which to generate stats from'
        )
        parser.add_argument(
            '--flagstat', nargs='?',
            help='Optional flagstat file, needed for generating run stats.'
        )
        parser.add_argument(
            '--build', nargs='?',
            help='Optional text file with build number used for alignment.'
        )
        parser.add_argument(
            '--outfile', nargs='?', help='Output file name prefix, if not\
            given the input file name will be used as the name prefix.',
            type=str
        )
        parser.add_argument(
            '--thresholds', nargs='*',
            default=[10, 20, 30, 50, 100],
            help='threshold values to calculate coverage for as comma\
                seperated integers (default: 10, 20, 30, 50, 100).'
        )
        parser.add_argument(
            '--cores', nargs='?', default=None,
            help='Number of cores to utilise, for larger numbers of genes this\
            will drastically reduce run time. If not given will use maximum\
            available'
        )

        args = parser.parse_args()

        if not args.outfile:
            # output file name not given
            args.outfile = Path(args.file).stem
            # remove extension if present (from annotate_bed.sh)
            args.outfile = args.outfile.replace("_annotated", "")
            args.outfile = args.outfile.replace("_markdup", "")

        return args


def main():
    """
    Main to generate single coverage statistics output files
    """
    # turns off chained assignment warning - not req. as
    # intentionally writing back to df
    pd.options.mode.chained_assignment = None

    single = singleCoverage()

    # parse arguments
    args = single.parse_args()

    # import data
    data, thresholds, flagstat, build = single.import_data(args)

    # get total cores available
    num_cores = multiprocessing.cpu_count()

    if args.cores is not None:
        # cores to use passed
        if int(args.cores) > num_cores:
            print(
                "Number cores given: {}, but only {} are available.\
                Only using total cores available.".format(
                    args.cores, num_cores
                )
            )
        else:
            num_cores = int(args.cores)

    # get list of genes in data
    genes = sorted(data.gene.unique().tolist())

    # split gene list equally for seperate processes
    gene_array = np.array_split(np.array(genes), num_cores)

    # split df into seperate dfs by genes in each list
    split_dfs = np.asanyarray(
        [data[data["gene"].isin(x)] for x in gene_array], dtype=object
    )

    with multiprocessing.Pool(num_cores) as pool:
        # use a pool to spawn multiple processes
        # uses number of cores defined and splits processing of df
        # slices, add each to pool with threshold values and
        # concatenates together when finished
        print("Generating per base exon stats")
        try:
            # map threshold arg to function since same is passed to all, then
            # use imap to pass each df to worker to generate stats
            stats_w_arg = partial(single.cov_stats, thresholds=thresholds)
            pool_output = pool.imap_unordered(stats_w_arg, split_dfs)
        except AssertionError:
            # more than one region for each exon found => exit
            pool.close()
            pool.terminate()
        else:
            cov_stats = pd.concat(pool_output, ignore_index=True)

            # imap_unordered() returns everything out of order (funnily enough)
            # sort by gene and exon to be nicely formatted
            cov_stats.exon = cov_stats.exon.astype(int)
            cov_stats = cov_stats.sort_values(['gene', 'exon'])

    # split up output coverage stats df for multiprocessing
    split_stats_dfs = np.asanyarray(
        [cov_stats[cov_stats["gene"].isin(x)] for x in gene_array],
        dtype=object
    )

    with multiprocessing.Pool(num_cores) as pool:
        print("Generating gene level summary stats")

        cov_summary = pd.concat(
            pool.starmap(
                single.summary_stats, map(
                    lambda e: (e, thresholds), split_stats_dfs
                )), ignore_index=True)

    cov_summary = cov_summary.sort_values(['gene'])
    cov_summary = cov_summary.drop(columns=["exon"])

    # write tables to output files
    single.write_outfiles(
        cov_stats, cov_summary, args.outfile, flagstat, build
    )


if __name__ == "__main__":

    main()
