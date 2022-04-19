"""
Script to generate single sample coverage report.
Takes single sample coverage stats files as input, along with the raw
coverage input file and an optional "low" coverage threshold (default 20).

Jethro Rainford 200722
"""

import argparse
import base64
import math
import matplotlib
# use agg instead of tkinter for pyplot backend
matplotlib.use('agg')
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import pandas as pd
import pandasql as pdsql
import plotly.graph_objs as go
import pybedtools as bedtools
import sys

from datetime import datetime
from io import BytesIO
from pathlib import Path
from plotly.subplots import make_subplots
from string import Template
from matplotlib.ticker import ScalarFormatter

import load_data


class generatePlots():
    """Functions to generate required plots"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.threshold_val = int(threshold.strip("x"))


    def img2str(self, plt):
        """
        Converts png of plot to HTML formatted string
        Args:
            - plt (image): png of plot
        Returns:
            - img (str): HTML formatted string of plot
        """
        buffer = BytesIO()
        plt.savefig(buffer, format='png', dpi=65, transparent=True)
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        graphic = base64.b64encode(image_png)
        data_uri = graphic.decode('utf-8')
        img_tag = (
            f"<img src=data:image/png;base64,{data_uri} style='max-width: "
            "100%; max-height: auto; object-fit: contain; ' />"
        )

        return img_tag


    def low_exon_plot(self, low_raw_cov):
        """
        Generate array of low exon plot values to pass into report

        Args:
            - low_raw_cov (df): df of raw coverage for exons with low
                                coverage

        Returns:
            - low_exon_plots (str): list of plot values in div tags
        """
        if len(low_raw_cov.index) == 0:
            # empty df passed, likely from difference in total cores and plots
            return

        # get list of tuples of transcripts and exons to define plots
        transcripts = low_raw_cov.drop_duplicates(
            ["tx", "exon"])[["tx", "exon"]].values.tolist()
        transcripts = [tuple(exon) for exon in transcripts]

        # sort list of genes/exons by gene and exon
        transcripts = sorted(
            transcripts, key=lambda element: (element[0], element[1])
        )

        low_raw_cov["exon_len"] =\
            low_raw_cov["exon_end"] - low_raw_cov["exon_start"]

        low_raw_cov["relative_position"] = low_raw_cov["exon_end"] - round(((
            low_raw_cov["cov_end"] + low_raw_cov["cov_start"]) / 2
        ))

        low_exon_plots = []  # array to add string data of plots to

        for tx in transcripts:
            # get rows for current gene and exon
            exon_cov = low_raw_cov.loc[(
                low_raw_cov["tx"] == tx[0]
            ) & (
                low_raw_cov["exon"] == tx[1]
            )]

            exon_cov = exon_cov.sort_values(by='cov_start', ascending=True)
            start = exon_cov.iloc[0]
            end = exon_cov.iloc[-1]

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

            # create empty df for unbinned data with same columns
            exon_cov_unbinned = exon_cov[0:0]

            for i, row in exon_cov.iterrows():
                for pos in range(row["cov_start"], row["cov_end"] + 1):
                    # unbin each row, set start & end to same value for each
                    # use +1 since range is non inclusive of final value
                    pos_row = row
                    pos_row["cov_start"] = pos
                    pos_row["cov_end"] = pos
                    exon_cov_unbinned = exon_cov_unbinned.append(
                        pos_row, ignore_index=True
                    )

            if sum(exon_cov_unbinned["cov"]) == 0:
                continue

            # build div str of plot data to pass to template
            x_vals = str(exon_cov_unbinned['cov_start'].tolist()).strip('[]')
            y_vals = str(exon_cov_unbinned['cov'].tolist()).strip('[]')
            title = f"{tx[0]} exon {tx[1]}"

            tx_data = (
                f"""'<div class="sub_plot">{title},{x_vals},{y_vals}</div>'"""
            )

            low_exon_plots.append(tx_data)

        low_exon_plots = ','.join(low_exon_plots)

        return low_exon_plots


    def all_gene_plots(self, raw_coverage):
        """
        Generate full plots for each gene / transcript

        Args:
            - raw_coverage (file): from args; bp coverage file used as
                                    input for coverage_stats_single.py
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - all-plots (str): str of lists of all plots with gene symbol
        """
        all_plots = ""

        if len(raw_coverage.index) == 0:
            # passed empty df, most likely because there were less genes
            # than processes => empty df passed with multiprocess
            return ""

        # get unique list of genes
        transcripts = raw_coverage.drop_duplicates(["tx"])["tx"].values.tolist()

        for tx in transcripts:
            # get coverage data for current gene
            tx_cov = raw_coverage.loc[(raw_coverage["tx"] == tx)]
            # get list of exons
            exons = tx_cov.drop_duplicates(["exon"])["exon"].tolist()

            # no. plot columns = no. of exons
            column_no = len(exons)

            gene = tx_cov.loc[0]['gene']

            # make subplot grid size of no. of exons, height variable
            # splits large genes to several rows and maintains height
            height = math.ceil(len(exons) / 30) * 4.5
            fig = plt.figure(figsize=(30, height))

            # generate grid with space for each exon
            # splits genes with >25 exons to multiple rows
            rows = math.ceil(len(exons) / 30)
            if column_no > 30:
                column_no = 30

            grid = fig.add_gridspec(rows, column_no, wspace=0)
            axs = grid.subplots(sharey=True)

            if column_no == 1:
                # handle single exon genes, axs needs turning into np
                # array to flatten
                axs = np.array([axs])

            axs = axs.flatten()

            fig.suptitle(f"{gene} ({tx})", fontweight="bold", fontsize=14)
            count = 0

            for exon in exons:
                # get coverage data for current exon
                exon_cov = raw_coverage.loc[(
                    raw_coverage["tx"] == tx
                ) & (
                    raw_coverage["exon"] == exon
                )]

                exon_cov = exon_cov.reset_index(drop=True)

                # sort and check coordinates are correct
                exon_cov = exon_cov.sort_values(by='cov_start', ascending=True)

                start = exon_cov.iloc[0]
                end = exon_cov.iloc[-1]

                if start["exon_start"] != start["cov_start"]:
                    # if cov_start is diff to tx start due to mosdepth
                    # binning, use tx start avoids wrongly estimating
                    # coverage by using wrong tx length
                    exon_cov.iloc[
                        0, exon_cov.columns.get_loc("cov_start")
                    ] = int(start["exon_start"])

                if end["exon_end"] != end["cov_end"]:
                    # same as start
                    exon_cov.loc[exon_cov.index[-1], "cov_end"] = int(
                        end["exon_end"]
                    )

                # check if coverage column empty
                if (exon_cov['cov'] == 0).all():
                    # no coverage, generate empty plot with just
                    # threshold line
                    axs[count].plot(
                        [0, 100], [self.threshold_val, self.threshold_val],
                        color='red', linestyle='-', linewidth=2, rasterized=True
                    )
                else:
                    axs[count].plot(
                        exon_cov["cov_start"].tolist(),
                        exon_cov["cov"].tolist()
                    )

                    # threshold line
                    axs[count].plot(
                        [exon_cov["exon_start"], exon_cov["exon_end"]],
                        [self.threshold_val, self.threshold_val], color='red',
                        linestyle='-', linewidth=1, rasterized=True
                    )

                # add labels
                xlab = str(
                    exon_cov["exon_end"].iloc[0] -
                    exon_cov["exon_start"].iloc[0]
                ) + " bp"

                if len(exons) > 20:
                    # drop bp to new line for better spacing
                    xlab = xlab.replace("bp", "\nbp")

                axs[count].title.set_text(exon)
                axs[count].set_xlabel(xlab, fontsize=13)

                count += 1

            # remove y ticks & label for all but first plot of lines
            for i in range(column_no * rows):
                if i in [x * column_no for x in range(rows)]:
                    # first plot of line, keep ticks and labels
                    axs[i].tick_params(axis='y', labelsize=12)
                    continue
                else:
                    axs[i].yaxis.set_ticks_position('none')

            # strip x axis ticks and labels
            plt.setp(plt.gcf().get_axes(), xticks=[])

            # adjust yaxis limits
            ymax = max(tx_cov["cov"].tolist()) + 10
            plt.ylim(bottom=0, top=ymax)

            # remove outer white margins
            fig.tight_layout(h_pad=1.4)

            # convert plot png to html string and append to one string
            img = self.img2str(plt)

            # add img to str list with gene symbol for filtering in table
            # expects to be a string of lists to write in report
            img_str = f'["{gene} ({tx})", "{img}" ], '
            all_plots += img_str

            plt.close(fig)

        return all_plots


    def summary_gene_plot(self, cov_summary):
        """
        Generate summary plot of all genes against threshold value

        Args:
            - cov_summary (df): df of gene coverage values
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - summary_plot (fig): plot of all genes
        """
        print("Generating summary plot")

        summary_data = cov_summary.copy()

        # define colours based on values
        summary_data["colours"] = 'green'
        summary_data.loc[
            summary_data[self.threshold] < 100, 'colours'] = 'orange'
        summary_data.loc[summary_data[self.threshold] < 90, 'colours'] = 'red'

        summary_data = summary_data.sort_values(
            by=[self.threshold], ascending=False
        )
        summary_plot, axs = plt.subplots(figsize=(25, 7.5))

        genes100pct = None

        if len(summary_data.index) > 100:
            # split off some of 100% covered genes to limit size of plot
            if len(summary_data[summary_data[self.threshold] < 100]) > 100:
                # over 100 sub threshold genes, remove all 100% genes
                genes100pct = len(
                    summary_data[summary_data[self.threshold] == 100]
                )
                summary_data = summary_data[summary_data[self.threshold] < 100]
            else:
                # split off bottom 100 genes, plot includes some 100% covered
                genes100pct = len(summary_data.iloc[:-100])
                summary_data = summary_data.iloc[-100:]

        # generate the plot
        plt.bar(
            summary_data["tx"],
            [int(x) for x in summary_data[self.threshold]],
            color=summary_data.colours
        )

        if genes100pct is not None:
            # more than 100 genes, add title inc. 100% covered not shown
            genes100pct = str(genes100pct)
            axs.set_title(
                r"$\bf{" + genes100pct + "}$" + " genes covered 100% at " +
                r"$\bf{" + self.threshold + "}$" +
                " were omitted from the plot due to the panel size", loc='left'
            )

        # threshold lines
        plt.axhline(y=99, linestyle='--', color="#565656", alpha=0.6)
        plt.axhline(y=95, linestyle='--', color="#565656", alpha=0.6)

        plt.text(1.005, 0.94, '99%', transform=axs.transAxes)
        plt.text(1.005, 0.91, '95%', transform=axs.transAxes)

        # plot formatting
        axs.tick_params(labelsize=8, length=0)
        plt.xticks(rotation=55, color="#565656")

        # adjust whole plot margins
        axs.autoscale_view(scaley=True)
        if len(cov_summary.index) < 10:
            # add wider margins for low no. genes to stop v. wide bars
            margin = 1 / len(cov_summary.index) / len(cov_summary.index)
            axs.margins(x=margin)
        else:
            axs.margins(x=0.01)

        # add legend
        green = mpatches.Patch(color='green', label='100%')
        orange = mpatches.Patch(color='orange', label='90-99.99%')
        red = mpatches.Patch(color='red', label='<90%')

        plt.legend(
            handles=[green, orange, red], loc='upper center',
            bbox_to_anchor=(0.5, -0.18),
            fancybox=True, shadow=True, ncol=12, fontsize=14
        )

        vals = np.arange(0, 110, 10).tolist()
        plt.yticks(vals, vals)

        if len(summary_data.index) > 250:
            # x axis label overlap on lots of genes => skip every third
            # shouldn't be a huge amount more than this to stop overlap
            axs.set_xticks(axs.get_xticks()[::3])
            axs.tick_params(axis='both', which='major', labelsize=10)
            plt.figtext(
                0.5, 0.1,
                "Some gene labels are not shown due to high number of genes",
                ha="center", fontsize=12
            )
        elif len(summary_data.index) > 125:
            # x axis label overlap on lots of genes => skip every second
            axs.set_xticks(axs.get_xticks()[::2])
            axs.tick_params(axis='both', which='major', labelsize=10)
            plt.figtext(
                0.5, 0.1,
                "Some gene labels are not shown due to high number of genes",
                ha="center", fontsize=12
            )
        else:
            axs.tick_params(axis='both', which='major', labelsize=10)

        plt.xlabel("")
        plt.ylabel(f"% coverage ({self.threshold})", fontsize=11)

        axs.yaxis.grid(linewidth=0.5, color="grey", linestyle="-.")
        plt.box(False)
        axs.set_axisbelow(True)
        plt.tight_layout()

        # convert image to html string to insert in report
        summary_plot = self.img2str(plt)

        return summary_plot

    def coverage_per_chromosome_plot(
        self, per_base_coverage: pd.DataFrame, nrows: int = 6,
        ncols: int = 4, sharey: bool = True
    ) -> str:
        """
        Produce plots of coverage per chromosome, given a data-frame
        of coverage from the per-base.bed.gz data output from mosdepth.
        A base64-encoded string of the plots is returned.

        Args:
            per_base_coverage: per_base.bed.gz data frame
            nrows: number of subplot rows
            ncols: number of subplot columns
            sharey: if true, all plots share the same y-axis limits.

        Returns:
           base64-encoded string representation of coverage plots
        """

        chr_index = [str(i) for i in range(1, 23)] + ["X"] + ["Y"]

        grouped_coverage = per_base_coverage.groupby("chrom")

        fig, axs = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=(25, 25),
                sharey=sharey,
                constrained_layout=True
            )

        fig.set_constrained_layout_pads(
                w_pad=0.6,
                h_pad=0.2,
                hspace=0.0,
                wspace=0.0
            )

        for i, chrom_name in enumerate(chr_index):

            # choose first row, iterate over all columns, go to next row etc.
            row_index = i // ncols # row index will be 0, 0, 0, 0, 1, 1, 1, 1 etc.
            col_index = i % ncols  # col index will be 0, 1, 2, 3, 0, 1, 2, 3 etc.

            # select plot
            ax = axs[row_index][col_index]

            # plot data
            ax.scatter(data=grouped_coverage.get_group(chrom_name),
                       x="start", y="cov", s=1)

            # set plot text parameters
            ax.set_title(f"chr{chrom_name}", fontsize=24, fontstyle='italic')
            ax.tick_params(axis='both', labelsize=20)

            # adjust size and format of the scientific notation label
            ax.xaxis.offsetText.set_fontsize(18)
            ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))

        # add a common y-axis label
        fig.text(x=0, y=0.5, s="depth", va="center", rotation="vertical", fontsize=34)

        return self.img2str(fig)


class styleTables():
    """Functions for styling tables for displaying in report"""

    def __init__(
        self, cov_stats: pd.DataFrame, cov_summary: pd.DataFrame,
        threshold: str, threshold_cols: list, vals: list
    ) -> None:
        self.cov_stats = cov_stats
        self.cov_summary = cov_summary
        self.threshold = threshold
        self.threshold_val = int(threshold.strip("x"))
        self.threshold_cols = threshold_cols
        self.vals = vals
        self.column_names = {
            "gene": "Gene",
            "tx": "Transcript",
            "chrom": "Chr",
            "exon": "Exon",
            "exon_len": "Length",
            "exon_start": "Start",
            "exon_end": "End",
            "min": "Min",
            "mean": "Mean",
            "max": "Max"
        }


    def style_sub_threshold(self):
        """
        Styling of sub threshold stats df for displaying in report

        Returns:
            - sub_threshold_stats (list): list of sub threshold coverage stats
            - low_exon_columns (list): list of column headers for report
            - gene_issues (int): total number of genes under threshold
            - exon_issues (int): total number of exons under threshold
        """
        column = [
            "gene", "tx", "chrom", "exon", "exon_len", "exon_start",
            "exon_end", "min", "mean", "max"
        ]
        column.extend(self.threshold_cols)

        dtypes = {
            'gene': str, 'tx': str, 'chrom': str, 'exon': int, 'exon_len': int,
            'exon_start': int, 'exon_end': int, 'min': int, 'mean': float,
            'max': int
        }

        for col in self.threshold_cols:
            # add dtype for threshold columns
            dtypes[col] = float

        sub_threshold = pd.DataFrame(columns=column)
        sub_threshold.astype(dtype=dtypes)

        # get all exons with <100% coverage at threshold
        for i, row in self.cov_stats.iterrows():
            if int(row[self.threshold]) < 100:
                sub_threshold = sub_threshold.append(row, ignore_index=True)

        if not sub_threshold.empty:
            # some low covered regions identified
            sub_threshold = sub_threshold.astype(dtypes)

            sub_threshold_stats = pd.pivot_table(sub_threshold, index=[
                "gene", "tx", "chrom", "exon",
                "exon_len", "exon_start", "exon_end"
            ], values=self.vals)

            # reset index to fix formatting
            sub_threshold_stats = sub_threshold_stats.reindex(self.vals, axis=1)
            sub_threshold_stats.reset_index(inplace=True)
            sub_threshold_stats.index = np.arange(
                1, len(sub_threshold_stats.index) + 1
            )

            gene_issues = len(list(set(sub_threshold_stats["gene"].tolist())))
            exon_issues = len(sub_threshold_stats["exon"])
        else:
            # if no low regions set to empty df with appropriate columns
            print("No low coverage regions, generating empty table")
            sub_threshold_stats = pd.DataFrame(columns=column)
            gene_issues = 0
            exon_issues = 0

        # rename columns to display properly
        sub_threshold_stats = sub_threshold_stats.rename(
            columns=self.column_names
        )

        # add index as column to have numbered rows in report
        sub_threshold_stats.insert(0, ' ', sub_threshold_stats.index)

        # threshold_cols -> list of strings, add mean to loop over for rounding
        round_cols = ['Mean'] + self.threshold_cols

        # limit to 2dp using math.floor, use of round() with
        # 2dp may lead to inaccuracy such as 99.99 => 100.00
        for col in round_cols:
            sub_threshold_stats[col] = sub_threshold_stats[col].map(
                lambda col: math.floor(col * 100) / 100
            )

        # generate list of dicts with column headers for styling
        low_exon_columns = []

        for col in sub_threshold_stats.columns:
            low_exon_columns.append({'title': col})

        # convert df to list of lists to pass to report
        sub_threshold_stats = sub_threshold_stats.values.tolist()

        return sub_threshold_stats, low_exon_columns, gene_issues, exon_issues


    def style_total_stats(self):
        """
        Styling of full gene-exon stats table for displaying in report
        Args:
            - cov_stats (df): df of exon stats
            - threshold_cols (list): list of threshold columns
            - vals (list): list of min, mean and max strs
        Returns:
            - total_stats (str): HTML formatted string of cov_stats df
        """
        # do some excel level formatting to make table more readable
        total_stats = pd.pivot_table(
            self.cov_stats,
            index=["gene", "tx", "chrom", "exon", "exon_len",
                   "exon_start", "exon_end"],
            values=self.vals
        )

        # reset index to fix formatting, set beginning to 1
        total_stats = total_stats.reindex(self.vals, axis=1)
        total_stats.reset_index(inplace=True)
        total_stats.index = np.arange(1, len(total_stats.index) + 1)
        total_stats.insert(0, 'index', total_stats.index)

        total_stats = total_stats.rename(columns=self.column_names)

        # limit to 2dp using math.floor, use of round() with
        # 2dp may lead to inaccuracy such as 99.99 => 100.00
        round_cols = ['Mean'] + self.threshold_cols

        for col in round_cols:
            total_stats[col] = total_stats[col].map(
                lambda col: math.floor(col * 100) / 100
            )

        # turn gene stats table into list of lists
        total_stats = total_stats.values.tolist()

        return total_stats


    def style_cov_summary(self):
        """
        Add styling to per gene coverage summary table
        Args:
            - cov_summary (df): df of gene coverage stats
            - threshold_cols (list): list of threshold values
        Returns:
            - gene_stats (list): HTML formatted str of gene summary df
            - total_genes (int): total number of genes
        """
        # rename columns for displaying in report
        gene_stats = self.cov_summary.copy()
        gene_stats = gene_stats.rename(columns=self.column_names)

        # get values to display in report
        total_genes = len(gene_stats["Gene"].tolist())

        # limit to 2dp using math.floor, use of round() with
        # 2dp may lead to inaccuracy such as 99.99 => 100.00
        round_cols = ['Mean'] + self.threshold_cols

        for col in round_cols:
            gene_stats[col] = gene_stats[col].map(
                lambda col: math.floor(col * 100) / 100
            )

        # reset index to start at 1
        gene_stats.index = np.arange(1, len(gene_stats.index) + 1)
        gene_stats.insert(0, 'index', gene_stats.index)

        # turn gene stats table into list of lists
        gene_stats = gene_stats.values.tolist()

        return gene_stats, total_genes


    def style_snps_cov(self, snps_cov):
        """
        Add styling to tables of SNPs covered above / beneath threshold
        Args:
            - snps_cov (df): df of snps above / below threshold
        Returns:
            - snps_cov (list): list of snps to render in report
            - total_snps (int): total number of snps in df
        """
        if not snps_cov.empty:
            # set uuid appropriate for passed df
            if all(
                x > self.threshold_val for x in snps_cov["Coverage"].tolist()
                ):
                uuid = "var_high_cov"
            else:
                uuid = "var_low_cov"

            snps_cov.index = np.arange(1, len(snps_cov.index) + 1)
            total_snps = len(snps_cov.index)

            # reset index to start at 1
            snps_cov.index = np.arange(1, len(snps_cov.index) + 1)
            snps_cov.insert(0, 'index', snps_cov.index)

            # turn gene stats table into list of lists
            snps_cov = snps_cov.values.tolist()
        else:
            snps_cov = []
            total_snps = 0

        return snps_cov, total_snps

    @staticmethod
    def style_snps_no_cov(snps_no_cov):
        """
        Add styling to table of snps that span exon boundaries => have
        coverage values
        Args:
            - snps_no_cov (df): df of snps with no coverage values
        Returns:
            - snps_no_cov (list): list of snps for rendering in report
            - snps_out_panel (int): total number snps with no cov
        """
        # if variants from vcf found that span exon boundaries
        if not snps_no_cov.empty:
            # get number of variants to display in report
            snps_out_panel = len(snps_no_cov.index)

            # reset index to start at 1
            snps_no_cov.index = np.arange(1, len(snps_no_cov.index) + 1)
            snps_no_cov.insert(0, 'index', snps_no_cov.index)

            # turn gene stats table into list of lists
            snps_no_cov = snps_no_cov.values.tolist()
        else:
            snps_no_cov = []
            snps_out_panel = 0

        return snps_no_cov, snps_out_panel


class calculateValues():
    """Functions to calculate values and write summary for report"""

    def __init__(self, threshold):
        self.threshold = threshold


    def panel_coverage(self, cov_stats):
        """
        Calculates mean coverage of all panel regions at given threshold,
        normalised against length of each gene

        Args:
            - cov_stats (df): df of coverage stats for each exon
            - threshold (int): threshold cut off for low coverage

        Returns:
            - panel_pct_coverage (str): % coverage of panel as str
        """
        print("Calculating panel average coverage")

        gene_stats = pd.DataFrame(
            columns=["gene", "gene_len", "coverage"])

        # make list of genes
        genes = sorted(list(set(cov_stats["gene"].tolist())))

        for gene in genes:
            # for each gene, calculate length and average % at threshold
            gene_cov = cov_stats.loc[cov_stats["gene"] == gene]

            length = sum(gene_cov["exon_len"])
            coverage = sum(
                gene_cov[self.threshold] * gene_cov["exon_len"] / length)

            gene_stats = gene_stats.append({
                "gene": gene,
                "gene_len": length,
                "coverage": coverage
            }, ignore_index=True)

        # calculate % panel coverage
        panel_coverage = sum(
            gene_stats["coverage"] * gene_stats["gene_len"] / sum(
                gene_stats["gene_len"]
            )
        )

        # round to 12 dp to account for limit of accuracy of float &
        # length of human genome
        panel_coverage = round(panel_coverage, 12)

        panel_pct_coverage = str(math.floor(panel_coverage * 100) / 100)

        return panel_pct_coverage


    def snp_coverage(self, snp_vcfs, raw_coverage):
        """
        Produces tables of coverage for variants inside of capture
        regions, and larger structural variants spanning region
        boundaries.

        Args:
            - snp_vcfs (str): list of vcf files used for SNP analysis
            - raw_coverage (df): raw bp coverage for each exon
            - threshold (int): threshold value passed from parse args

        Returns:
            - snps_low_cov (df): variants with lower coverage than threshold
            - snps_high_cov (df): variants with higher coverage than threshold
            - snps_no_cov (df): variants that span exon boundaries (i.e SVs)
        """
        print("Calculating coverage of given SNPs")

        bedFile = raw_coverage[
            ["chrom", "exon_start", "exon_end"]].drop_duplicates()
        coverageFile = raw_coverage[
            ["chrom", "cov_start", "cov_end", "cov"]].drop_duplicates()

        # turn dfs into BedTools objects
        bed = bedtools.BedTool.from_dataframe(bedFile)
        cov = bedtools.BedTool.from_dataframe(coverageFile)

        # empty df to add all SNP info to
        snp_df = pd.DataFrame(columns=[
            'VCF', 'chrom', 'pos', 'id', 'ref', 'alt', 'info'
        ])

        for vcf in snp_vcfs:
            # read vcf into BedTools object
            v = bedtools.BedTool(vcf)

            # get vcf name to add to table, req. for multiple VCFS and
            # recording variant source VCF
            name = Path(vcf).stem.split("_")[0]

            # use bedtools intersect to get SNPs in capture region
            snps = bed.intersect(v, wb=True)

            for row in snps:
                # get data from returned BedTools object, add to df
                snp_data = str(row).split()
                snp_df = snp_df.append({
                    'VCF': name, 'chrom': snp_data[3],
                    'pos': snp_data[4], 'ref': snp_data[6],
                    'alt': snp_data[7], 'info': snp_data[10]
                }, ignore_index=True)

        snp_df = snp_df[
            ['VCF', 'chrom', 'pos', 'ref', 'alt', 'info']].drop_duplicates()

        # reset index
        raw_coverage = raw_coverage.reset_index(drop=True)

        # use pandasql to intersect SNPs against coverage df to find the
        # coverage at each SNP position
        coverage_sql = """
            SELECT snp_df.VCF, snp_df.chrom, snp_df.pos, snp_df.ref,
            snp_df.alt, snp_df.info, raw_coverage.gene, raw_coverage.exon,
            raw_coverage.cov_start, raw_coverage.cov_end, raw_coverage.cov
            FROM snp_df
            LEFT JOIN raw_coverage on snp_df.CHROM=raw_coverage.chrom
            WHERE snp_df.POS > raw_coverage.cov_start AND
            snp_df.POS <= raw_coverage.cov_end
            """

        snp_cov = pdsql.sqldf(coverage_sql, locals())

        # get SNPs that won't have coverage data but do intersect panel
        # regions (i.e. large deletions that span a region)
        snps_no_cov = snp_df.merge(snp_cov, how='outer', indicator=True).loc[
            lambda x: x['_merge'] == 'left_only']

        snps_no_cov = snps_no_cov[[
            "VCF", "chrom", "pos", "ref", "alt", "info"
        ]].reset_index(drop=True)

        # get required columns for SNP tables
        snps_cov = snp_cov[
            ["VCF", "gene", "exon", "chrom", "pos", "ref", "alt", "cov"]
        ].drop_duplicates(subset=[
            "VCF", "chrom", "pos", "ref", "alt"]).reset_index(drop=True)

        # rename columns for displaying in report
        snps_cov.columns = ["VCF", "Gene", "Exon", "Chromosome", "Position",
                            "Ref", "Alt", "Coverage"]

        snps_no_cov.columns = [
            "VCF", "Chromosome", "Position", "Ref", "Alt", "Info"
        ]

        # remove <> from DELs to stop being interpreted as HTML tags
        snps_no_cov["Alt"] = snps_no_cov["Alt"].str.strip("<>")

        snps_cov["Coverage"] = snps_cov["Coverage"].astype(int)

        # sort no_cov table by chrom & pos, as pos is str first define
        # order to sort by
        order = [str(x) for x in range(0, 23)]
        order.extend(["X", "Y", "MT"])
        snps_no_cov["Chromosome"] = pd.Categorical(
            snps_no_cov["Chromosome"], order
        )

        snps_cov = snps_cov.sort_values(by=["Gene", "Exon", "Position"])
        snps_no_cov = snps_no_cov.sort_values(by=["Chromosome", "Position"])

        # split SNPs by coverage against threshold
        threshold_int = int(self.threshold.strip("x"))
        snps_low_cov = snps_cov.loc[snps_cov["Coverage"] < threshold_int]
        snps_high_cov = snps_cov.loc[snps_cov["Coverage"] >= threshold_int]

        return snps_low_cov, snps_high_cov, snps_no_cov

    @staticmethod
    def calculate_snp_vals(snps_covered, snps_not_covered, snps_out_panel):
        """
        Calculate % values for SNP totals
        Args:
            - snps_covered (int): total number snps covered at threshold
            - snps_not_covered (int): total number snps not covered at
                threshold
            - snps_out_panel (int): total number snps spanning exon
                boundaries
        Returns:
            - total_snps (int): sum of all snps
            - snps_pct_covered (float): % value of snps_covered
            - snps_pct_not_covered (float): % value of snps_not_covered
            - snps_pct_out_panel (float): % value of snps_out_panel
        """
        total_snps = str(snps_covered + snps_not_covered + snps_out_panel)

        # calculate % SNPs covered vs. not, limit to 2dp with math.floor
        if snps_covered != 0:
            snps_pct_covered = int(snps_covered) / int(total_snps) * 100
            snps_pct_covered = math.floor(snps_pct_covered * 100) / 100
        else:
            snps_pct_covered = 0

        if snps_not_covered != 0:
            snps_pct_not_covered = int(
                snps_not_covered) / int(total_snps) * 100
            snps_pct_not_covered = math.floor(snps_pct_not_covered * 100) / 100
        else:
            snps_pct_not_covered = 0

        if snps_out_panel != 0:
            snps_pct_out_panel = int(
                snps_out_panel) / int(total_snps) * 100
            snps_pct_out_panel = math.floor(snps_pct_out_panel * 100) / 100
        else:
            snps_pct_out_panel = 0

        return total_snps, snps_pct_covered, snps_pct_not_covered,\
            snps_pct_out_panel


class generateReport():
    """Functions to combine variables and generate report"""

    def __init__(self, threshold):
        self.threshold = threshold

    def write_summary(self, cov_summary, threshold, panel_pct_coverage):
        """
        Write summary paragraph with sequencing details and list of
        genes / transcripts used in panel.

        Args:
            - cov_summary (df): df of gene coverage values
            - threshold (int): defined threshold level (default: 20)
            - panel_pct_coverage (str): % coverage of panel as str
        Returns:
            - summary_text (str): summary text with req. HTML markup
        """
        pct_cov = str(math.floor(float(panel_pct_coverage)))

        # summary text paragraph with div styling
        summary_text = """
        <li>Clinical report summary:</li>
        <div style="background-color:aliceblue; margin-top: 15px;
        border-radius: 15px; padding-left:25px; overflow-y: auto;
        max-height:500px;"><div id="summary_text" style="font-size: 14px;
        padding-bottom: 15px; padding-top:10px">"""

        sub90 = ""

        for idx, gene in cov_summary.iterrows():
            # build string of each gene, trascript and coverage at
            # threshold to display in summary
            summary = "{} ({}); ".format(gene["gene"], gene["tx"])
            summary_text += summary

            if gene[self.threshold] < 90:
                # build string of genes with <90% coverage at threshold
                sub90 += (
                    f'{gene["gene"]} ({gene["tx"]}) '
                    f'{math.floor(gene[self.threshold] * 100)/100.0}%; '
                )

        summary_text = summary_text.strip(" ;") + "."

        if not sub90:
            # all genes >90% at threshold
            sub90 = "<b>None</b>"
        else:
            sub90 = sub90.strip(" ;") + "."

        summary_text += (
            f"<br></br><b>Genes with coverage at {self.threshold} "
            f"less than 90%:</b> {sub90}"
        )

        summary_text += (
            f"<br></br>{pct_cov} % of this panel was sequenced to a depth "
            f"of {self.threshold} or greater.<br></div></div>"
        )

        return summary_text

    def generate_report(self, cov_stats, cov_summary, snps_low_cov,
                        snps_high_cov, snps_no_cov, fig, all_plots,
                        coverage_per_chromosome_fig,
                        summary_plot, html_template, args, build, panel, vcfs,
                        panel_pct_coverage, bootstrap, version, summary_text
                        ):
        """
        Generate single sample report from coverage stats

        Args:
            - cov_stats (df): df of coverage stats for each exon
            - cov_summary (df): df of gene level coverage
            - snps_low_cov (df): SNPs with lower coverage than threshold
            - snps_high_cov (df): SNPs with higher coverage than threshold
            - snps_no_cov (df): variants that span exon boundaries (i.e SVs)
            - fig (figure): plots of low coverage regions
            - all-plots (figure): grid of all full gene- exon plots
            - coverage_per_chromosome_fig: coverage-per-chromosome figure
            - summary_plot (figure): gene summary plot - % at threshold
            - html_template (str): string of HTML template
            - args (args): passed cmd line arguments
            - build (str): build number used for alignment
            - panel (str): panes(s) / gene(s) included in report
            - vcfs (str): vcfs(s) passed for SNP analysis
            - panel_pct_coverage (str): total % coverage of panel
            - bootstrap (str): bootstrap to store directly in html
            - version (str): version of Athena, used to add to report

        Returns: None

        Outputs:
            - coverage_report.html (file): HTML coverage report
        """
        print("Generating report")

        # format threshold val & select threshold columns
        threshold_cols = list(cov_stats.filter(regex='[0-9]+x', axis=1))
        vals = ["min", "mean", "max"]
        vals.extend(threshold_cols)

        # generate html formatted list of table headings for tables
        gene_table_headings = [" ", "Gene", "Transcript"] + vals
        exon_table_headings = gene_table_headings.copy()
        exon_table_headings[3:3] = ['Chr', 'Exon', 'Length', 'Start', 'End']

        gene_table_headings = "\n".join(
            [f"<th>{x.capitalize()}</th>" for x in gene_table_headings]
        )
        exon_table_headings = "\n".join(
            [f"<th>{x.capitalize()}</th>" for x in exon_table_headings]
        )

        styling = styleTables(
            cov_stats, cov_summary, self.threshold, threshold_cols, vals
        )
        calculate = calculateValues(self.threshold)

        # apply styling to tables for displaying in report
        sub_threshold_stats, low_exon_columns, gene_issues,\
            exon_issues = styling.style_sub_threshold()

        total_stats = styling.style_total_stats()
        gene_stats, total_genes = styling.style_cov_summary()

        snps_low_cov, snps_not_covered = styling.style_snps_cov(snps_low_cov)
        snps_high_cov, snps_covered = styling.style_snps_cov(snps_high_cov)
        snps_no_cov, snps_out_panel = styling.style_snps_no_cov(snps_no_cov)

        # get values to display in report
        fully_covered_genes = total_genes - gene_issues

        total_snps, snps_pct_covered, snps_pct_not_covered,\
            snps_pct_out_panel = calculate.calculate_snp_vals(
                snps_covered, snps_not_covered, snps_out_panel
            )

        # add values to dict to pass into report
        report_vals = {}
        report_vals["summary_text"] = summary_text
        report_vals["name"] = str(args.sample_name).replace("_", " ")
        report_vals["total_genes"] = str(total_genes)
        report_vals["fully_covered_genes"] = str(fully_covered_genes)
        report_vals["gene_issues"] = str(gene_issues)
        report_vals["threshold"] = self.threshold
        report_vals["gene_table_headings"] = gene_table_headings
        report_vals["exon_table_headings"] = exon_table_headings
        report_vals["exon_issues"] = str(exon_issues)
        report_vals["build"] = build
        report_vals["panel"] = panel
        report_vals["vcfs"] = vcfs
        report_vals["version"] = version
        report_vals["panel_pct_coverage"] = panel_pct_coverage
        report_vals["total_snps"] = total_snps
        report_vals["snps_covered"] = str(snps_covered)
        report_vals["snps_not_covered"] = str(snps_not_covered)
        report_vals["snps_pct_covered"] = str(snps_pct_covered)
        report_vals["snps_pct_not_covered"] = str(snps_pct_not_covered)
        report_vals["snps_out_panel"] = str(snps_out_panel)
        report_vals["snps_pct_out_panel"] = str(snps_pct_out_panel)

        # add tables & plots to template
        html_string = self.build_report(
            html_template, total_stats, gene_stats, sub_threshold_stats,
            low_exon_columns, snps_low_cov, snps_high_cov, snps_no_cov, fig,
            coverage_per_chromosome_fig,
            all_plots, summary_plot, report_vals, bootstrap
        )

        # write output report to file
        self.write_report(html_string, args.output)


    def build_report(self, html_template, total_stats, gene_stats,
                     sub_threshold_stats, low_exon_columns, snps_low_cov, snps_high_cov,
                     snps_no_cov, fig, coverage_per_chromosome_fig, all_plots, summary_plot, report_vals,
                     bootstrap
                     ):
        """
        Build report from template and variables to write to file

        Args:
            - html_template (str): string of HTML template file
            - total_stats (df): total stats table of all genes & exons
            - gene_stats (df): stats table of whole gene
            - sub_threshold_stats (df): table of exons with < threshold
            - snps_low_cov (df): table of snps with cov < threshold
            - snsp_high_cov (df): table of snps with cov > threshold
            - snps_no_cov (df): variants that span exon boundaries (i.e SVs)
            - fig (figure): grid of low coverage exon plots (plotly)
            - coverage_per_chromsome_fig (figure): coverage-per-chromosome plot
            - all-plots (figure): grid of all full gene-exon plots
            - summary_plot (figure): gene summary plot - % at threshold
            - report_vals (dict): values to display in report text
        Returns:
            - single_report (str): HTML string of filled report
        """
        # convert logo image into string to pass in to template
        logo = str(os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "../data/static/images/logo.png"
        ))

        data_uri = base64.b64encode(open(logo, 'rb').read()).decode('utf-8')
        logo = '<img height="25" width="22" src=data:image/png;base64,{0}\
            alt="" style="vertical-align:middle; padding-bottom:3px">'.format(
            data_uri)

        t = Template(html_template)

        date = datetime.today().strftime('%Y-%m-%d')

        single_report = t.safe_substitute(
            bootstrap=bootstrap,
            logo=logo,
            total_genes=report_vals["total_genes"],
            threshold=report_vals["threshold"],
            summary_text=report_vals["summary_text"],
            exon_issues=report_vals["exon_issues"],
            gene_issues=report_vals["gene_issues"],
            fully_covered_genes=report_vals["fully_covered_genes"],
            name=report_vals["name"],
            sub_threshold_stats=sub_threshold_stats,
            low_exon_columns=low_exon_columns,
            low_cov_plots=fig,
            coverage_per_chromosome_fig=coverage_per_chromosome_fig,
            all_plots=all_plots,
            summary_plot=summary_plot,
            gene_stats=gene_stats,
            gene_table_headings=report_vals["gene_table_headings"],
            exon_table_headings=report_vals["exon_table_headings"],
            total_stats=total_stats,
            snps_high_cov_data=snps_high_cov,
            snps_low_cov_data=snps_low_cov,
            snps_no_cov_data=snps_no_cov,
            total_snps=report_vals["total_snps"],
            snps_covered=report_vals["snps_covered"],
            snps_pct_covered=report_vals["snps_pct_covered"],
            snps_not_covered=report_vals["snps_not_covered"],
            snps_pct_not_covered=report_vals["snps_pct_not_covered"],
            snps_out_panel=report_vals["snps_out_panel"],
            snps_pct_out_panel=report_vals["snps_pct_out_panel"],
            date=date,
            build=report_vals["build"],
            vcfs=report_vals["vcfs"],
            panel=report_vals["panel"],
            panel_pct_coverage=report_vals["panel_pct_coverage"],
            version=report_vals["version"]
        )

        return single_report


    def write_report(self, html_string, output_name):
        """
        Write HTML string of populated report to output file
        Args:
            - html_string (str): HTML formatted string of report
            - output_name (str): file name prefix from args.output

        Returns: None

        Outputs:
            - {output_name}_coverage_report.html (file): generated report
        """

        # write report
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        out_dir = os.path.join(bin_dir, "../output/")
        outfile = os.path.join(out_dir, output_name)

        file = open(outfile, 'w')
        file.write(html_string)
        file.close()
        print(f"Output report written to {outfile}")


def load_files(load, threshold, exon_stats, gene_stats, raw_coverage, snp_vcfs, panel, per_base_coverage=None):
    """
    Load in raw coverage data, coverage stats file and template.

    Args:
        - threshold (int): threshold value passed from parse_args
        - exon_stats (file): exon stats file (from args;
                            generated by coverage_stats_single.py)
        - gene_stats (file): gene stats file (from args;
                            generated by coverage_stats_single.py)
        - raw_coverage (file): from args; bp coverage file used as
                            input for coverage_stats_single.py
        - snp_vcfs (list): VCFs of SNPs passed from args
        - panel (file): panel bed file used for annotation, used to
                        display panel name in report if passed
        - per_base_coverage: per base coverage from mosdepth

    Returns:
        - cov_stats (df): df of coverage stats for each exon
        - cov_summary (df): df of gene level coverage
        - raw_coverage (df): raw bp coverage for each exon
        - html_template (str): string of HTML report template
        - flagstat (dict): flagstat metrics, from gene_stats header
        - build (str): ref build used, from gene_stats header
        - panel (str): panes(s) / gene(s) included in report
        - vcfs (str): list of vcf names used for SNP analysis
        - version (str): version of Athena, used to add to report
        - per_base_coverage (df): df of per-base coverage
    """
    print("Reading in files")

    # read in required files
    cov_stats = load.read_exon_stats(exon_stats)
    cov_summary = load.read_gene_stats(gene_stats)
    raw_coverage = load.read_raw_coverage(raw_coverage)
    bootstrap = load.read_bootstrap()
    html_template = load.read_template()
    if per_base_coverage:
        per_base_coverage = load.read_coverage_data(per_base_coverage)

    # get other required attributes
    low_raw_cov = load.get_low_coverage_regions(
        cov_stats, raw_coverage, threshold
    )
    flagstat, build = load.get_build_and_stats(gene_stats)
    version = load.get_athena_ver()
    panel = load.get_panel_name(panel)
    vcfs = load.get_snp_vcfs(snp_vcfs)

    # check if given low coverage threshold is valid
    threshold = load.check_threshold(threshold, cov_stats, cov_summary)

    return cov_stats, cov_summary, raw_coverage, low_raw_cov,\
        html_template, flagstat, build, panel, vcfs, bootstrap, version, \
        per_base_coverage


def parse_args():
    """
    Parse cmd line arguments

    Args: None

    Returns:
        - args (arguments): args passed from cmd line
    """

    parser = argparse.ArgumentParser(
        description='Generate coverage report for a single sample.'
    )
    parser.add_argument(
        '-e', '--exon_stats',
        help='exon stats file (from coverage_stats_single.py)',
        type=argparse.FileType('r'), required=True
    )
    parser.add_argument(
        '-g', '--gene_stats',
        help='gene stats file (from coverage_stats_single.py)',
        required=True
    )
    parser.add_argument(
        '-r', '--raw_coverage',
        help='raw coverage bed file used to generate stats',
        required=True
    )
    parser.add_argument(
        '-b', '--per_base_coverage',
        help='Per-base coverage bed file from mosdepth. If not\
            submitted, plots displaying global coverage per\
            chromosome will not be displayed.',
        default=None
    )
    parser.add_argument(
        '-s', '--snps', nargs='*',
        help='Optional; check coverage of VCF(s) of SNPs.',
        required=False
    )
    parser.add_argument(
        '-t', '--threshold', nargs='?',
        default=20, type=int,
        help="threshold to define low coverage (int), if not\
            given 20 will be used as default. Must be one of\
            the thresholds in the input file.",
        required=False
    )
    parser.add_argument(
        '-n', '--sample_name', nargs='?',
        help="Name of sample to display in report, if not\
            specified this will be the prefix of the\
            gene_stats input file.",
        required=False
    )
    parser.add_argument(
        '-o', '--output', nargs='?',
        help='Output report name, if not specified the sample\
        name from the report will be used.',
        required=False
    )
    parser.add_argument(
        '-p', '--panel', nargs='?',
        help='(Optional) Panel bed file used from annotation, if passed\
        name of file will be displayed in report to show what\
        panel(s) / gene(s) were included.',
        required=False
    )
    parser.add_argument(
        '-l', '--limit', nargs='?',
        help="Number of genes at which to limit including full gene plots,\
        large numbers of genes takes a long time to generate the plots.",
        default=-1,
        required=False
    )
    parser.add_argument(
        '-m', '--summary',
        help="If passed, a short paragraph will be included in the\
        summary section. This includes details on the sequencing and the\
        genes/transcripts used in the panel.",
        default=False, action='store_true'
    )
    parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )

    args = parser.parse_args()

    if not args.sample_name:
        # sample name not given, use input file name
        args.sample_name = Path(args.gene_stats).stem
        if "_" in args.sample_name:
            # if named X1000_ take prefix
            args.sample_name = args.sample_name.split("_", 1)[0]

    if not args.output:
        # output file name not given, using sample name
        args.output = args.sample_name + "_coverage_report.html"
    else:
        args.output = args.output + "_coverage_report.html"

    args.threshold = str(args.threshold) + "x"

    return args


def main():
    """
    Main function to generate coverage report
    """
    args = parse_args()
    load = load_data.loadData()
    calculate = calculateValues(args.threshold)
    plots = generatePlots(args.threshold)
    report = generateReport(args.threshold)

    # read in files
    cov_stats, cov_summary, raw_coverage, low_raw_cov, html_template,\
        flagstat, build, panel, vcfs, bootstrap, version,\
        per_base_coverage = load_files(
            load,
            args.threshold,
            args.exon_stats,
            args.gene_stats,
            args.raw_coverage,
            args.snps,
            args.panel,
            args.per_base_coverage
        )

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

    if args.snps:
        # if SNP VCF(s) have been passed
        snps_low_cov, snps_high_cov, snps_no_cov = calculate.snp_coverage(
            args.snps, raw_coverage
        )
    else:
        # set to empty dfs
        snps_low_cov, snps_high_cov, snps_no_cov = pd.DataFrame(),\
            pd.DataFrame(), pd.DataFrame()

    # calculate mean panel coverage
    panel_pct_coverage = calculate.panel_coverage(cov_stats)

    # generate summary plot
    summary_plot = plots.summary_gene_plot(cov_summary)

    if len(cov_summary.index) < int(args.limit) or int(args.limit) == -1:
        # generate plots of each full gene
        print("Generating full gene plots")
        if num_cores == 1:
            # specified one core, generate plots slowly
            all_plots = plots.all_gene_plots(raw_coverage)
        else:
            raw_coverage = raw_coverage.sort_values(
                ["gene", "tx", "exon"], ascending=[True, True]
            )

            # get unique list of transcripts
            transcripts = raw_coverage.drop_duplicates(
                ["gene"])["gene"].values.tolist()

            # split gene list equally for seperate processes
            tx_array = np.array_split(np.array(transcripts), num_cores)

            # split df into seperate dfs by genes in each list
            split_dfs = np.asanyarray(
                [raw_coverage[raw_coverage["tx"].isin(x)] for x in tx_array],
                dtype=object
            )

            with multiprocessing.Pool(num_cores) as pool:
                # use a pool to spawn multiple processes
                # uses number of cores defined and splits processing of df
                # slices, add each to pool with threshold values
                all_plots = pool.map(plots.all_gene_plots, split_dfs)
                all_plots = "".join(all_plots)
    else:
        all_plots = "<br><b>Full gene plots have been omitted from this report\
            due to the high number of genes in the panel.</b></br>"

    if len(low_raw_cov.index) > 0:
        # some low covered regions, generate plots
        print("Generating plots of low covered regions")

        # get unique list of genes
        transcripts = low_raw_cov.drop_duplicates(["tx"])["tx"].values.tolist()
        print(f"Plots for {len(transcripts)} to generate")

        # split gene list equally for seperate processes
        tx_array = np.array_split(np.array(transcripts), num_cores)

        # split df into seperate dfs by genes in each list
        split_dfs = np.asanyarray(
            [low_raw_cov[low_raw_cov["tx"].isin(x)] for x in tx_array],
            dtype=object
        )

        with multiprocessing.Pool(num_cores) as pool:
            # use a pool to spawn multiple processes
            # uses number of cores defined and splits processing of df
            # slices, add each to pool with threshold values
            fig = pool.map(plots.low_exon_plot, split_dfs)

            # can return None => remove before joining
            fig = [fig_str for fig_str in fig if fig_str]
            fig = ",".join(fig)
    else:
        fig = "<br></br><b>All regions in panel above threshold, no plots\
                to show.</b><br></br>"

    if args.summary:
        # summary text to be included
        summary_text = report.write_summary(
            cov_summary, args.threshold, panel_pct_coverage
        )
    else:
        summary_text = ""

    # generate coverage-per-chromosome figure
    if args.per_base_coverage:
        coverage_per_chromosome_fig = plots.coverage_per_chromosome_plot(per_base_coverage)
    else:
        coverage_per_chromosome_fig = 0

    # generate report
    report.generate_report(
        cov_stats, cov_summary, snps_low_cov, snps_high_cov, snps_no_cov, fig,
        all_plots, coverage_per_chromosome_fig, summary_plot, html_template, args, build, panel, vcfs,
        panel_pct_coverage, bootstrap, version, summary_text
    )


if __name__ == "__main__":

    main()
