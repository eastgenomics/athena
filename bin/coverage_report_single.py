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

import load_data as load


class generatePlots():
    """Functions to generate required plots"""

    def __init__(self, threshold):
        self.threshold = threshold


    def low_exon_plot(self, low_raw_cov):
        """
        Plot bp coverage of exon, used for those where coverage is given
        threshold

        Args:
            - low_raw_cov (df): df of raw coverage for exons with low
                                coverage
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - fig (figure): plots of low coverage regions
        """
        print("Generating plots of low covered regions")

        # get list of tuples of genes and exons to define plots
        genes = low_raw_cov.drop_duplicates(
            ["gene", "exon"])[["gene", "exon"]].values.tolist()
        genes = [tuple(exon) for exon in genes]

        if len(genes) == 0:
            # everything above threshold, don't generate plots
            fig = "<br></br><b>All regions in panel above threshold, no plots\
                to show.</b><br></br>"

            return fig

        # sort list of genes/exons by gene and exon
        genes = sorted(genes, key=lambda element: (element[0], element[1]))

        plot_titles = [str(x[0]) + " exon: " + str(int(x[1])) for x in genes]

        low_raw_cov["exon_len"] =\
            low_raw_cov["exon_end"] - low_raw_cov["exon_start"]

        low_raw_cov["relative_position"] = low_raw_cov["exon_end"] - round(((
            low_raw_cov["cov_end"] + low_raw_cov["cov_start"]) / 2
        ))

        # set no. rows to no. of plots / no of columns to define grid
        columns = 4
        rows = math.ceil(len(genes) / 4)

        # variable height depeendent on no. of plots
        v_space = (1 / rows) * 0.25

        # define grid to add plots to
        fig = make_subplots(
            rows=rows, cols=columns, print_grid=False,
            horizontal_spacing=0.04, vertical_spacing=v_space,
            subplot_titles=plot_titles
        )

        # counter for grid
        row_no = 1
        col_no = 1

        for gene in genes:
            # make plot for each gene / exon

            # counter for grid, by gets to 5th entry starts new row
            if row_no // 5 == 1:
                col_no += 1
                row_no = 1

            # get rows for current gene and exon
            exon_cov = low_raw_cov.loc[(
                low_raw_cov["gene"] == gene[0]
            ) & (
                low_raw_cov["exon"] == gene[1]
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

            # build list of first and last point for threshold line
            xval = [x for x in range(
                exon_cov_unbinned["cov_start"].iloc[0],
                exon_cov_unbinned["cov_end"].iloc[-1]
            )]
            xval = xval[::len(xval) - 1]
            yval = [self.threshold] * 2

            # info field for hovering on plot line
            label = '<i>position: </i>%{x}<br>coverage: %{y}<extra></extra>'

            # generate plot and threshold line to display
            if sum(exon_cov_unbinned["cov"]) != 0:
                plot = go.Scatter(
                    x=exon_cov_unbinned["cov_start"], y=exon_cov_unbinned["cov"],
                    mode="lines",
                    hovertemplate=label
                )
            else:
                # if any plots have no coverage, just display empty plot
                # very hacky way by making data point transparent but
                # ¯\_(ツ)_/¯
                plot = go.Scatter(
                    x=exon_cov_unbinned["cov_start"], y=exon_cov_unbinned["cov"],
                    mode="markers", marker={"opacity": 0}
                )

            threshold_line = go.Scatter(
                x=xval, y=yval, hoverinfo='skip', mode="lines",
                line=dict(color='rgb(205, 12, 24)', width=1)
            )

            # add to subplot grid
            fig.add_trace(plot, col_no, row_no)
            fig.add_trace(threshold_line, col_no, row_no)

            row_no = row_no + 1

        # set height of grid by no. rows and scale value of 325
        height = (rows * 300) + 150

        # update plot formatting
        fig.update_xaxes(nticks=3, ticks="", showgrid=True, tickformat=',d')
        fig.update_yaxes(title='coverage', title_standoff=0)
        fig.update_xaxes(title='exon position', color='#FFFFFF')
        fig["layout"].update(
            height=height, showlegend=False, margin=dict(l=50, r=0)
        )

        # write plots to html string
        fig = fig.to_html(full_html=False)

        return fig


    def all_gene_plots(self, raw_coverage):
        """
        Generate full plots for each gene

        Args:
            - raw_coverage (file): from args; bp coverage file used as
                                    input for coverage_stats_single.py
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - all-plots (figure): grid of all full gene- exon plots
        """

        all_plots = ""

        if len(raw_coverage.index) == 0:
            # passed empty df, most likely because there were less genes
            # than processes => empty df passed with multiprocess
            return ""

        # get unique list of genes
        genes = raw_coverage.drop_duplicates(["gene"])["gene"].values.tolist()

        for gene in genes:

            # get coverage data for current gene
            gene_cov = raw_coverage.loc[(raw_coverage["gene"] == gene)]
            # get list of exons
            exons = gene_cov.drop_duplicates(["exon"])["exon"].values.tolist()

            # no. plot columns = no. of exons
            column_no = len(exons)

            # make subplot grid size of no. of exons, height variable
            # splits large genes to several rows and maintains height
            height = math.ceil(len(exons) / 30) * 4
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

            fig.suptitle(gene, fontweight="bold")
            count = 0

            for exon in exons:
                # get coverage data for current exon
                exon_cov = raw_coverage.loc[(
                    raw_coverage["gene"] == gene
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
                        [0, 100], [self.threshold, self.threshold],
                        color='red', linestyle='-', linewidth=2
                    )
                else:
                    axs[count].plot(exon_cov["cov_start"], exon_cov["cov"])

                    # threshold line
                    axs[count].plot(
                        [exon_cov["exon_start"], exon_cov["exon_end"]],
                        [self.threshold, self.threshold], color='red',
                        linestyle='-', linewidth=1
                    )

                # add labels
                xlab = str(
                    exon_cov["exon_end"].iloc[0] -
                    exon_cov["exon_start"].iloc[0]
                ) + "\nbp"
                axs[count].title.set_text(exon)
                axs[count].set_xlabel(xlab)

                count += 1

            # remove y ticks & label for all but first plot of lines
            for i in range(column_no * rows):
                if i in [x * column_no for x in range(rows)]:
                    # first plot of line, keep ticks and labels
                    continue
                else:
                    axs[i].yaxis.set_ticks_position('none')

            # strip x axis ticks and labels
            plt.setp(plt.gcf().get_axes(), xticks=[])

            # adjust yaxis limits
            ymax = max(gene_cov["cov"].tolist()) + 10
            plt.ylim(bottom=0, top=ymax)

            # remove outer white margins
            fig.tight_layout(h_pad=1.2)

            # convert image to html string and append to one really long
            # string to insert in report
            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)
            image_png = buffer.getvalue()
            buffer.close()
            graphic = base64.b64encode(image_png)
            data_uri = graphic.decode('utf-8')
            img_tag = "<img src=data:image/png;base64,{0} style='max-width:\
                100%; max-height: auto; object-fit: contain; ' />".format(
                data_uri
            )

            all_plots = all_plots + img_tag + "<br></br>"

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

        threshold = str(self.threshold) + "x"

        summary_data = cov_summary.copy()

        # define colours based on values
        summary_data["colours"] = 'green'
        summary_data.loc[
            summary_data[threshold] < 100, 'colours'] = 'orange'
        summary_data.loc[summary_data[threshold] < 90, 'colours'] = 'red'

        summary_data = summary_data.sort_values(
            by=[threshold], ascending=False
        )
        summary_plot, axs = plt.subplots(figsize=(25, 7.5))

        if len(summary_data.index) > 100:
            # split off some of 100% covered genes to limit size of plot
            if len(summary_data[summary_data[threshold] < 100]) > 100:
                # over 100 sub threshold genes, remove all 100% genes
                genes100pct = len(
                    summary_data[summary_data[threshold] == 100]
                )
                summary_data = summary_data[summary_data[threshold] < 100]
            else:
                # split off bottom 100 genes, plot includes some 100% covered
                genes100pct = len(summary_data.iloc[:-100])
                summary_data = summary_data.iloc[-100:]

        plt.bar(
            summary_data["gene"],
            [int(x) for x in summary_data[threshold]],
            color=summary_data.colours
        )

        if "genes100pct" in locals():
            genes100pct = str(genes100pct)
            # more than 100 genes, add title inc. 100% covered not shown
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
        axs.tick_params(labelsize=6, length=0)
        plt.xticks(rotation=55, color="#565656")

        # adjust whole plot marins
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
            bbox_to_anchor=(0.5, -0.1),
            fancybox=True, shadow=True, ncol=12, fontsize=12
        )

        vals = np.arange(0, 110, 10).tolist()
        plt.yticks(vals, vals)
        axs.tick_params(axis='both', which='major', labelsize=8)

        plt.xlabel("")
        plt.ylabel("% coverage ({})".format(self.threshold), fontsize=11)

        axs.yaxis.grid(linewidth=0.5, color="grey", linestyle="-.")
        plt.box(False)
        axs.set_axisbelow(True)
        plt.tight_layout()

        # convert image to html string to insert in report
        buffer = BytesIO()
        plt.savefig(buffer, format='png')
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        graphic = base64.b64encode(image_png)
        data_uri = graphic.decode('utf-8')
        summary_plot = "<img src=data:image/png;base64,{0} style='max-width:\
            100%; max-height: auto; object-fit: contain; ' />".format(
            data_uri
        )

        return summary_plot


class styleTables():
    """Functions for styling tables for displaying in report"""

    def __init__(self, cov_stats, cov_summary, threshold, threshold_cols, vals):
        self.cov_stats = cov_stats
        self.cov_summary = cov_summary
        self.threshold = threshold
        self.threshold_cols = threshold_cols
        self.vals = vals


    def style_sub_threshold(self):
        """
        Styling of sub threshold stats df for displaying in report

        Args:
            - cov_stats (df): df of per exon coverage stats
            - threshold (str): low coverage threshold value
            - threshold_cols (list): threshold values for coverage
        Returns:
            - sub_threshold_stats (str): HTML formatted str of cov stats
                table
            - gene_issues (int): total number of genes under threshold
            - exon_issues (int): total numbner of exons under threshold
        """
        column = [
            "gene", "tx", "chrom", "exon", "exon_len", "exon_start",
            "exon_end", "min", "mean", "max"
        ]

        column.extend(self.threshold_cols)

        sub_threshold = pd.DataFrame(columns=column)

        # get all exons with <100% coverage at threshold
        for i, row in self.cov_stats.iterrows():
            if int(row[self.threshold]) < 100:
                sub_threshold = sub_threshold.append(row, ignore_index=True)

        # pandas is terrible and forces floats, change back to int
        dtypes = {
            'chrom': str,
            'exon': int,
            'exon_len': int,
            'exon_start': int,
            'exon_end': int,
            'min': int,
            'max': int
        }

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

            gene_issues = len(list(set(sub_threshold_stats["gene"].tolist())))
            exon_issues = len(sub_threshold_stats["exon"])
        else:
            # if no low regions set to empty df with appropriate columns
            print("No low coverage regions, generating empty table")
            sub_threshold_stats = pd.DataFrame(columns=column)
            gene_issues = 0
            exon_issues = 0

        # rename columns to display properly
        sub_threshold_stats = sub_threshold_stats.rename(columns={
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
        })

        # reindex & set to begin at 1
        sub_threshold_stats.index = np.arange(
            1, len(sub_threshold_stats.index) + 1
        )

        # creat slices of sub_threshold stats df to add styling to
        slice_ranges = {
            "x0": (10, 0), "x10": (30, 10), "x30": (50, 30), "x50": (70, 50),
            "x70": (90, 70), "x90": (95, 90), "x95": (99, 95), "x99": (101, 99)
        }

        sub_slice = {}

        for key, val in slice_ranges.items():
            sub_slice[key] = pd.IndexSlice[sub_threshold_stats.loc[(
                sub_threshold_stats[self.threshold] < val[0]
            ) & (
                sub_threshold_stats[
                    self.threshold] >= val[1])].index, self.threshold]

        # df column index of threshold
        col_idx = sub_threshold_stats.columns.get_loc(self.threshold)

        # make dict for rounding coverage columns to 2dp
        rnd = {}
        for col in list(sub_threshold_stats.columns[10:]):
            rnd[col] = '{0:.2f}%'

        # set threshold column widths as a fraction of 40% table width
        t_width = str(40 / len(self.threshold_cols)) + "%"

        # apply colours to coverage cell based on value, 0 is given solid red
        s = sub_threshold_stats.style.apply(lambda x: [
            "background-color: #b30000" if x[self.threshold] == 0 and
            idx == col_idx else "" for idx, v in enumerate(x)
        ], axis=1)\
            .bar(subset=sub_slice["x0"], color='#b30000', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x10"], color='#990000', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x30"], color='#C82538', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x50"], color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x70"], color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x90"], color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x95"], color='#FFBF00', vmin=0, vmax=100)\
            .bar(subset=sub_slice["x99"], color='#007600', vmin=0, vmax=100)\
            .format(rnd)\
            .set_table_attributes('table border="1"\
                class="dataframe table table-hover table-bordered"')\
            .set_uuid("low_exon_table")\
            .set_properties(**{'font-size': '0.85vw', 'table-layout': 'auto'})\
            .set_properties(subset=self.threshold_cols, **{'width': t_width})\

        sub_threshold_stats["Mean"] = sub_threshold_stats["Mean"].apply(
            lambda x: int(x)
        )

        sub_threshold_stats = s.render()

        return sub_threshold_stats, gene_issues, exon_issues


    def style_total_stats(self):
        """
        Styling of full gene-exon stats table for displaying in report
        Args:
            - cov_stats (df): df of exon stats
            - threshold_cols (list): list of threshold columns
            - vals (list): list of min, mean and max strs
        Returns:
            -
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

        total_stats = total_stats.rename(columns={
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
        })

        # limit to 2dp using math.floor, use of round() with
        # 2dp may lead to inaccuracy such as 99.99 => 100.00
        round_cols = ['Mean'] + self.threshold_cols

        for col in round_cols:
            total_stats[col] = total_stats[col].map(
                lambda col: math.floor(col * 100) / 100
            )
        # CSS table class for styling tables
        style = (
            '<table border="1" class="dataframe">',
            '<table class="table table-striped" style="font-size: 0.85vw;" >'
        )

        total_stats = total_stats.to_html(justify='left').replace(
            style[0], style[1]
        )

        return total_stats


    def style_cov_summary(self):
        """
        Add styling to per gene coverage summary table
        Args:
            - cov_summary (df): df of gene coverage stats
            - threshold_cols (list): list of threshold values
        Returns:
            - gene_stats (str): HTML formatted str of gene summary df
            - total_genes (int): total number of genes
        """
        # rename columns for displaying in report
        cov_summary = self.cov_summary.drop(columns=["exon"])
        cov_summary = cov_summary.rename(columns={
            "gene": "Gene",
            "tx": "Transcript",
            "min": "Min",
            "mean": "Mean",
            "max": "Max"
        })

        # get values to display in report
        total_genes = len(cov_summary["Gene"].tolist())

        # limit to 2dp using math.floor, use of round() with
        # 2dp may lead to inaccuracy such as 99.99 => 100.00
        round_cols = ['Mean'] + self.threshold_cols

        for col in round_cols:
            cov_summary[col] = cov_summary[col].map(
                lambda col: math.floor(col * 100) / 100
            )

        # reset index to start at 1
        cov_summary.index = np.arange(1, len(cov_summary.index) + 1)

        # CSS table class for styling tables
        style = (
            '<table border="1" class="dataframe">',
            '<table class="table table-striped" style="font-size: 0.85vw;" >'
        )

        # generate HTML strings from table objects to write to report
        gene_stats = cov_summary.to_html(justify='left').replace(
            style[0], style[1]
        )

        return gene_stats, total_genes


    @staticmethod
    def style_snps_low_cov(snps_low_cov):
        """
        Add styling to table of snps under coverage threshold
        Args:
            - snps_low_cov (df): df of snps under covegrage threshold
        Returns:
            - snps_low_cov (str): HTML formatted str of low covered snps
            - snps_not_covered (int): total number snps not covered at
                threshold
        """
        # get snps values and format dfs to display
        if not snps_low_cov.empty:
            # format low coverage SNPs table
            snps_low_cov.index = np.arange(1, len(snps_low_cov.index) + 1)
            snps_not_covered = len(snps_low_cov.index)
            snps_low_cov = snps_low_cov.style\
                .set_table_attributes(
                    'class="dataframe table table-striped"')\
                .set_uuid("var_low_cov")\
                .set_properties(**{
                    'font-size': '0.80vw', 'table-layout': 'auto'
                })\
                .set_properties(subset=["VCF", "Gene"], **{'width': '10%'})\
                .set_properties(subset=["Exon"], **{'width': '7.5%'})\
                .set_properties(subset=["Chromosome"], **{'width': '10%'})\
                .set_properties(subset=["Position"], **{'width': '12.5%'})\
                .set_properties(subset=["Ref"], **{'width': '20%'})\
                .set_properties(subset=["Alt"], **{'width': '20%'})\
                .set_properties(subset=["Coverage"], **{'width': '10%'})

            snps_low_cov = snps_low_cov.render()
        else:
            snps_low_cov = "<b>No low covered SNPs</b>"
            snps_not_covered = 0

        return snps_low_cov, snps_not_covered


    @staticmethod
    def style_snps_high_cov(snps_high_cov):
        """
        Add styling to table of SNPs covered above threshold
        Args:
            - snps_high_cov (df): df of snps covered above threshold
        Returns:
            - snps_high_cov (str): HTML formatted str of covered snps
            - snps_covered (int): total number of snps covered
        """

        if not snps_high_cov.empty:
            # format high coverage SNPs table
            snps_high_cov.index = np.arange(1, len(snps_high_cov.index) + 1)

            snps_covered = len(snps_high_cov.index)

            snps_high_cov = snps_high_cov.style\
                .set_table_attributes(
                    'class="dataframe table table-striped"')\
                .set_uuid("var_high_cov")\
                .set_properties(**{
                    'font-size': '0.80vw', 'table-layout': 'auto'
                })\
                .set_properties(subset=["VCF", "Gene"], **{'width': '10%'})\
                .set_properties(subset=["Exon"], **{'width': '7.5%'})\
                .set_properties(subset=["Chromosome"], **{'width': '10%'})\
                .set_properties(subset=["Position"], **{'width': '12.5%'})\
                .set_properties(subset=["Ref"], **{'width': '20%'})\
                .set_properties(subset=["Alt"], **{'width': '20%'})\
                .set_properties(subset=["Coverage"], **{'width': '10%'})

            snps_high_cov = snps_high_cov.render()
        else:
            snps_high_cov = "<b>No covered SNPs</b>"
            snps_covered = 0

        return snps_high_cov, snps_covered


    @staticmethod
    def style_snps_no_cov(snps_no_cov):
        """
        Add styling to table of snps that span exon boundaries => have
        coverage values
        Args:
            - snps_no_cov (df): df of snps with no coverage values
        Returns:
            - snps_no_cov (str): HTML formatted str of snps with no cov
            - snps_out_panel (int): total number snps with no cov
        """
        # if variants from vcf found that span exon boundaries
        if not snps_no_cov.empty:
            # manually add div and styling around rendered table, allows
            # to be fully absent from the report if the table is empty
            snps_no_cov.index = np.arange(1, len(snps_no_cov) + 1)

            # get number of variants to display in report
            snps_out_panel = len(snps_no_cov.index)

            html_string = snps_no_cov.style\
                .set_table_attributes(
                    'class="dataframe table table-striped"')\
                .set_uuid("var_no_cov")\
                .set_properties(**{
                    'font-size': '0.80vw', 'table-layout': 'auto'
                })\
                .set_properties(subset=["VCF"], **{
                    'width': '7.5%'
                })\
                .set_properties(subset=[
                    "Chromosome", "Position", "Ref", "Alt"
                ], **{'width': '10%'})

            html_string = html_string.render()

            snps_no_cov = """
                <br> Variants included in the first table below either fully\
                    or partially span panel region(s). These are most likely\
                    large structural variants and as such do not have\
                    coverage data available. See the "info" column for details\
                    on the variant.
                </br>
                <br> Table of variants spanning panel regions(s) &nbsp
                <button class="btn btn-info collapsible btn-sm">Show /\
                     hide table</button>
                <div class="content">
                    <table>
                        {}
                    </table>
                </div></br>
                """.format(html_string)
        else:
            snps_no_cov = ""
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

        # threshold column to check at
        threshold = str(self.threshold) + "x"

        gene_stats = pd.DataFrame(
            columns=["gene", "gene_len", "coverage"])

        # make list of genes
        genes = sorted(list(set(cov_stats["gene"].tolist())))

        for gene in genes:
            # for each gene, calculate length and average % at threshold
            gene_cov = cov_stats.loc[cov_stats["gene"] == gene]

            length = sum(gene_cov["exon_len"])
            coverage = sum(
                gene_cov[threshold] * gene_cov["exon_len"] / length)

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
        snps_low_cov = snps_cov.loc[snps_cov["Coverage"] < self.threshold]
        snps_high_cov = snps_cov.loc[snps_cov["Coverage"] >= self.threshold]

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
        threshold = str(threshold) + "x"

        pct_cov = str(math.floor(float(panel_pct_coverage)))

        # summary text paragraph with div styling
        summary_text = """
        <li>Clinical report summary:</li>
        <div style="background-color:aliceblue; margin-top: 15px;
        border-radius: 15px; padding-left:25px;">
        <div id="summary_text" style="font-size: 14px;
        padding-bottom: 15px; padding-top:10px">"""

        sub90 = ""

        for idx, gene in cov_summary.iterrows():
            # build string of each gene, trascript and coverage at
            # threshold to display in summary
            summary = "{} ({}); ".format(gene["gene"], gene["tx"])
            summary_text += summary

            if gene[threshold] < 90:
                # build string of genes with <90% coverage at threshold
                sub90 += "{} ({}); ".format(gene["gene"], gene["tx"])

        summary_text = summary_text.strip(" ;") + "."
        sub90 = sub90.strip(" ;") + "."

        summary_text += """
            <br></br>Genes with coverage at {} less than 90%:
            {}""".format(threshold, sub90)

        summary_text += """
            <br></br>{} % of this panel was sequenced to a depth of {} or
            greater.<br>""".format(pct_cov, threshold)

        # add closing div and copy button for summary text
        summary_text += """</div><div style="padding-bottom:15px;">
        <button class="btn-info btn-sm summarybtn" onclick=
        "CopyToClipboard('summary_text')";return false; style="font-size: 14px;
        padding:5px 10px; border-radius: 10px;">Copy summary text
        </button></div></div>"""

        return summary_text


    def generate_report(self, cov_stats, cov_summary, snps_low_cov,
                        snps_high_cov, snps_no_cov, fig, all_plots,
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
        threshold = str(args.threshold) + "x"
        threshold_cols = list(cov_stats.filter(regex='[0-9]+x', axis=1))
        vals = ["min", "mean", "max"]
        vals.extend(threshold_cols)

        styling = styleTables(
            cov_stats, cov_summary, threshold, threshold_cols, vals
        )
        calculate = calculateValues(threshold)

        # apply styling to tables for displaying in report
        sub_threshold_stats, gene_issues,\
            exon_issues = styling.style_sub_threshold()

        total_stats = styling.style_total_stats()
        gene_stats, total_genes = styling.style_cov_summary()

        snps_low_cov, snps_not_covered = styling.style_snps_low_cov(
            snps_low_cov
        )

        snps_high_cov, snps_covered = styling.style_snps_high_cov(
            snps_high_cov
        )

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
        report_vals["threshold"] = threshold
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
            snps_low_cov, snps_high_cov, snps_no_cov, fig, all_plots,
            summary_plot, report_vals, bootstrap
        )

        # write output report to file
        self.write_report(html_string, args.output)


    def build_report(self, html_template, total_stats, gene_stats,
                     sub_threshold_stats, snps_low_cov, snps_high_cov,
                     snps_no_cov, fig, all_plots, summary_plot, report_vals,
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
            - all-plots (figure): grid of all full gene- exon plots
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
            low_cov_plots=fig,
            all_plots=all_plots,
            summary_plot=summary_plot,
            gene_stats=gene_stats,
            total_stats=total_stats,
            snps_high_cov=snps_high_cov,
            snps_low_cov=snps_low_cov,
            snps_no_cov=snps_no_cov,
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


def load_files(
        threshold, exon_stats, gene_stats, raw_coverage, snp_vcfs, panel):
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
    """
    print("Reading in files")

    # read in required files
    cov_stats = load.read_exon_stats(exon_stats)
    cov_summary = load.read_gene_stats(gene_stats)
    raw_coverage = load.read_raw_coverage(raw_coverage)
    bootstrap = load.read_bootstrap()
    html_template = load.read_template()

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
        html_template, flagstat, build, panel, vcfs, bootstrap, version


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

    return args


def main():
    """
    Main function to generate coverage report
    """
    args = parse_args()
    calculate = calculateValues(args.threshold)
    plots = generatePlots(args.threshold)
    report = generateReport()

    # read in files
    cov_stats, cov_summary, raw_coverage, low_raw_cov, html_template,\
        flagstat, build, panel, vcfs, bootstrap, version = load_files(
            args.threshold,
            args.exon_stats,
            args.gene_stats,
            args.raw_coverage,
            args.snps,
            args.panel
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

    # generate plot of sub optimal regions
    fig = plots.low_exon_plot(low_raw_cov)

    if len(cov_summary.index) < int(args.limit) or int(args.limit) == -1:
        # generate plots of each full gene
        print("Generating full gene plots")

        raw_coverage = raw_coverage.sort_values(
            ["gene", "exon"], ascending=[True, True]
        )

        # get unique list of genes
        genes = raw_coverage.drop_duplicates(["gene"])["gene"].values.tolist()

        # split gene list equally for seperate processes
        gene_array = np.array_split(np.array(genes), num_cores)

        # split df into seperate dfs by genes in each list
        split_dfs = np.asanyarray(
            [raw_coverage[raw_coverage["gene"].isin(x)] for x in gene_array],
            dtype=object
        )

        with multiprocessing.Pool(num_cores) as pool:
            # use a pool to spawn multiple processes
            # uses number of cores defined and splits processing of df
            # slices, add each to pool with threshold values and
            # concatenates together when finished
            all_plots = ''.join(pool.map(plots.all_gene_plots, split_dfs))
    else:
        all_plots = "<br><b>Full gene plots have been omitted from this report\
            due to the high number of genes in the panel.</b></br>"

    if args.summary:
        # summary text to be included
        summary_text = report.write_summary(
            cov_summary, args.threshold, panel_pct_coverage
        )
    else:
        summary_text = ""

    # generate report
    report.generate_report(
        cov_stats, cov_summary, snps_low_cov, snps_high_cov, snps_no_cov, fig,
        all_plots, summary_plot, html_template, args, build, panel, vcfs,
        panel_pct_coverage, bootstrap, version, summary_text
    )


if __name__ == "__main__":

    main()
