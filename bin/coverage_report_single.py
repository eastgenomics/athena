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

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly.graph_objs as go
import sys
import tempfile

from datetime import datetime
from io import BytesIO
from plotly.subplots import make_subplots
from string import Template


class singleReport():

    def load_files(self, threshold, exon_stats,
                   gene_stats, raw_coverage, snp_vcfs):
        """
        Load in raw coverage data, coverage stats file and template.

        Args:
            - exon_stats (file): exon stats file (from args;
                                generated by coverage_stats_single.py)
            - gene_stats (file): gene stats file (from args;
                                generated by coverage_stats_single.py)
            - raw_coverage (file): from args; bp coverage file used as
                                input for coverage_stats_single.py
            - snp_vcfs (list):

        Returns:
            - cov_stats (df): df of coverage stats for each exon
            - cov_summary (df): df of gene level coverage
            - raw_coverage (df): raw bp coverage for each exon
            - html_template (str): string of HTML report template
            - flagstat (dict): flagstat metrics, from gene_stats header
            - build (str): ref build used, from gene_stats header
        """

        print("Reading in files")

        # read in single sample report template
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(bin_dir, "../data/templates/")
        single_template = os.path.join(template_dir, "single_template.html")

        with open(single_template, 'r') as temp:
            html_template = temp.read()

        # read in exon stats file
        with open(exon_stats.name) as exon_file:
            cov_stats = pd.read_csv(exon_file, sep="\t", comment='#')

        # read in gene stats file
        with open(gene_stats) as gene_file:
            cov_summary = pd.read_csv(gene_file, sep="\t", comment='#')

        flagstat = {}
        # read in flagstat and build from header of gene stats file
        with open(gene_stats) as gene_file:
            for ln in gene_file:
                if ln.startswith("#"):
                    if "build" in ln:
                        # get build number
                        build = ln.split(":")[1]
                    else:
                        # read in flagstat from header
                        key = ln.split(":")[0].strip("#")
                        val = ln.split(":")[1]
                        flagstat[key] = val

        column = [
            "chrom", "exon_start", "exon_end",
            "gene", "tx", "exon", "cov_start",
            "cov_end", "cov"
        ]

        # read in raw coverage stats file
        with open(raw_coverage) as raw_file:
            raw_coverage = pd.read_csv(raw_file, sep="\t", names=column)

        if snp_vcfs:
            # SNP vcfs(s) passed
            # read in all VCF(s) and concatenate into one df
            header = ["chrom", "snp_pos", "ref", "alt"]
            snp_df = pd.concat((pd.read_csv(
                f, sep="\t", usecols=[0, 1, 3, 4], comment='#',
                low_memory=False, header=None, names=header) for f in snp_vcfs
            ))
        else:
            snp_df = None

        # check given threshold is in the stats files
        threshold = str(threshold) + "x"

        if threshold not in list(cov_stats) and\
                threshold not in list(cov_summary):
            print("--threshold must be one of the gene and exon\
                    stats coverage thresholds. Exiting now.")
            sys.exit()

        return cov_stats, cov_summary, snp_df, raw_coverage,\
            html_template, build


    def build_report(self, html_template, total_stats, gene_stats,
                     sub_thrshld_stats, snps_low_cov, snps_high_cov, fig,
                     all_plots, summary_plot, report_vals
                     ):
        """
        Build report from template and variables to write to file

        Args:
            - html_template (str): string of HTML template
            - total_stats (df): total stats table of all genes & exons
            - gene_stats (df): stats table of whole gene
            - sub_thrshld_stats (df): table of exons with < threshold
            - snps_low_cov (df): table of snps with cov < threshold
            - snsp_high_cov (df): table of snps with cov > threshold
            - fig (figure): grid of low coverage exon plots (plotly)
            - all-plots (figure): grid of all full gene- exon plots
            - summary_plot (figure): gene summary plot - % at threshold
            - report_vals (dict): values to display in report text
        Returns:
            - single_report (str): HTML string of filled report
        """

        t = Template(html_template)

        date = datetime.today().strftime('%Y-%m-%d')

        single_report = t.safe_substitute(
            total_genes=report_vals["total_genes"],
            threshold=report_vals["threshold"],
            exon_issues=report_vals["exon_issues"],
            gene_issues=report_vals["gene_issues"],
            covered_genes=report_vals["covered_genes"],
            name=report_vals["name"],
            sub_thrshld_stats=sub_thrshld_stats,
            low_cov_plots=fig,
            all_plots=all_plots,
            summary_plot=summary_plot,
            gene_stats=gene_stats,
            total_stats=total_stats,
            snps_high_cov=snps_high_cov,
            snps_low_cov=snps_low_cov,
            total_snps=report_vals["total_snps"],
            snps_covered=report_vals["snps_covered"],
            snps_pct_covered=report_vals["snps_pct_covered"],
            snps_not_covered=report_vals["snps_not_covered"],
            snps_pct_not_covered=report_vals["snps_pct_not_covered"],
            date=date,
            build=report_vals["build"]
        )

        return single_report


    def snp_coverage(self, snp_df, raw_coverage, threshold):
        """
        Produces table of coverage for SNPs inside of capture regions.

        Args:
            - snp_df (df): df of all SNPs from input VCF(s)
            - raw_coverage (df): raw bp coverage for each exon

        Returns:
            - snps_low_cov (df): SNPs with lower coverage than threshold
            - snps_high_cov (df): SNPs with higher coverage than threshold
        """
        print("Calculating coverage of given SNPs")

        # reset indexes
        snp_df = snp_df.reset_index(drop=True)
        raw_coverage = raw_coverage.reset_index(drop=True)

        # select unique exons coordinates, coverage seperated due to size
        exons = raw_coverage[["chrom", "exon_start", "exon_end"]]\
            .drop_duplicates().reset_index(drop=True)

        exons_cov = raw_coverage[[
            "gene", "exon", "chrom", "exon_start", "exon_end", "cov"
        ]].drop_duplicates().reset_index(drop=True)

        exons["chrom"] = exons["chrom"].astype(str)
        exons_cov["chrom"] = exons_cov["chrom"].astype(str)

        # intersect all SNPs against exons to find those SNPs in capture
        snps = exons.merge(snp_df, on='chrom', how='left')
        snps = snps[
            (snps.snp_pos >= snps.exon_start) & (snps.snp_pos <= snps.exon_end)
        ]

        snps = snps[["chrom", "snp_pos", "ref", "alt"]].reset_index(drop=True)

        # add coverage data back to df of snps in capture
        # uses less ram than performing in one go
        snp_cov = snps.merge(exons_cov, on='chrom', how='left')

        snps_cov = snp_cov[
            ["gene", "exon", "chrom", "snp_pos", "ref", "alt", "cov"]
        ].drop_duplicates(
            subset=["chrom", "snp_pos", "ref", "alt"]
        ).reset_index(drop=True)

        # rename columns for displaying in report
        snps_cov.columns = ["Gene", "Exon", "Chromosome", "Position",
                            "Ref", "Alt", "Coverage"]

        snps_cov["Coverage"] = snps_cov["Coverage"].astype(int)

        # split SNPs by coverage against threshold
        snps_low_cov = snps_cov.loc[snps_cov["Coverage"] < threshold]
        snps_high_cov = snps_cov.loc[snps_cov["Coverage"] >= threshold]

        return snps_low_cov, snps_high_cov


    def low_coverage_regions(self, cov_stats, raw_coverage, threshold):
        """
        Get regions where coverage at given threshold is <100%

        Args:
            - cov_stats (df): df of coverage stats for each exon
            - raw_coverage (df): raw bp coverage for each exon
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - low_raw_cov (df): df of raw bp values for each region with
                                coverage less than 100% at threshold
        """
        # threshold column to check at
        threshold = str(threshold) + "x"

        # get threshold columns and add to column names
        threshold_cols = list(cov_stats.filter(regex='[0-9]+x', axis=1))

        column = [
            "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
            "min", "mean", "max"
        ]

        column.extend(threshold_cols)

        # empty df
        low_stats = pd.DataFrame(columns=column)

        # get all exons with <100% coverage at given threshold
        for i, row in cov_stats.iterrows():
            if int(row[threshold]) < 100:
                low_stats = low_stats.append(row, ignore_index=True)

        # pandas is terrible and forces floats, change back to int
        dtypes = {
            'chrom': int,
            'exon': int,
            'exon_start': int,
            'exon_end': int,
            'min': int,
            'max': int
        }

        low_stats = low_stats.astype(dtypes)

        # get list of tuples of genes and exons with low coverage to
        # select out raw coverage
        low_exon_list = low_stats.reset_index()[['gene',
                                                'exon']].values.tolist()
        low_exon_list = [tuple(exon) for exon in low_exon_list]

        # get raw coverage for low coverage regions to plot
        low_raw_cov = raw_coverage[raw_coverage[['gene', 'exon']].apply(
            tuple, axis=1).isin(low_exon_list)].reset_index()

        return low_raw_cov


    def low_exon_plot(self, low_raw_cov, threshold):
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

        # sort list of genes/exons by gene and exon
        genes = sorted(genes, key=lambda element: (element[0], element[1]))

        plot_titles = [str(x[0]) + " exon: " + str(x[1]) for x in genes]

        low_raw_cov["exon_len"] =\
            low_raw_cov["exon_end"] - low_raw_cov["exon_start"]

        low_raw_cov["relative_position"] = low_raw_cov["exon_end"] - round(((
            low_raw_cov["cov_end"] + low_raw_cov["cov_start"]) / 2
        ))

        # set no. rows to no. of plots / no of columns to define grid
        columns = 4
        rows = math.ceil(len(genes) / 4)

        # variable height dependent on no. of plots
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

            # build list of first and last point for line
            xval = [x for x in range(
                exon_cov["cov_start"].iloc[0],
                exon_cov["cov_end"].iloc[-1]
            )]
            xval = xval[::len(xval) - 1]
            yval = [threshold] * 2

            # generate plot and threshold line to display
            if sum(exon_cov["cov"]) != 0:
                plot = go.Scatter(
                    x=exon_cov["cov_start"], y=exon_cov["cov"],
                    mode="lines",
                    hovertemplate='<i>position: </i>%{x}' +
                                  '<br>coverage: %{y}<br>',
                )
            else:
                # if any plots have no coverage, just display empty plot
                # very hacky way by making data point transparent but
                # ¯\_(ツ)_/¯
                plot = go.Scatter(
                    x=exon_cov["cov_start"], y=exon_cov["cov"],
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

        # set height of grid by no. rows and scale value of 275
        height = rows * 275

        # update plot formatting
        fig["layout"].update(height=height, showlegend=False)
        fig.update_xaxes(nticks=3, ticks="", showgrid=True, tickformat=',d')
        fig.update_yaxes(title='coverage')
        fig.update_xaxes(title='exon position', color='#FFFFFF')

        # write plots to html string
        fig = fig.to_html(full_html=False)

        return fig


    def all_gene_plots(self, raw_coverage, threshold):
        """
        Generate full plots for each gene

        Args:
            - raw_coverage (file): from args; bp coverage file used as
                                    input for coverage_stats_single.py
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - all-plots (figure): grid of all full gene- exon plots
        """
        print("Generating full gene plots")

        raw_coverage = raw_coverage.sort_values(
            ["gene", "exon"], ascending=[True, True]
        )
        genes = raw_coverage.drop_duplicates(["gene"])["gene"].values.tolist()

        all_plots = ""

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
            fig = plt.figure(figsize=(20, height))

            if column_no == 1:
                # handle genes with single exon and not using subplots
                plt.plot(exon_cov["cov_start"], exon_cov["cov"])
                plt.plot(
                    [exon_cov["exon_start"], exon_cov["exon_end"]],
                    [threshold, threshold], color='red', linestyle='-',
                    linewidth=1
                )
                plt.xticks([])

                ymax = max(gene_cov["cov"].tolist()) + 10
                plt.ylim(bottom=0, top=ymax)

                xlab = str(
                    exon_cov["exon_end"].iloc[0] -
                    exon_cov["exon_start"].iloc[0]
                ) + "bp"

                plt.xlabel(xlab)

                title = gene + "; exon " + str(exon)
                fig.suptitle(title)

            else:
                # generate grid with space for each exon
                # splits genes with >25 exons to multiple rows
                rows = math.ceil(len(exons) / 30)
                if column_no > 30:
                    column_no = 30

                grid = fig.add_gridspec(rows, column_no, wspace=0)
                axs = grid.subplots(sharey=True)
                axs = axs.flatten()

                fig.suptitle(gene)
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
                    exon_cov = exon_cov.sort_values(
                        by='cov_start', ascending=True
                    )

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
                        # no coverage, generate empty plot
                        axs[count].plot([0, 0], [0, 0])
                    else:
                        axs[count].plot(
                            exon_cov["cov_start"], exon_cov["cov"]
                        )

                    # threshold line
                    axs[count].plot(
                        [exon_cov["exon_start"], exon_cov["exon_end"]],
                        [threshold, threshold], color='red', linestyle='-',
                        linewidth=1
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


    def summary_gene_plot(self, cov_summary, threshold):
        """
        Generate summary plot of all genes against threshold value

        Args:
            - cov_summary (df): df of gene coverage values
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - summary_plot (fig): plot of all genes
        """

        print("Generating summary plot")

        thrshld = str(threshold) + "x"

        # define colours based on values
        cov_summary["colours"] = 'green'
        cov_summary.loc[cov_summary[thrshld] < 100, 'colours'] = 'orange'
        cov_summary.loc[cov_summary[thrshld] < 90, 'colours'] = 'red'

        cov_summary = cov_summary.sort_values(by=[thrshld], ascending=False)

        summary_plot, axs = plt.subplots(figsize=(18, 10))
        plt.bar(
            cov_summary["gene"], [int(x) for x in cov_summary[thrshld]],
            color=cov_summary.colours
        )

        # threshold lines
        plt.axhline(y=99, linestyle='--', color="#565656", alpha=0.6)
        plt.axhline(y=95, linestyle='--', color="#565656", alpha=0.6)

        plt.text(1.02, 0.93, '99%', transform=axs.transAxes)
        plt.text(1.02, 0.90, '95%', transform=axs.transAxes)

        # plot formatting
        axs.tick_params(labelsize=6, length=0)
        plt.xticks(rotation=35, color="#565656")

        vals = np.arange(0, 110, 10).tolist()
        plt.yticks(vals, vals)
        axs.tick_params(axis='both', which='major', labelsize=10)

        plt.xlabel("")
        plt.ylabel("% coverage >= {}".format(thrshld))

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

        cov_summary = cov_summary.drop(columns=['colours'])
        cov_summary = cov_summary.sort_values(["gene"], ascending=[True])
        cov_summary = cov_summary.reset_index()

        return summary_plot, cov_summary


    def generate_report(self, cov_stats, cov_summary, snps_low_cov,
                        snps_high_cov, fig, all_plots, summary_plot,
                        html_template, args, build
                        ):
        """
        Generate single sample report from coverage stats

        Args:
            - cov_stats (df): df of coverage stats for each exon
            - cov_summary (df): df of gene level coverage
            - snps_low_cov (df): SNPs with lower coverage than threshold
            - snps_high_cov (df): SNPs with higher coverage than threshold
            - fig (figure): plots of low coverage regions
            - all-plots (figure): grid of all full gene- exon plots
            - summary_plot (figure): gene summary plot - % at threshold
            - html_template (str): string of HTML template
            - args (args): passed cmd line arguments
            - build (str): build number used for alignment

        Returns: None

        Outputs:
            - coverage_report.html (file): HTML coverage report
        """
        print("Generating report")

        # str of threshold for selecting df columns etc.
        thrshld = str(args.threshold) + "x"

        # get threshold columns and add to column names
        threshold_cols = list(cov_stats.filter(regex='[0-9]+x', axis=1))

        column = [
            "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
            "min", "mean", "max"
        ]

        column.extend(threshold_cols)

        sub_thrshld = pd.DataFrame(columns=column)

        # get all exons with <100% coverage at threshold
        for i, row in cov_stats.iterrows():
            if int(row[thrshld]) < 100:
                sub_thrshld = sub_thrshld.append(row, ignore_index=True)

        # pandas is terrible and forces floats, change back to int
        dtypes = {
            'chrom': int,
            'exon': int,
            'exon_len': int,
            'exon_start': int,
            'exon_end': int,
            'min': int,
            'max': int
        }

        sub_thrshld = sub_thrshld.astype(dtypes)

        vals = ["min", "mean", "max"]
        vals.extend(threshold_cols)

        # do some excel level formatting to make table more readable
        total_stats = pd.pivot_table(
            cov_stats,
            index=["gene", "tx", "chrom", "exon", "exon_len",
                   "exon_start", "exon_end"],
            values=vals
        )

        sub_thrshld_stats = pd.pivot_table(
            sub_thrshld,
            index=["gene", "tx", "chrom", "exon", "exon_len",
                   "exon_start", "exon_end"],
            values=vals
        )

        # reset index to fix formatting
        total_stats = total_stats.reindex(vals, axis=1)
        sub_thrshld_stats = sub_thrshld_stats.reindex(vals, axis=1)
        total_stats.reset_index(inplace=True)
        sub_thrshld_stats.reset_index(inplace=True)

        # rename columns to display properly
        sub_thrshld_stats = sub_thrshld_stats.rename(columns={
            "gene": "Gene",
            "tx": "Transcript",
            "chrom": "Chromosome",
            "exon": "Exon",
            "exon_len": "Length",
            "exon_start": "Start",
            "exon_end": "End",
            "min": "Min",
            "mean": "Mean",
            "max": "Max"
        })

        cov_summary = cov_summary.drop(columns=["index", "exon"])
        cov_summary = cov_summary.rename(columns={
            "gene": "Gene",
            "tx": "Transcript",
            "min": "Min",
            "mean": "Mean",
            "max": "Max"
        })

        total_stats = total_stats.rename(columns={
            "gene": "Gene",
            "tx": "Transcript",
            "chrom": "Chromosome",
            "exon": "Exon",
            "exon_len": "Length",
            "exon_start": "Start",
            "exon_end": "End",
            "min": "Min",
            "mean": "Mean",
            "max": "Max"
        })

        # get values to display in report
        total_genes = len(cov_summary["Gene"])
        gene_issues = len(list(set(sub_thrshld_stats["Gene"].tolist())))
        exon_issues = len(sub_thrshld_stats["Exon"])
        covered_genes = total_genes - gene_issues

        # empty dict to add values for displaying in report text
        report_vals = {}

        report_vals["name"] = str(args.sample_name)
        report_vals["total_genes"] = str(total_genes)
        report_vals["covered_genes"] = str(covered_genes)
        report_vals["gene_issues"] = str(gene_issues)
        report_vals["threshold"] = thrshld
        report_vals["exon_issues"] = str(exon_issues)
        report_vals["build"] = build

        # set ranges for colouring cells
        x0 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 10
        ) & (
            sub_thrshld_stats[thrshld] > 0)].index, thrshld]
        x10 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 30
        ) & (
            sub_thrshld_stats[thrshld] >= 10)].index, thrshld]
        x30 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 50
        ) & (
            sub_thrshld_stats[thrshld] >= 30)].index, thrshld]
        x50 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 70
        ) & (
            sub_thrshld_stats[thrshld] >= 50)].index, thrshld]
        x70 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 90
        ) & (
            sub_thrshld_stats[thrshld] >= 70)].index, thrshld]
        x90 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] < 95
        ) & (
            sub_thrshld_stats[thrshld] >= 90)].index, thrshld]
        x95 = pd.IndexSlice[sub_thrshld_stats.loc[(
            sub_thrshld_stats[thrshld] >= 95)].index, thrshld]

        # df column index of threshold
        col_idx = sub_thrshld_stats.columns.get_loc(thrshld)

        # make dict for rounding coverage columns to 2dp
        rnd = {}
        for col in list(sub_thrshld_stats.columns[10:15]):
            rnd[col] = '{0:.2f}%'

        # apply colours to coverage cell based on value, 0 is given solid red
        s = sub_thrshld_stats.style.apply(lambda x: [
            "background-color: #d70000" if x[thrshld] == 0 and idx == col_idx
            else "" for idx, v in enumerate(x)
        ], axis=1)\
            .bar(subset=x0, color='red', vmin=0, vmax=100)\
            .bar(subset=x10, color='#990000', vmin=0, vmax=100)\
            .bar(subset=x30, color='#C82538', vmin=0, vmax=100)\
            .bar(subset=x50, color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=x70, color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=x90, color='#45731E', vmin=0, vmax=100)\
            .bar(subset=x95, color='#007600', vmin=0, vmax=100)\
            .format(rnd)\
            .set_table_attributes('table border="1"\
                class="dataframe table table-hover table-bordered"')

        sub_thrshld_stats["Mean"] = sub_thrshld_stats["Mean"].apply(
            lambda x: int(x)
        )

        # CSS table class for styling tables
        style = (
            '<table border="1" class="dataframe">',
            '<table class="table table-striped">'
        )

        # generate HTML strings from table objects to write to report
        gene_stats = cov_summary.to_html(justify='left').replace(
            style[0], style[1]
        )
        total_stats = total_stats.to_html(justify='left').replace(
            style[0], style[1]
        )
        sub_thrshld_stats = s.render()

        if snps_low_cov is not None:
            snps_not_covered = len(snps_low_cov.index)
            snps_low_cov = snps_low_cov.to_html().replace(style[0], style[1])
        else:
            snps_low_cov = "<b>No SNPs present</b>"
            snps_not_covered = 0

        if snps_high_cov is not None:
            snps_covered = len(snps_high_cov.index)
            snps_high_cov = snps_high_cov.to_html().replace(style[0], style[1])
        else:
            snps_high_cov = "<b>No SNPs present</b>"
            snps_covered = 0

        total_snps = str(snps_covered + snps_not_covered)

        snps_pct_covered = int(snps_covered) - int(total_snps) * 100
        snps_pct_not_covered = int(snps_not_covered) - int(total_snps) * 100

        report_vals["total_snps"] = total_snps
        report_vals["snps_covered"] = str(snps_covered)
        report_vals["snps_not_covered"] = str(snps_not_covered)
        report_vals["snps_pct_covered"] = str(snps_pct_covered)
        report_vals["snps_pct_not_covered"] = str(snps_pct_not_covered)

        # add tables & plots to template
        html_string = self.build_report(
            html_template, total_stats, gene_stats, sub_thrshld_stats,
            snps_low_cov, snps_high_cov, fig, all_plots, summary_plot,
            report_vals
        )

        # write report
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        out_dir = os.path.join(bin_dir, "../output/")
        outfile = os.path.join(out_dir, args.output)

        file = open(outfile, 'w')
        file.write(html_string)
        file.close()


    def parse_args(self):
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

        args = parser.parse_args()

        if not args.sample_name:
            # sample name not given, use input file name
            args.sample_name = args.gene_stats.rsplit(".")[0]
            # input file includes full path, just get file name
            args.sample_name = args.sample_name.rsplit("/")[1]

        if not args.output:
            # output file name not given, using sample name
            args.output = args.sample_name + "_coverage_report.html"

        return args


def main():
    """
    Main function to generate coverage report
    """
    report = singleReport()

    args = report.parse_args()

    # read in files
    cov_stats, cov_summary, snp_df, raw_coverage,\
        html_template, build = report.load_files(
            args.threshold,
            args.exon_stats,
            args.gene_stats,
            args.raw_coverage,
            args.snps
        )

    if args.snps:
        # if SNP VCF(s) have been passed
        snps_low_cov, snps_high_cov = report.snp_coverage(
            snp_df, raw_coverage, args.threshold
        )
    else:
        snps_low_cov, snps_high_cov = None, None

    # generate summary plot
    summary_plot, cov_summary = report.summary_gene_plot(
        cov_summary, args.threshold
    )

    # get regions with low coverage
    low_raw_cov = report.low_coverage_regions(
        cov_stats, raw_coverage, args.threshold
    )

    # generate plot of sub optimal regions
    fig = report.low_exon_plot(low_raw_cov, args.threshold)

    # generate plots of each full gene
    all_plots = report.all_gene_plots(raw_coverage, args.threshold)

    # generate report
    report.generate_report(
        cov_stats, cov_summary, snps_low_cov, snps_high_cov,
        fig, all_plots, summary_plot, html_template, args, build
    )


if __name__ == "__main__":

    main()
