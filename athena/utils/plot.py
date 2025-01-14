"""Plotting related functions"""

from __future__ import annotations
from base64 import b64encode
from io import BytesIO
import math
import numpy as np
import pathlib
from timeit import default_timer as timer
from typing import List

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import polars as pl

from utils import log_handle
from .util_functions import call_in_parallel, format_timer


def to_html(plot: matplotlib.figure.Figure) -> str:
    """
    Converts matplotlib figure to HTML formatted string

    Parameters
    ----------
    plot : matplotlib.figure.Figure

    Returns
    -------
    str
        HTML formatted string of plot
    """
    buffer = BytesIO()
    plot.savefig(buffer, format="png", dpi=50, transparent=True)

    buffer.seek(0)
    graphic = b64encode(buffer.getvalue())
    buffer.close()

    data_uri = graphic.decode("utf-8")
    img_tag = (
        f"<img src=data:image/png;base64,{data_uri} style='max-width: "
        "100%; max-height: auto; object-fit: contain;' />"
    )

    return img_tag


def all_regions(
    coverage_data: pl.DataFrame, threshold: int
) -> List[List[str, str]]:
    """
    Generates the data for plotting all regions in the report.

    Data are returned as a list of all depths for each region of each
    transcript, with plotting happening on the fly using Plotly in the
    report. This is formatted as:


    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
    threshold : int
        threshold for low coverage

    Returns
    -------
    list
    """
    unique_transcripts = coverage_data["transcript"].unique().to_list()

    log_handle.debug(
        "Generating plot data for all %s transcripts", len(unique_transcripts)
    )

    start = timer()
    matplotlib.style.use("fast")

    coverage_data = coverage_data.with_columns(
        (pl.col("region_end") - pl.col("region_start")).alias("length")
    )

    plot_data = call_in_parallel(
        single_gene,
        unique_transcripts,
        progress=True,
        coverage_data=coverage_data,
        threshold=threshold,
    )

    log_handle.debug(
        "Generated all plot data in %s", format_timer(start=start, end=timer())
    )

    return plot_data


def single_gene(
    transcript: str, coverage_data: pl.DataFrame, threshold: int
) -> List[str, str]:
    """
    Generate the plot for a single gene.

    This generates subplots in a maximum of 20 plots per row, one per
    region (i.e. exon / intron). It uses matplotlib which is relatively
    slow and is the longest step of the generating the report.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of coverage data
    transcript : str
        transcript to generate plot for
    threshold : int
        threshold for low coverage

    Returns
    -------
    str
        gene and transcript of the plot
    str
        HTML string of the plot
    """
    transcript_filter = coverage_data.filter(
        pl.col("transcript") == transcript
    )

    total_regions = transcript_filter["region"].unique().shape[0]

    columns = min(total_regions, 20)
    rows = math.ceil(total_regions / 20)

    max_y_cols = [transcript_filter.select(pl.max("depth")).item(), threshold]

    if "normal_mean" in transcript_filter.columns:
        max_y_cols.append(transcript_filter.select("mean_+_std").max().item())

    max_y = max(max_y_cols) * 1.05

    fig = plt.figure(figsize=(40, rows * 4.5))

    grid = fig.add_gridspec(rows, columns, wspace=0)
    axs = grid.subplots(sharey=True)
    plt.subplots_adjust(
        left=0.02, right=0.98, top=0.78 + (rows * 0.05), bottom=0.05
    )

    gene = transcript_filter["gene"][0]
    if total_regions == 1:
        # handle single exon genes
        axs = np.array([axs])

    axs = axs.flatten()
    plt.setp(axs, xticks=[])

    for idx, region in enumerate(
        transcript_filter["region"].unique().to_list()
    ):
        region_filter = transcript_filter.filter(pl.col("region") == region)

        if "normal_mean" in region_filter.columns:
            # normal values have been provided => plot them
            axs[idx].plot(
                region_filter["position"].to_list(),
                region_filter["mean_-_std"].to_list(),
                color="#64e764",
                rasterized=True,
                markevery=None,
            )

            axs[idx].plot(
                region_filter["position"].to_list(),
                region_filter["mean_+_std"].to_list(),
                color="#64e764",
                rasterized=True,
                markevery=None,
            )

            axs[idx].fill_between(
                x=region_filter["position"].to_list(),
                y1=region_filter["mean_+_std"].to_list(),
                y2=region_filter["mean_-_std"].to_list(),
                color="#90ee90",
                rasterized=True,
            )

            axs[idx].plot(
                region_filter["position"].to_list(),
                region_filter["normal_mean"].to_list(),
                color="#64e764",
                rasterized=True,
                markevery=None,
            )

        if region_filter["depth"].unique().to_list() == [0]:
            axs[idx].plot(
                [0, 100],
                [threshold, threshold],
                color="red",
                linestyle="-",
                linewidth=2,
                rasterized=True,
                markevery=None,
            )
        else:
            axs[idx].plot(
                region_filter["position"].to_list(),
                region_filter["depth"].to_list(),
            )

            axs[idx].plot(
                [region_filter["position"][0], region_filter["position"][-1]],
                [threshold, threshold],
                color="red",
                linestyle="-",
                linewidth=1,
                rasterized=True,
                markevery=None,
            )

        fig.suptitle(
            f"{gene} ({transcript})",
            fontweight="bold",
            fontsize=14,
        )

        # remove y ticks & label for all but first plot of lines
        if idx == 0 or idx % 20 == 0:
            axs[idx].tick_params(axis="y", labelsize=12)
        else:
            axs[idx].yaxis.set_ticks_position("none")

        axs[idx].title.set_text(region)
        axs[idx].set_xlabel(f"{region_filter['length'][0]} bp", fontsize=13)
        plt.ylim(bottom=0, top=max_y)

    plot_html = to_html(fig)
    plt.cla()
    plt.clf()
    plt.close(fig)

    return [f"{gene}_{transcript}", plot_html]


def sub_threshold_regions(coverage_data: pl.DataFrame, threshold: int) -> str:
    """
    Generates the HTML formatted data of all exons with at least one base
    beneath given threshold depth for displaying in the low covered
    regions plots in the report.

    For each region a string is returned to add into the report with the
    title, start position and depth per base in the region. The positions
    are then generated from the length of depth bases when plotting to
    reduce the amount of data stored in the report.

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data
    threshold : int
        threshold beneath which regions are considered low coverage

    Returns
    -------
    str
        HTML formatted string representation of plot data
    """
    log_handle.debug("Generating data for low coverage regions plots")
    start = timer()

    # get the rows where depth for any position in the region under threshold
    low_coverage = (
        coverage_data.group_by("transcript", "region")
        .agg(
            pl.col("depth")
            .filter(pl.col("depth") < threshold)
            .alias("sub_threshold")
        )
        .filter(pl.col("sub_threshold") != [])
        .join(coverage_data, on=("transcript", "region"), how="left")
        .sort(by="gene", descending=False)
    )

    # format as a HTML string with transcript, first position, and depth
    low_coverage = (
        (
            low_coverage.group_by("transcript", "region", maintain_order=True)
            .agg(
                pl.first("gene"),
                pl.first("position"),
                pl.col("depth").str.join(",").alias("depths"),
            )
            .select(
                pl.format(
                    "<div class='sub_plot'>{} ({}) - {},{},{}</div>",
                    "gene",
                    "transcript",
                    "region",
                    "position",
                    "depths",
                ).alias("data")
            )
        )
        .get_column("data")
        .to_list()
    )

    low_coverage = ",".join([f'"{x}"' for x in low_coverage])

    log_handle.debug(
        "Generated low coverage regions plot data in %s",
        format_timer(start=start, end=timer()),
    )

    return low_coverage


def gene_summary(gene_coverage: pl.DataFrame, threshold: int) -> str:
    """
    Generate per gene coverage summary plot for top of report.

    Parameters
    ----------
    gene_coverage : pl.DataFrame
        DataFrame of summarised per transcript coverage values
    threshold : int
        low coverage cut off threshold value

    Returns
    -------
    str
        HTML formatted string of plot
    """
    threshold = f"{threshold}x"

    gene_coverage = (
        gene_coverage.select("gene", "transcript", threshold)
        .with_columns(
            pl.when(pl.col(threshold) < 90)
            .then(pl.lit("red"))
            .when(pl.col(threshold) < 99)
            .then(pl.lit("orange"))
            .otherwise(pl.lit("green"))
            .alias("colour")
        )
        .sort(by=[threshold, "gene"], descending=[True, False])
    )

    gene_coverage = gene_coverage.with_columns(
        pl.format("{} ({})", pl.col("gene"), pl.col("transcript")).alias(
            "label"
        )
    )

    summary_plot, axs = plt.subplots(figsize=(45, 20))
    total_genes = gene_coverage.height

    # limit the number of genes we plot for large panels for readability
    if total_genes > 100:
        if gene_coverage.filter(pl.col(threshold) < 100).height > 100:
            gene_coverage = gene_coverage.filter(pl.col(threshold) < 100)
        else:
            gene_coverage = gene_coverage.tail(100)

        total_omitted_genes = total_genes - gene_coverage.height
        axs.set_title(
            f"{total_omitted_genes} genes covered 100% at {threshold} were"
            " omitted from the plot due to the panel size",
            loc="left",
        )

    plt.bar(
        gene_coverage["label"],
        gene_coverage[threshold],
        color=gene_coverage["colour"],
    )

    # threshold lines
    plt.axhline(y=99, linestyle="--", color="#565656", alpha=0.6)
    plt.axhline(y=95, linestyle="--", color="#565656", alpha=0.6)

    plt.text(1.005, 0.94, "99%", transform=axs.transAxes)
    plt.text(1.005, 0.91, "95%", transform=axs.transAxes)

    # plot formatting
    axs.tick_params(labelsize=12, length=0)
    plt.xticks(
        rotation=55,
        color="#565656",
        ha="right",
        rotation_mode="anchor",
        weight="bold",
        fontsize=18,
    )

    # set appropriate margins
    axs.autoscale_view(scaley=True)
    if gene_coverage.height < 10:
        axs.margins(x=1 / gene_coverage.height**2)
    else:
        axs.margins(x=0.01)

    plt.legend(
        handles=[
            mpatches.Patch(color="green", label="100%"),
            mpatches.Patch(color="orange", label="90-99.99%"),
            mpatches.Patch(color="red", label="<90%"),
        ],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.4),
        fancybox=True,
        shadow=True,
        ncol=12,
        fontsize=18,
    )

    # set x tick label frequency to prevent overlap
    if gene_coverage.height > 125:
        x_tick_frequency = 2

        if gene_coverage.height > 250:
            x_tick_frequency = 3

        axs.set_xticks(axs.get_xticks()[::x_tick_frequency])
        axs.tick_params(axis="both", which="major", labelsize=10)
        plt.figtext(
            0.505,
            0.01,
            "Some gene labels are not shown due to high number of genes",
            ha="center",
            fontsize=12,
        )

    plt.xlabel("")
    plt.ylabel(f"% coverage ({threshold})", fontsize=11)
    plt.yticks(ticks=range(0, 110, 10), labels=range(0, 110, 10))

    axs.yaxis.grid(linewidth=0.5, color="grey", linestyle="-.")
    axs.set_axisbelow(True)

    plt.box(False)
    plt.tight_layout()

    return to_html(summary_plot)


def all_chromosomes(raw_coverage: pl.DataFrame) -> str:
    """
    Generates grid of plots of coverage per chromosome

    Parameters
    ----------
    raw_coverage : pl.DataFrame
        DataFrame of raw unbinned coverage data

    Returns
    -------
    str
        HTML string of generated plots
    """
    log_handle.debug("Generating full chromosome plots")
    start = timer()

    raw_coverage = raw_coverage.filter(pl.col("depth") > 50)

    chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    ncols, nrows = 24, 1

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(50, 10),
        sharey=True,
        constrained_layout=True,
    )

    fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, hspace=0, wspace=0)

    for plot_idx, chrom in enumerate(chroms):
        chrom_data = raw_coverage.filter(pl.col("chrom") == chrom)

        ax = axs[plot_idx]

        ax.scatter(
            x=chrom_data.get_column("depth_bin_start"),
            y=chrom_data.get_column("depth"),
            s=1,
        )

        ax.set_title(f"chr{chrom}", fontsize=24, fontstyle="italic")
        ax.xaxis.offsetText.set_fontsize(18)
        ax.xaxis.set_ticks_position("none")

        ax.set_yscale("log", base=2)

    axs[0].set_ylabel("Depth", fontsize=28)
    plt.setp(axs, xticks=[], yticks=[])

    log_handle.debug(
        "Generated full chromosome plots in %s",
        format_timer(start=start, end=timer()),
    )

    return to_html(fig)
