"""Plotting related functions"""

from __future__ import annotations
from base64 import b64encode
from functools import reduce
from io import BytesIO
import pathlib
from timeit import default_timer as timer

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import polars as pl

from .util_functions import format_timer
from utils import log_handle


def to_html(plot: matplotlib.figure.Figure):
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
    plot.savefig(buffer, format="png", dpi=65, transparent=True)

    buffer.seek(0)
    graphic = b64encode(buffer.getvalue())
    buffer.close()

    data_uri = graphic.decode("utf-8")
    img_tag = (
        f"<img src=data:image/png;base64,{data_uri} style='max-width: "
        "100%; max-height: auto; object-fit: contain;' />"
    )

    return img_tag


def all_regions(coverage_data: pl.DataFrame) -> list(dict):
    """
    Generates the data for plotting all regions in the report.

    Data are returned as a list of all depths for each region of each
    transcript, with plotting happening on the fly using Plotly in the
    report. This is formatted as:

    {
        'NM_000123.4': [
            (1, [36, 39, 47, 51, 43, 52, ...]),
            (2, [33, 32, 38, 44, 40, 41, ...]),
            ...
        ],
        'NM_000567.8': ...
    }

    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame of per base coverage data

    Returns
    -------
    list
        list of dicts of data per transcript
    """
    log_handle.debug("Generating plot data for all regions")

    plot_data = (
        coverage_data.group_by("transcript")
        .agg(
            pl.max("depth").alias("max_depth"),
        )
        .join(coverage_data, on="transcript", how="left")
    )

    plot_data = (
        plot_data.with_columns(
            pl.format("{} ({})", "gene", "transcript").alias("title"),
            (pl.col("region_end") - pl.col("region_start")).alias("length"),
        )
        .group_by("title", "region")
        .agg(
            pl.first("max_depth"),
            pl.first("position").alias("start"),
            pl.col("depth").alias("depths"),
            pl.first("length"),
        )
        .select(
            "title",
            "region",
            "start",
            "length",
            "max_depth",
            "depths",
        )
    ).rows_by_key(key="title", named=True)

    plot_data = [
        f"<div class='gene_sub_plot' title='{title}'>{data}</div>"
        for title, data in plot_data.items()
    ] * 10

    # print(plot_data)

    return plot_data


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
                pl.concat_str(
                    pl.first("transcript"), pl.first("region"), separator=" "
                ).alias("title"),
                pl.first("position"),
                pl.col("depth").str.join(",").alias("depths"),
            )
            .select(
                pl.format(
                    "<div class='sub_plot'>{},{},{}</div>",
                    "title",
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
        "Generated plot data in %s", format_timer(start=start, end=timer())
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

    summary_plot, axs = plt.subplots(figsize=(25, 10))
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
    axs.tick_params(labelsize=8, length=0)
    plt.xticks(
        rotation=55,
        color="#565656",
        ha="right",
        rotation_mode="anchor",
        weight="bold",
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
        fontsize=14,
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

    axs.tick_params(axis="both", which="major", labelsize=10)

    plt.xlabel("")
    plt.ylabel(f"% coverage ({threshold})", fontsize=11)
    plt.yticks(ticks=range(0, 110, 10), labels=range(0, 110, 10))

    axs.yaxis.grid(linewidth=0.5, color="grey", linestyle="-.")
    axs.set_axisbelow(True)

    plt.box(False)
    plt.tight_layout()

    return to_html(summary_plot)


def all_chromosomes(coverage_file: pathlib.Path) -> str:
    """
    Generates plots of depth across all chromosomes.

    Requires the full raw coverage data output from samtools / mosdepth.

    Parameters
    ----------
    coverage_file : pathlib.Path
        _description_

    Returns
    -------
    str
        _description_
    """
    pass
