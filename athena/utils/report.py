from datetime import datetime
from pathlib import Path
from string import Template
from timeit import default_timer as timer
from typing import List

import polars as pl

from . import style
from .io import read_file, read_image
from .log import get_logger
from .util_functions import compress_and_encode, format_timer

log_handle = get_logger("athena")


def generate_summary_text(
    gene_df: pl.DataFrame,
    threshold: int,
    panel_coverage_pct: float,
    indication: str,
) -> str:
    """
    Generates the optional clinical report summary text

    Parameters
    ----------
    gene_coverage : pl.DataFrame
        DataFrame of summarised per transcript coverage values
    threshold : int
        threshold at which to calculate percent coverage
    panel_coverage_pct : float
        total percentage coverage of all regions
    clinical_indication : str
        clinical indication the panel relates to

    Returns
    -------
    str
        HTML formatted text of report summary
    """
    summary_text = ""

    if indication:
        summary_text += f"Gene panel(s): {indication}<br></br>"

    summary_text += "; ".join(
        gene_df.select("gene", "transcript")
        .select(
            pl.format("{} ({})", pl.col("gene"), pl.col("transcript")).alias(
                "text"
            )
        )
        .get_column("text")
        .to_list()
    )

    summary_text += (
        f"<br></br><b>Genes with coverage at {threshold}x less than 90%: </b>"
    )

    sub_90_genes = "; ".join(
        gene_df.filter(pl.col(f"{threshold}x") < 90)
        .select("gene", "transcript", f"{threshold}x")
        .select(
            pl.format(
                "{} ({}) {}%",
                pl.col("gene"),
                pl.col("transcript"),
                pl.col(f"{threshold}x").round(2),
            ).alias("text")
        )
        .get_column("text")
        .to_list()
    )

    if sub_90_genes:
        summary_text += sub_90_genes
    else:
        summary_text += "<b>None</b>"

    summary_text += (
        f"<br></br>{panel_coverage_pct} % of this panel was sequenced to a"
        f" depth of {threshold}x or greater.<br>"
    )

    return summary_text


def generate_panel_filters(filters: List[str]) -> str:
    """
    Generate HTML formatted filters for drop down menu for full
    gene plots in the report.

    Parameters
    ----------
    filters : list
        list of panel -> gene filters from input

    Returns
    -------
    str
        HTML formatted option list to pass to the report

    Raises
    ------
    ValueError
        Raised when all filter strings do not contain exactly 1 colon
    """
    if not all([x.count(":") == 1 for x in filters]):
        raise ValueError("Invalid filter string(s) provided")

    return "".join(
        f'<option value="{x.split(":")[1]}">{x.split(":")[0]}</option>'
        for x in filters
    )


def get_sub_threshold_regions(
    region_df: pl.DataFrame, threshold: int
) -> pl.DataFrame:
    """
    Returns the regions where the target threshold is not covered to 100%

    Parameters
    ----------
    region_df : pl.DataFrame
        DataFrame of summarised per region coverage

    Returns
    -------
    pl.DataFrame
        DataFrame of regions under 100% at given threshold

    Raises
    ------
    ValueError
        Raised when provided threshold not in DataFrame columns
    """
    threshold = f"{threshold}x"

    if threshold not in region_df.columns:
        raise ValueError(
            f"provided threshold column '{threshold}' not in DataFrame"
        )

    return region_df.filter(pl.col(threshold) < 100)


def get_total_unique_regions(gene_df: pl.DataFrame) -> tuple((int, int)):
    """
    Get total unique number of genes and transcripts

    Parameters
    ----------
    gene_df : pl.DataFrame
        DataFrame of per gene coverage data

    Returns
    -------
    int
        Total number of unique genes
    int
        Total number of unique transcripts
    """
    return (
        gene_df.select("gene").unique().height,
        gene_df.select("transcript").unique().height,
    )


def get_total_fully_covered_genes(
    gene_df: pl.DataFrame, threshold: int
) -> int:
    """
    Get the total number of genes with 100% coverage at given threshold

    Parameters
    ----------
    gene_df : pl.DataFrame
        DataFrame of per gene coverage data
    threshold : int
        Threshold for low coverage

    Returns
    -------
    int
        Total number of genes covered at 100%

    Raises
    ------
    ValueError
        Raised when provided threshold not in DataFrame columns
    """
    threshold = f"{threshold}x"

    if threshold not in gene_df.columns:
        raise ValueError(
            f"provided threshold column '{threshold}' not in DataFrame"
        )

    return (
        gene_df.group_by("gene")
        .agg(pl.col(threshold))
        .filter(pl.col(threshold) == [100.0])
        .height
    )


def get_total_sub_threshold_genes_and_regions(
    region_df: pl.DataFrame, threshold: int
) -> tuple((int, int)):
    """
    Get the total number of genes and exons that are under the given threshold

    Parameters
    ----------
    region_df : pl.DataFrame
        DataFrame of summarised per region coverage
    threshold : int
        Threshold for low coverage

    Returns
    -------
    int
        Total number of genes under 100% coverage at threshold
    int
        Total number of regions under 100% coverage at threshold

    Raises
    ------
    ValueError
        Raised when provided threshold not in DataFrame columns
    """
    threshold = f"{threshold}x"

    if threshold not in region_df.columns:
        raise ValueError(
            f"provided threshold column '{threshold}' not in DataFrame"
        )

    sub_threshold_genes = (
        region_df.group_by("gene")
        .agg(pl.col(threshold).unique())
        .filter(pl.col(threshold) != [100.0])
        .height
    )

    sub_threshold_regions = (
        region_df.group_by("gene", "region")
        .agg(pl.col(threshold).unique())
        .filter(pl.col(threshold) != [100.0])
        .height
    )

    return sub_threshold_genes, sub_threshold_regions


def populate_template(
    summary_text: str,
    gene_df: pl.DataFrame,
    region_df: pl.DataFrame,
    sub_threshold_plot_data: list,
    all_region_plots: list,
    summary_plot: str,
    chromosome_plots: str,
    threshold: int,
    sample: str,
    build: str,
    panel: str,
    panel_coverage_pct: float,
    panel_filters: list,
    version: str,
) -> str:
    """
    Populate the HTML template with all data for the report

    Parameters
    ----------
    summary_text : str
        clinical report summary text
    gene_df : pl.DataFrame
        DataFrame of summarised per transcript coverage values
    region_df : pl.DataFrame
        DataFrame of summarised per region coverage values
    sub_threshold_plot_data : list
        Data for sub threshold plots
    all_region_plots : list
        HTML string plots of all regions
    summary_plot : str
        Summary plot for top level of report
    chromosome_plots : str
        Per chromosome plots
    threshold : int
        Threshold for defining low coverage
    sample : str
        name of sample report is generated for
    build : str
        Reference build used for data
    panel : str
        Panel the genes belong to
    panel_coverage_pct : float
        Percentage covered of the panel
    panel_filters : list
        list of colon separated panel names to gene list, used for
        preset filters
    version : str
        current version of Athena to write to footer

    Returns
    -------
    str
        String representation of the populated report
    """
    log_handle.debug("Populating report template")
    start = timer()

    template_path = (
        Path(__file__)
        .absolute()
        .parent.parent.joinpath("data/templates/template.html")
    )
    template_contents = read_file(file=template_path)
    template = Template(template_contents)

    logo_path = (
        Path(__file__)
        .absolute()
        .parent.parent.joinpath("data/images/logo_small.png")
    )
    logo = read_image(file=logo_path)

    total_genes, total_transcripts = get_total_unique_regions(gene_df=gene_df)
    total_covered_genes = get_total_fully_covered_genes(
        gene_df=gene_df, threshold=threshold
    )
    total_sub_threshold_genes, total_sub_threshold_regions = (
        get_total_sub_threshold_genes_and_regions(
            region_df=region_df, threshold=threshold
        )
    )

    sub_threshold_df = get_sub_threshold_regions(
        region_df=region_df, threshold=threshold
    )

    sub_threshold_data, sub_threshold_columns = style.dataframe_for_html(
        coverage_df=sub_threshold_df, sort_by=("Transcript", "Region")
    )
    region_data, region_columns = style.dataframe_for_html(
        coverage_df=region_df, sort_by=("Transcript", "Region")
    )
    gene_data, gene_columns = style.dataframe_for_html(
        coverage_df=gene_df, sort_by=("Transcript",)
    )

    if panel_filters and panel_filters != [""]:
        panel_filters = generate_panel_filters(panel_filters)

    sub_threshold_data = compress_and_encode(sub_threshold_data)
    gene_data = compress_and_encode(gene_data)
    region_data = compress_and_encode(region_data)
    sub_threshold_plot_data = compress_and_encode(sub_threshold_plot_data)

    report_data = template.safe_substitute(
        name=sample,
        threshold=threshold,
        summary_text=summary_text,
        panel=panel,
        panel_pct_coverage=panel_coverage_pct,
        total_genes=total_genes,
        total_transcripts=total_transcripts,
        total_sub_threshold_regions=total_sub_threshold_regions,
        total_sub_threshold_genes=total_sub_threshold_genes,
        fully_covered_genes=total_covered_genes,
        low_exon_columns=sub_threshold_columns,
        sub_threshold_stats=sub_threshold_data,
        gene_table_headings=gene_columns,
        gene_stats=gene_data,
        region_table_headings=region_columns,
        region_stats=region_data,
        summary_plot=summary_plot,
        sub_threshold_plots=sub_threshold_plot_data,
        all_region_plots=all_region_plots,
        coverage_per_chromosome_fig=chromosome_plots,
        panel_filters=panel_filters,
        date=datetime.today().strftime("%H:%M %Y-%m-%d"),
        build=build,
        version=version,
        logo=logo,
    )

    log_handle.debug(
        "Populated report in %s", format_timer(start=start, end=timer())
    )

    return report_data
