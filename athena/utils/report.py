from datetime import datetime
from pathlib import Path
from string import Template
from timeit import default_timer as timer

import polars as pl

# from athena import VERSION
from .io import read_file, read_image
from .util_functions import format_timer
from utils import log_handle
from . import style


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
    indication : str
        clinical indication the panel relates to

    Returns
    -------
    str
        HTML formatted text of report summary
    """
    summary_text = ""

    return summary_text


def get_sub_threshold_regions(
    region_df: pl.DataFrame, threshold: int
) -> pl.DataFrame:
    """
    Returns the regions where the target threshold is not covered to 100%

    Parameters
    ----------
    region_df : pl.DataFrame
        _description_

    Returns
    -------
    pl.DataFrame
        _description_
    """
    return region_df.filter(pl.col(f"{threshold}x") < 100)


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
    _summary_

    Parameters
    ----------
    gene_df : pl.DataFrame
        DataFrame of per gene coverage data
    threshold : int
        threshold for low coverage

    Returns
    -------
    int
        _description_
    """
    return (
        gene_df.filter(pl.col(f"{threshold}x") == 100)
        .select("gene")
        .unique()
        .height
    )


def get_total_sub_threshold_regions(
    region_df: pl.DataFrame, threshold: int
) -> tuple((int, int)):
    """
    Get the total number of genes and exons that are under the given threshold

    Parameters
    ----------
    region_df : _type_
        _description_
    int : _type_
        _description_

    Returns
    -------
    int
        _
    int
        _
    """
    sub_threshold_genes = (
        region_df.filter(pl.col(f"{threshold}x") < 100)
        .select("gene")
        .unique()
        .height
    )
    sub_threshold_regions = (
        region_df.filter(pl.col(f"{threshold}x") < 100)
        .select("gene", "region")
        .unique()
        .height
    )

    return sub_threshold_genes, sub_threshold_regions


def populate_template(
    per_base_df: pl.DataFrame,
    gene_df: pl.DataFrame,
    region_df: pl.DataFrame,
    sub_threshold_plot_data: list,
    all_regions_plot_data: list,
    summary_plot: str,
    chromosome_plot: str,
    threshold: int,
    sample: str,
    build: str,
    panel: str,
    panel_coverage_pct: float,
):
    """
    Populate the HTML template with all data for the report

    Parameters
    ----------
    per_base_df : pl.DataFrame
        _description_
    gene_df : pl.DataFrame
        DataFrame of summarised per transcript coverage values
    region_df : pl.DataFrame
        _description_
    sub_threshold_plot_data : list
        _description_
    all_regions_plot_data : list
        _description_
    summary_plot : str
        _description_
    chromosome_plot : str
        _description_
    threshold : int
        _desctipion_
    sample : str
        name of sample report is generated for
    build : str
        _description_
    panel : str
        _description_
    panel_coverage_pct : float
        _description_
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
        .parent.parent.joinpath("data/images/logo.png")
    )
    logo = read_image(file=logo_path)

    total_genes, total_transcripts = get_total_unique_regions(gene_df=gene_df)
    total_covered_genes = get_total_fully_covered_genes(
        gene_df=gene_df, threshold=threshold
    )
    total_sub_threshold_genes, total_sub_threshold_regions = (
        get_total_sub_threshold_regions(
            region_df=region_df, threshold=threshold
        )
    )

    summary_text = generate_summary_text(
        gene_df=gene_df,
        threshold=threshold,
        panel_coverage_pct=panel_coverage_pct,
        indication=None,
    )

    sub_threshold_df = get_sub_threshold_regions(
        region_df=region_df, threshold=threshold
    )

    sub_threshold_data, sub_threshold_columns = style.dataframe_for_html(
        coverage_df=sub_threshold_df, sort_by=("Transcript", "Region")
    )
    region_df, region_df_columns = style.dataframe_for_html(
        coverage_df=region_df, sort_by=("Transcript", "Region")
    )
    gene_df, gene_df_columns = style.dataframe_for_html(
        coverage_df=gene_df, sort_by=("Transcript",)
    )

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
        gene_table_headings=gene_df_columns,
        gene_stats=gene_df,
        region_table_headings=region_df_columns,
        region_stats=region_df,
        summary_plot=summary_plot,
        sub_threshold_plots=sub_threshold_plot_data,
        all_plots=all_regions_plot_data,
        # coverage_per_chromosome_fig=coverage_per_chromosome_fig,
        panel_filters=None,
        hide_filter=False,
        hide_plots=False,
        date=datetime.today().strftime("%Y-%m-%d"),
        build=build,
        version=1,
        logo=logo,
    )

    with open("report.html", "w") as fh:
        fh.write(report_data)

    log_handle.debug(
        "Populated report in %s", format_timer(start=start, end=timer())
    )
