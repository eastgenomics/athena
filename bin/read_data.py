
"""Functions to import data and get required attributes"""
from pathlib import Path
import os
import sys

import pandas as pd


def read_exon_stats(exon_stats):
    """
    Read exon stats file from coverage_single_stats into df

    Args:
        - exon_stats (file): file with per exon coverage stats
    Returns:
        - cov_stats (df): df of exon_stats file
    """
    with open(exon_stats.name) as exon_file:
        dtypes = {
            "chrom": str, "exon_start": int, "exon_end": int, "gene": str,
            "tx": str, "exon": int, "min": int, "mean": float, "max": int,
            r'[0-9]*x': float, "exon_len": int
        }

        cov_stats = pd.read_csv(
            exon_file, sep="\t", comment='#', dtype=dtypes
        )

        # strip chr from chrom in cases of diff. formatted bed
        cov_stats["chrom"] = cov_stats["chrom"].apply(
            lambda x: str(x).replace("chr", "")
        )

    return cov_stats

def read_gene_stats(gene_stats):
    """
    Read gene stats file from coverage_single_stats into df

    Args:
        - gene_stats (file): file with per gene coverage stats
    Returns:
        - cov_summary (df): df of gene_stats file
    """
    with open(gene_stats) as gene_file:
        dtypes = {
            "gene": str, "tx": str, "min": int,
            "mean": float, "max": int, r'[0-9]*x': float
        }

        cov_summary = pd.read_csv(
            gene_file, sep="\t", comment='#', dtype=dtypes
        )

    return cov_summary


def read_raw_coverage(raw_coverage):
    """
    Read in raw coverage data (annotated bed file) from single stats

    Args:
        - raw_coverage (file): file of annotated bed file
    Returns:
        - raw_coverage (df): df of raw coverage
    """
    column = [
        "chrom", "exon_start", "exon_end", "gene", "tx", "exon",
        "cov_start", "cov_end", "cov"
    ]

    dtypes = {
        "chrom": str, "exon_start": int, "exon_end": int, "gene": str,
        "tx": str, "exon": int, "cov_start": int, "cov_end": int,
        "cov": int
    }

    # read in raw coverage stats file
    with open(raw_coverage) as raw_file:
        raw_coverage = pd.read_csv(
            raw_file, sep="\t", names=column, dtype=dtypes
        )
        # strip chr from chrom in cases of diff. formatted bed
        raw_coverage["chrom"] = raw_coverage["chrom"].apply(
            lambda x: str(x).replace("chr", "")
        )

    return raw_coverage


def read_bootstrap():
    """
    Read in bootstrap for styling report

    Args: None
    Returns:
        - bootstrap (str): str of bootstrap file to store in report
    """
    bs = str(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "../data/static/css/bootstrap.min.css"
    ))
    with open(bs) as bs:
        bootstrap = bs.read()

    return bootstrap


def read_template():
    """
    Read in HTML template for report

    Args: None

    Returns:
        - html_template (str): report template
    """
    bin_dir = os.path.dirname(os.path.abspath(__file__))
    template_dir = os.path.join(bin_dir, "../data/templates/")
    single_template = os.path.join(template_dir, "single_template.html")

    with open(single_template, 'r') as template:
        html_template = template.read()

    return html_template


def get_low_coverage_regions(cov_stats, raw_coverage, threshold):
    """
    Get regions where coverage at given threshold is <100% for
    generating low coverage plots

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
        'chrom': str,
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


def get_build_and_stats(gene_stats):
    """
    Get flagstats (if present) and reference build number used for
    alignment from header of gene_stats file.
    Args:
        - gene_stats (file): gene stats file from single_stats
    Returns:
        - flagstat (dict): flagstat values
        - build (str): HTML formatted string of reference build
    """
    flagstat = {}
    # read in flagstat and build from header of gene stats file
    with open(gene_stats) as gene_file:
        for ln in gene_file:
            if ln.startswith("#"):
                if "build" in ln:
                    # get build number
                    reference = ln.split(":")[1]
                    # add build to string to display
                    if "37" in reference:
                        build = "<li>Reference build used for aligment<b>\
                            {}</b></li>".format(reference)
                    if "38" in build:
                        build = "<li>Reference build used for aligment<b>\
                            {}</b></li>".format(reference)
                else:
                    # read in flagstat from header
                    key = ln.split(":")[0].strip("#")
                    val = ln.split(":")[1]
                    flagstat[key] = val

    if "build" not in locals():
        # build no. not included in gene_stats file
        build = ""

    return flagstat, build


def get_panel_name(panel):
    """
    Get panel name from panel bed file for displaying in the report
    Args:
        - panel (file): panel bed file
    Returns:
        - panel (str): HTML formatted str for displaying in report
    """
    if panel is not None:
        panel_name = Path(panel).stem

        # format according to output of
        # https://github.com/eastgenomics/eggd_generate_bed
        panel_name = [x.strip("_") for x in panel_name.split("&&") if x]
        panel_name = [
            x.strip("_b37").strip("_b38") for x in panel_name if x
        ]
        panel_name = [x.replace("_", " ") for x in panel_name if x]
        panel_name = ",&nbsp".join(panel_name)
        panel = "<li>Panel(s) / gene(s) included in report: <b>{}</b>\
            </li>".format(panel_name)
    else:
        panel = ""

    return panel


def get_snp_vcfs(snp_vcfs):
    """
    Get names of SNP VCFs (if used) to display in report summary

    Args:
        - snp_vcfs (array): array of file names of SNP VCFs

    Returns:
        - vcfs (str): HTML formatted string of SNP VCFs
    """
    if snp_vcfs:
        # get names of SNP vcfs used to display in report
        vcfs = ", ".join([Path(x).stem for x in snp_vcfs])
        vcfs = "<br>VCF(s) of known variants included in report: <b>{}</b>\
            </br>".format(vcfs)
    else:
        vcfs = ""

    return vcfs


def get_athena_ver():
    """
    Attempt to get version of Athena from dir name to display in
    report footer, will only work for zip/tar

    Args: None
    Returns:
        - version (str): version of Athena
    """
    bin_dir = os.path.dirname(os.path.abspath(__file__))

    try:
        path = str(os.path.join(bin_dir, "../")).split("/")
        version = [s for s in path if "athena" in s][0].split("-")[1]
        version = "({})".format(version)
    except Exception:
        print("Error getting version from dir name, continuing.")
        version = ""
        pass

    return version


def check_threshold(threshold, cov_stats, cov_summary):
    """
    Check the given low coverage threshold is one of the threshold
    columns in both stats files

    Args:
        - threshold (int): given low threshold value
        - cov_stats (df): df of exon_stats file
        - cov_summary (df): df of gene_stats file
    Returns:
        - threshold (str): given low threshold value
    """
    if "x" not in str(threshold):
        threshold = str(threshold) + "x"

    if threshold not in list(cov_stats) and\
            threshold not in list(cov_summary):
        print("""--threshold must be one of the gene and exon
            stats coverage thresholds. Exiting now.""")
        sys.exit()

    return str(threshold)
