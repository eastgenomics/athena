"""Main entrypoint to control all running of Athena"""

import argparse
from timeit import default_timer as timer

from utils import calculate
from utils import log_handle
from utils.annotate import call_bedtools_intersect
from utils.arguments import parse_args
from utils.io import (
    read_annotated_bed,
    read_hsmetrics,
    read_sample_files,
    write_file,
    write_multi_sample_coverage,
)
from utils import plot
from utils.report import generate_summary_text, populate_template
from utils.util_functions import (
    call_in_parallel,
    format_timer,
    pair_up_sample_files,
    strip_html_markup,
    unbin,
)

import polars as pl

pl.enable_string_cache()


def generate_report(args: argparse.Namespace) -> None:
    """
    Call all methods for calculating coverage and generating report

    Parameters
    ----------
    args : argparse.Namespace
        Command line argument Namespace object
    """
    print("Beginning generating coverage stats and coverage report")
    start = timer()

    if args.annotated_bed:
        annotated_bed_file = args.annotated_bed
    else:
        annotated_bed_file = call_bedtools_intersect(
            regions=args.regions,
            coverage=args.coverage,
            overwrite=args.force,
            build=args.build,
        )

    per_base_df = read_annotated_bed(annotated_bed=annotated_bed_file)
    per_base_df = unbin(coverage_data=per_base_df)

    # generate stats
    gene_df, exon_df = calculate.region_coverage(
        coverage_data=per_base_df, thresholds=args.thresholds
    )
    panel_coverage_pct = calculate.total_pct_coverage(
        coverage_data=per_base_df, threshold=args.minimum
    )

    if args.normal_coverage:

        from utils.constants import NORM_VALUE

        hsmetrics_df = read_hsmetrics(hsmetrics_file=args.hsmetrics)
        normal_coverage_df = pl.read_csv(
            source=args.normal_coverage,
            separator="\t",
            comment_prefix="#",
            schema={
                "chrom": pl.Categorical,
                "position": pl.UInt32,
                "mean": pl.Float64,
                "std": pl.Float64,
            },
        )

        sample_bases = hsmetrics_df.select(
            pl.col("ON_TARGET_BASES").cast(pl.Int32)
            * pl.col("PCT_USABLE_BASES_ON_TARGET").cast(pl.Float64)
        ).item()

        norm_factor = sample_bases / NORM_VALUE

        normal_coverage_df = normal_coverage_df.with_columns(
            (pl.col("mean") * norm_factor).alias("normal_mean"),
            (pl.col("std") * norm_factor).alias("normal_std"),
        ).drop("mean", "std")

        normal_coverage_df = normal_coverage_df.with_columns(
            (pl.col("normal_mean") - (pl.col("normal_std")) * 3).alias(
                "mean_-_std"
            ),
            (pl.col("normal_mean") + (pl.col("normal_std")) * 3).alias(
                "mean_+_std"
            ),
        )

        per_base_df = per_base_df.join(
            normal_coverage_df, how="left", on=["chrom", "position"]
        )

        # print(per_base_df)
        # print(per_base_df.columns)
        # exit()

    # generate plots
    # sub_threshold_plot_data = plot.sub_threshold_regions(
    #     coverage_data=per_base_df, threshold=args.minimum
    # )
    sub_threshold_plot_data = ""
    all_region_plots = ""
    all_region_plots = plot.all_regions(
        coverage_data=per_base_df, threshold=args.minimum
    )
    summary_plot = plot.gene_summary(
        gene_coverage=gene_df, threshold=args.minimum
    )

    if args.summary:
        summary_text = generate_summary_text(
            gene_df=gene_df,
            threshold=args.minimum,
            panel_coverage_pct=panel_coverage_pct,
            indication=args.clinical_indication,
        )
    else:
        summary_text = ""

    populated_report = populate_template(
        summary_text=summary_text,
        per_base_df=exon_df,
        gene_df=gene_df,
        region_df=exon_df,
        sub_threshold_plot_data=sub_threshold_plot_data,
        all_region_plots=all_region_plots,
        summary_plot=summary_plot,
        chromosome_plot=None,
        threshold=args.minimum,
        sample=args.output,
        build=args.build,
        panel=args.panel,
        panel_coverage_pct=panel_coverage_pct,
        panel_filters=args.panel_filters,
    )

    output_file = f"{args.output}_coverage_report.html"
    write_file(file=output_file, contents=populated_report)

    per_base_df.write_csv(
        file=f"{args.output}.coverage.bed", include_header=True, separator="\t"
    )
    gene_df.write_csv(
        file=f"{args.output}.gene_coverage.tsv",
        include_header=True,
        separator="\t",
    )
    exon_df.write_csv(
        file=f"{args.output}.region_coverage.tsv",
        include_header=True,
        separator="\t",
    )

    if args.summary_file:
        write_file(
            file=f"{args.output}_summary.txt",
            contents=strip_html_markup(summary_text),
        )

    print(
        "Completed all steps in"
        f" {format_timer(start=start, end=timer())}. Report written to"
        f" {output_file}",
    )


def generate_multi_sample_coverage(args: argparse.Namespace) -> None:
    """
    Calculates mean per base coverage and standard deviation from the
    mean for all provided samples. This is to generate a file to define
    'normal' coverage for adding context to whole gene plots in the
    output report.

    Parameters
    ----------
    args : argparse.Namespace
        Command line argument Namespace object
    """
    log_handle.info(
        "Calculating multi sample coverage from %s samples", len(args.coverage)
    )

    annotated_beds = call_in_parallel(
        call_bedtools_intersect,
        items=args.coverage,
        progress=True,
        regions=args.regions,
        build=args.build,
        overwrite=True,
    )

    sample_files = pair_up_sample_files(
        hsmetrics_files=args.hsmetrics, coverage_files=annotated_beds
    )

    sample_dfs = call_in_parallel(read_sample_files, sample_files.values())

    normalised_coverage_df = calculate.multi_sample_mean_and_std_dev(
        sample_dfs=sample_dfs
    )

    write_multi_sample_coverage(
        filename=f"{args.output}.tsv", coverage_df=normalised_coverage_df
    )

    log_handle.info("Completed calculating multi sample coverage.")


def main():
    args = parse_args()

    if args.verbose:
        log_handle.setLevel("DEBUG")

    if args.mode == "report":
        generate_report(args=args)
    else:
        generate_multi_sample_coverage(args=args)


if __name__ == "__main__":
    main()
