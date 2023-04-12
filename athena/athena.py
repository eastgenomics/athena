"""
Main script to generate coverage statistics and HTML report for given sample

Jethro Rainford
03/07/2021
"""
import argparse
from functools import partial
import multiprocessing
import numpy as np
import os
import pandas as pd
from pathlib import Path
import pybedtools as bedtools

from utils.load import loadData
from utils.annotate import annotateBed
# from coverage_stats_single import singleCoverage
# from coverage_report_single import calculateValues, generatePlots, \
# generateReport, styleTables, load_files


def annotate_bed(args):
    """
    Calls functions to annotate panel bed file with transcript info and coverage

    Parameters
    ----------

    Returns
    -------

    Outputs
    -------
    """
    if not args.output_name:
        # output name not defined, use sample identifier from coverage file
        args.output_name = Path(args.coverage_file).name.split('_')[0]

    # set dir for writing to
    bin_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(bin_dir, "../../output/")
    outfile_name = f"{args.output_name}_annotated.bed"
    outfile = os.path.join(out_dir, outfile_name)

    # add transcript info
    bed_w_transcript = annotateBed().add_transcript_info(
        panel_bed_df, transcript_info_df
    )

    # add coverage
    if args.chunk_size:
        # per-base coverage split to multiple dfs to limit memory usage
        bed_w_coverage = annotateBed().add_coverage(
            bed_w_transcript, pb_coverage_df, chunks=True
        )
    else:
        bed_w_coverage = annotateBed().add_coverage(
            bed_w_transcript, pb_coverage_df, chunks=False
        )

    # sense check generated file isn't empty, should be caught earlier
    assert len(bed_w_coverage.index) > 0, (
        'An error has occured: annotated bed file is empty. This is likely ',
        'due to an error in regions defined in bed file (i.e. different ',
        'transcripts to those in the transcripts file). Start debugging by ',
        'intersecting files manually...'
    )

    return bed_w_coverage


def call_generate_stats(args, data, thresholds, flagstat, build):
    """
    Calls functions to generate per exon and per gene stats

    Args:

    Returns:

    Outputs:
    """
    stats = singleCoverage()

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
            stats_w_arg = partial(stats.cov_stats, thresholds=thresholds)
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
                stats.summary_stats, map(
                    lambda e: (e, thresholds), split_stats_dfs
                )), ignore_index=True)

    cov_summary = cov_summary.sort_values(['gene'])
    cov_summary = cov_summary.drop(columns=["exon"])

    # write tables to output files
    stats.write_outfiles(
        cov_stats, cov_summary, args.outfile, flagstat, build
    )

    return cov_stats, cov_summary


def call_generate_report(args, load):
    """
    Calls functions to generate HTML report

    Args:

    Returns:

    Outputs:
    """
    calculate = calculateValues(args.threshold)
    plots = generatePlots(args.threshold)
    report = generateReport(args.threshold)
    styling = styleTables()

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
        genes = low_raw_cov.drop_duplicates(["gene"])["gene"].values.tolist()
        print(f"Plots for {len(genes)} to generate")

        # split gene list equally for seperate processes
        gene_array = np.array_split(np.array(genes), num_cores)

        # split df into seperate dfs by genes in each list
        split_dfs = np.asanyarray(
            [low_raw_cov[low_raw_cov["gene"].isin(x)] for x in gene_array],
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

    # generate report
    report.generate_report(
        cov_stats, cov_summary, snps_low_cov, snps_high_cov, snps_no_cov, fig,
        all_plots, summary_plot, html_template, args, build, panel, vcfs,
        panel_pct_coverage, bootstrap, version, summary_text
    )


def parse_args():
    """
    Parse cmd line arguments

    Args:

    Returns:
    """
    parser = argparse.ArgumentParser(
        description='Generate coverage stats and HTML report.'
    )

    parser.add_argument(
        '--panel_bed', '-p',
        help='panel bed file'
    )
    parser.add_argument(
        '--transcript_file', '-t',
        help='file with gene and exon information'
    )
    parser.add_argument(
        '--coverage_file', '-c',
        help='per base coverage data file'
    )
    parser.add_argument(
        '--chunk_size', '-s', type=int,
        help='number lines to read per-base coverage file in one go'
    )
    parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
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
        '-n', '--sample_name', nargs='?',
        help="Name of sample to display in report, if not\
            specified this will be the prefix of the\
            gene_stats input file.",
        required=False
    )
    parser.add_argument(
        '-o', '--output', nargs='?',
        help='name preifx for output file, if none will use coverage file',
        required=False
    )
    parser.add_argument(
        '--panel', nargs='?',
        help='(Optional) Panel bed file used from annotation, if passed\
        name of file will be displayed in report to show what\
        panel(s) / gene(s) were included.',
        required=False
    )
    parser.add_argument(
        '--thresholds', nargs='*',
        default=[10, 20, 30, 50, 100],
        help='threshold values to calculate coverage for as comma\
            seperated integers (default: 10, 20, 30, 50, 100).'
    )
    parser.add_argument(
        '--sub_threshold', nargs='?',
        default=20, type=int,
        help="threshold to define low coverage (int). Must be one defined \
            by --thresholds",
        required=False
    )

    # sub commands for running seperate parts
    subparsers = parser.add_subparsers(
        dest='development', help='sub-commands for development and testing'
    )

    # parser to annotate panel bed file
    annotate_parser = subparsers.add_parser(
        'annotate', help='annotate bed file with transcript & coverage info, \
            requires panel bed, transcript and coverage args'
    )
    annotate_parser.add_argument(
        '--panel_bed', '-p', required=True,
        help='panel bed file'
    )
    annotate_parser.add_argument(
        '--transcript_file', '-t', required=True,
        help='file with gene and exon information'
    )
    annotate_parser.add_argument(
        '--coverage_file', '-c', required=True,
        help='per base coverage data file'
    )
    annotate_parser.add_argument(
        '--chunk_size', '-s', type=int,
        help='number lines to read per-base coverage file in one go'
    )
    annotate_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )

    # parser to generate coverage stats
    stats_parser = subparsers.add_parser(
        'stats', help='generate per gene and per exon coverage stats'
    )
    stats_parser.add_argument(
        'annotated_bed',
        help='annotated bed file from annotated_bed'
    )
    stats_parser.add_argument(
        '--flagstat', nargs='?',
        help='Optional flagstat file, needed for generating run stats.'
    )
    stats_parser.add_argument(
        '--build', nargs='?',
        help='Optional text file with build number used for alignment.'
    )
    stats_parser.add_argument(
        '--outfile', nargs='?', help='Output file name prefix, if not\
        given the input file name will be used as the name prefix.',
        type=str
    )
    stats_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )


    # parser to generate report
    report_parser = subparsers.add_parser(
        'report', help='generate report from gene and exon stats'
    )
    report_parser.add_argument(
        'annotated_bed',
        help='annotated bed file from annotated_bed'
    )
    report_parser.add_argument(
        'exon_stats',
        help='exon stats file (output/{sample_id}_exon_stats.tsv)'
    )
    report_parser.add_argument(
        'gene_stats',
        help='gene stats file (output/{sample_id}_gene_stats.tsv)'
    )
    report_parser.add_argument(
        '-s', '--snps', nargs='*',
        help='Optional; check coverage of VCF(s) of SNPs.'
    )
    report_parser.add_argument(
        '-t', '--threshold', nargs='?',
        default=20, type=int,
        help="threshold to define low coverage (int), if not\
            given 20 will be used as default. Must be one of\
            the thresholds in the input file.",
        required=False
    )
    report_parser.add_argument(
        '-n', '--sample_name', nargs='?',
        help="Name of sample to display in report, if not\
            specified this will be the prefix of the\
            gene_stats input file.",
        required=False
    )
    report_parser.add_argument(
        '-o', '--output', nargs='?',
        help='Output report name, if not specified the sample\
        name from the report will be used.',
        required=False
    )
    report_parser.add_argument(
        '-p', '--panel', nargs='?',
        help='(Optional) Panel bed file used from annotation, if passed\
        name of file will be displayed in report to show what\
        panel(s) / gene(s) were included.',
        required=False
    )
    report_parser.add_argument(
        '-l', '--limit', nargs='?',
        help="Number of genes at which to limit including full gene plots,\
        large numbers of genes takes a long time to generate the plots.",
        default=-1,
        required=False
    )
    report_parser.add_argument(
        '-m', '--summary',
        help="If passed, a short paragraph will be included in the\
        summary section. This includes details on the sequencing and the\
        genes/transcripts used in the panel.",
        default=False, action='store_true'
    )
    report_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )

    args = parser.parse_args()

    return args


def main():
    """
    Main function to do all things Athena
    """
    args = parse_args()
    load = loadData()
    # stats = singleCoverage()


    if not args.development:
        # read in required files
        panel_bed_df = load.read_panel_bed(args.panel_bed)
        transcript_info_df = load.read_transcript_info(args.transcript_file)
        pb_coverage_df = load.read_coverage_data(
            args.coverage_file, args.chunk_size
        )
        html_template = load.read_template()

        # normal analysis from raw data to report
        annotated_bed = annotate_bed(args, load)

        cov_stats, cov_summary = call_generate_stats(
            args, annotated_bed, args.thresholds, flagstat, build
        )

        call_generate_report(
            coverage_stats, coverage_summary, annotated_bed, args.threshold,
            args.snps, args.panel
        )
        # read in files
        cov_stats, cov_summary, raw_coverage, low_raw_cov, html_template,\
            flagstat, build, panel, vcfs, bootstrap, version = load_files(
                load,
                args.threshold,
                args.exon_stats,
                args.gene_stats,
                args.raw_coverage,
                args.snps,
                args.panel
            )

    else:
        # running sections individually
        if args.development == 'annotate':
            # generating annotated bed file
            annotated_bed = annotate_bed(args)

        if args.development == 'stats':
            # generating exon and gene stats
            # import data
            data, thresholds, flagstat, build = stats.import_data(args)
            call_generate_stats(args, data, thresholds, flagstat, build)

        if args.report == 'report':
            # generating report from stats
            call_generate_report(args, load)

    annotated_bed = annotate(args)


if __name__ == "__main__":
    main()
