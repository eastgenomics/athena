import argparse
import pathlib
import sys

from utils import log_handle


def parse_args() -> argparse.Namespace:
    """
    Parse provided command line arguments

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments
    """
    main_parser = argparse.ArgumentParser(add_help=False)

    subparsers = main_parser.add_subparsers(
        help="Select mode to run", dest="mode", required=True
    )

    report_parser = subparsers.add_parser(
        "report",
        help="Calculate all coverage for sample and generate report",
    )

    report_parser.add_argument(
        "-r",
        "--regions",
        required=True,
        help="Bed file of target regions to provide coverage data for",
    )

    report_parser.add_argument(
        "-c",
        "--coverage",
        required=True,
        help="Bed file of raw coverage data output from samtools / mosdepth",
    )
    report_parser.add_argument(
        "--normal_coverage",
        required=False,
        help=(
            "tsv of previously calculated normal coverage values for n"
            " samples. Requires hsmetrics file for sample providing to"
            " --hsmetrics."
        ),
    )
    report_parser.add_argument(
        "--hsmetrics",
        required="--normal_coverage" in sys.argv,
        help=(
            "hsmetrics file for current sample, required for --normal_coverage"
        ),
    )
    report_parser.add_argument("-a", "--annotated_bed", required=False)

    report_parser.add_argument(
        "-t",
        "--thresholds",
        type=int,
        nargs="+",
        default=[10, 20, 30, 50, 100],
        help="Thresholds at which to calculate percent coverage",
    )

    report_parser.add_argument(
        "-m",
        "--minimum",
        type=int,
        default=20,
        help=(
            "Minimum threshold value to use as cut off for defining as low"
            " coverage region. Must be one of --threshold values."
        ),
    )

    report_parser.add_argument(
        "--panel",
        type=str,
        required=False,
        help="Name of sequencing panel the report is for",
    )

    report_parser.add_argument(
        "--clinical_indication",
        type=str,
        required=False,
        help="Clinical indication the report is for",
    )

    report_parser.add_argument(
        "-b",
        "--build",
        type=int,
        help="Reference build of sample data",
    )

    report_parser.add_argument(
        "-o",
        "--output",
        required=False,
        help=(
            "Prefix for naming output files. Defaults to prefix of coverage"
            " bed file."
        ),
    )

    report_parser.add_argument(
        "--panel_filters",
        type=str,
        nargs="+",
        required=False,
        help=(
            "Preset filters of genes / transcripts to set for the full gene"
            " plots, these will be presented in a drop down menu for filtering"
            " the plots. These should be passed as key:value pairs of panel"
            " name to display in the drop down and a comma separated list of"
            " gene symbols to filter with. Example: 'Cancer:BRCA1,BRCA2'"
            " 'Cardiac:MYH7,TNNT2'"
        ),
    )

    report_parser.add_argument(
        "--summary",
        action="store_true",
        required=False,
        help=(
            "Display summary of genes / transcripts in report in summary"
            " section"
        ),
    )

    report_parser.add_argument(
        "--summary_file",
        action="store_true",
        required=False,
        help="Output text in summary section to a text file",
    )

    report_parser.add_argument(
        "--limit",
        default=-1,
        type=int,
        help=(
            "Number of genes at which to skip full gene plot generation. For"
            " large panels this significantly increases the report file size."
        ),
    )

    report_parser.add_argument(
        "--plot_chromosomes",
        action="store_true",
        default=False,
        required=False,
        help="Generates full chromosome plots of each chromosome",
    )

    report_parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        required=False,
        help="Force overwriting of existing files with same output name",
    )

    report_parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Increase logging verbosity to DEBUG level",
    )

    normal_coverage_parser = subparsers.add_parser(
        "calculate_normal",
        help=(
            "Calculate normalised mean per base coverage from multiple samples"
        ),
    )

    normal_coverage_parser.add_argument(
        "--regions",
        required=True,
        help="Bed file of target regions to calculate coverage data for",
    )
    normal_coverage_parser.add_argument(
        "--coverage",
        nargs="+",
        required=True,
        help=(
            "Bed files of coverage data output from samtools / mosdepth for"
            " all samples"
        ),
    )
    normal_coverage_parser.add_argument(
        "--hsmetrics",
        nargs="+",
        required=True,
        help="hsmetrics files for all samples",
    )
    normal_coverage_parser.add_argument(
        "-b",
        "--build",
        type=int,
        choices=[37, 38],
        required=True,
        help="Reference build of sample data",
    )
    normal_coverage_parser.add_argument(
        "--output", required=True, type=str, help="prefix for output file name"
    )
    normal_coverage_parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="Increase logging verbosity to DEBUG level",
    )

    args = main_parser.parse_args()

    if args.mode == "report":
        if not args.output:
            args.output = set_default_output_name(pathlib.Path(args.coverage))

        if not args.panel:
            args.panel = set_default_panel_name(pathlib.Path(args.regions))

        # TODO - abstract this into a set of checking functions
        if args.minimum not in args.thresholds:
            raise ValueError("--minimum must be one of --threshold values")

    return args


def set_default_output_name(coverage_file: pathlib.Path) -> str:
    """
    Sets the default output file name for the report from the given
    coverage data file prefix

    Parameters
    ----------
    coverage_file : pathlib.Path
        Path to coverage file

    Returns
    -------
    str
        Name for output report prefix
    """
    return coverage_file.name.replace(
        "".join(coverage_file.suffixes), ""
    ).replace("_markdup", "")


def set_default_panel_name(regions_bed_file: pathlib.Path) -> str:
    """
    Set the panel name to default from the regions bed file if not passed

    Parameters
    ----------
    regions_bed_file : pathlib.Path
        bed file for panel

    Returns
    -------
    str
        prefix of bed file
    """
    return regions_bed_file.name.replace(
        "".join(regions_bed_file.suffixes), ""
    )
