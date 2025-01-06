import argparse
import pathlib


def parse_args() -> argparse.Namespace:
    """
    Parse provided command line arguments

    Returns
    -------
    argparse.Namespace
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r",
        "--regions",
        required=False,
        help="Bed file of target regions to provide coverage data for",
    )

    parser.add_argument(
        "-c",
        "--coverage",
        required=False,
        help="Bed file of coverage data output from samtools / mosdepth",
    )
    parser.add_argument("-a", "--annotated_bed", required=False)

    parser.add_argument(
        "-t",
        "--thresholds",
        type=int,
        nargs="+",
        default=[10, 20, 30, 50, 100],
        help="Thresholds at which to calculate percent coverage",
    )

    parser.add_argument(
        "-m",
        "--minimum",
        type=int,
        default=20,
        help=(
            "Minimum threshold value to use as cut off for defining as low"
            " coverage region. Must be one of --threshold values."
        ),
    )
    parser.add_argument(
        "--panel",
        type=str,
        required=False,
        help="Name of sequencing panel the report is for",
    )
    parser.add_argument(
        "--clinical_indication",
        type=str,
        required=False,
        help="Clinical indication the report is for",
    )

    parser.add_argument(
        "-b",
        "--build",
        choices=[37, 38],
        type=int,
        help="Reference build of sample data",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=False,
        help=(
            "Prefix for naming output files. Defaults to prefix of coverage"
            " bed file."
        ),
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Increase logging verbosity to DEBUG level",
    )

    args = parser.parse_args()

    if not args.output:
        args.output = set_default_output_name(pathlib.Path(args.coverage))

    if not args.panel:
        args.panel = set_default_panel_name(pathlib.Path(args.regions))

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
