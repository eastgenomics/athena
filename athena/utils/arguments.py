import argparse


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
        required=True,
        help="bed file of target regions to provide coverage data for",
    )

    parser.add_argument(
        "-c",
        "--coverage",
        required=True,
        help="bed file of coverage data output from mosdepth",
    )

    parser.add_argument(
        "-t",
        "--thresholds",
        type=int,
        nargs="+",
        default=[10, 20, 30, 50, 100],
        help="thresholds at which to calculate percent coverage",
    )

    parser.add_argument(
        "-m",
        "--minimum",
        type=int,
        default=20,
        help=(
            "minimum threshold value to use as cut off for defining as low"
            " coverage region. Must be one of --threshold values."
        ),
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
            "prefix for naming output files. Defaults to prefix of coverage"
            " bed file."
        ),
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="increase logging verbosity to DEBUG level",
    )

    return parser.parse_args()
