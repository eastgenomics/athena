"""Functions to handle bed file annotation"""

from pathlib import Path
import re
import subprocess
from timeit import default_timer as timer

from utils import log_handle
from .util_functions import format_timer


def call_bedtools_intersect(
    coverage: str, regions: str, build: int, overwrite: bool
) -> str:
    """
    Calls bedtools intersect via subshell to annotate the `regions` bed
    file with per base coverage data from the `coverage` bed file.

    Parameters
    ----------
    coverage : str
        bed file of per base coverage data
    regions : str
        bed file of regions to annotate
    build : int
        reference build the data is for, used for specifying genome
        file to bedtools intersect
    overwrite : bool
        forces overwriting of output file if exists

    Returns
    -------
    str
        file name of annotated regions bed file

    Outputs
    -------
    file
        annotated regions bed file

    Raises
    ------
    ValueError
        Raised when build is not one of 37 or 38
    FileNotFoundError
        Raised if either regions or coverage bed file does not exist
    FileExistsError
        Raised when a file with the proposed filename already exists
    subprocess.CalledProcessError
        Raised if a non-zero exit code returned from subprocess.run
    """
    log_handle.debug(
        "Annotating regions bed file via bedtools intersect with coverage data"
    )
    start = timer()

    if build == 37:
        genome = Path(__file__).parent.parent.joinpath(
            "data/genomes/human.hg19.genome"
        )
    elif build == 38:
        genome = Path(__file__).parent.parent.joinpath(
            "data/genomes/human.hg38.genome"
        )
    else:
        raise ValueError("build must be one of 37 or 38")

    if not Path(regions).exists():
        raise FileNotFoundError(f"regions bed file does not exist: {regions}")

    if not Path(coverage).exists():
        raise FileNotFoundError(
            f"coverage bed file does not exist: {coverage}"
        )

    output_file = re.sub(rf"{''.join(Path(coverage).suffixes)}$", "", coverage)
    output_file += ".coverage.bed.gz"

    if Path(output_file).exists() and not overwrite:
        raise FileExistsError(
            f"Output file {output_file} already exists, stopping now to not"
            " overwrite."
        )

    try:
        subprocess.run(
            f"bedtools intersect -sorted -nonamecheck -g {genome} -wa -wb -a"
            f" {regions} -b {coverage} | cut -f7 --complement | gzip >"
            f" {output_file}",
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError as err:
        raise subprocess.CalledProcessError(
            returncode=err.returncode,
            cmd=err.cmd,
            output=err.output,
            stderr=(
                f"Error in calling bedtools intersect: {err.stderr.decode()}"
            ),
        ) from err

    log_handle.debug(
        "Annotated regions bed file completed in %s, written to %s",
        format_timer(start=start, end=timer()),
        output_file,
    )

    return output_file
