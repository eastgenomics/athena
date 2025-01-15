"""Functions to handle bed file annotation"""

from pathlib import Path
import re
import subprocess
from timeit import default_timer as timer

from utils import log_handle
from .io import read_file, read_first_column
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
    RuntimeError
        Raised if any error emitted to stderr from bedtools intersect call
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

    # test if all chromosomes are in genome file to allow using -sorted
    defined_chromosomes = get_defined_chromosomes(genome)
    sample_chromosomes = get_coverage_chromosomes(coverage)
    undefined_chromosomes = set(sample_chromosomes) - defined_chromosomes
    sorted_arg = f"-sorted -g {genome}"

    if undefined_chromosomes:
        sorted_arg = ""
        log_handle.warning(
            "One or more chromosomes from coverage file not present in genome"
            " file: %s.\nWill use slower non sorted bedtools intersect",
            undefined_chromosomes,
        )

    try:
        proc = subprocess.run(
            f"bedtools intersect {sorted_arg} -wa -wb -a {regions} -b"
            f" {coverage} | cut -f7 --complement | gzip >"
            f" {output_file}",
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
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

    # subprocess doesn't seem to catch the non-zero exit code here,
    # manually check for error being dumped to stderr
    if "error" in proc.stderr.lower():
        raise RuntimeError(
            f"Error in calling bedtools intersect: {proc.stderr}"
        )

    log_handle.debug(
        "Annotated regions bed file completed in %s, written to %s",
        format_timer(start=start, end=timer()),
        output_file,
    )

    return output_file


def get_defined_chromosomes(genome_file: Path) -> list:
    """
    Reads the list of unique chromosomes defined in the genome file

    Parameters
    ----------
    genome_file : Path
        Path to genome file to read from

    Returns
    -------
    list
        List of unique chromosomes in the genome file
    """
    contents = read_file(file=genome_file)
    return set([x.split("\t")[0] for x in contents.splitlines()])


def get_coverage_chromosomes(coverage_file: Path) -> list:
    """
    Gets unique list of chromosomes from the given coverage file

    Parameters
    ----------
    coverage_file : Path
        Path to coverage file

    Returns
    -------
    list
        List of unique chromosomes in the coverage file
    """
    return (
        read_first_column(coverage_file, "chrom")
        .unique()
        .get_column("chrom")
        .to_list()
    )
