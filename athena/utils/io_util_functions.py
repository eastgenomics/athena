from pathlib import Path
from timeit import default_timer as timer
from typing import Dict, List

from .io import read_hsmetrics
from .util_functions import format_timer


from .log import get_logger

log_handle = get_logger("athena")


def filter_samples_by_minimum_reads(
    sample_files: dict, min_reads: int
) -> Dict[List[Path], List[Path]]:
    """
    Filter the provided sample files by having `minimum_reads` as read
    from the hsmetrics file. This is to ensure no low sequenced samples
    are included in our calculations.

    Parameters
    ----------
    sample_files : dict
        Mapping of sample ID -> coverage bed file and hsmetrics file paths
    min_reads : int
        Minimum n TOTAL_READS to filter out samples below

    Returns
    -------
    Dict[List[Path], List[Path]]
        Same mapping as `sample_files`, with samples < min_reads removed
    """
    log_handle.debug(
        "Filtering %s samples for having > %s total reads",
        len(sample_files),
        f"{min_reads:,}",
    )
    start = timer()

    filtered_sample_files = {}

    for sample, files in sample_files.items():
        metrics = read_hsmetrics(files[1])

        if int(metrics["TOTAL_READS"].item()) < min_reads:
            log_handle.warning(
                "%s has less than %s TOTAL_READS (%s), will be excluded",
                Path(sample).name,
                min_reads,
                metrics["TOTAL_READS"].item(),
            )
        else:
            filtered_sample_files[sample] = files

    log_handle.debug(
        "Completed filtering samples in %s.\n%s samples left after filtering"
        " for minimum reads",
        format_timer(start=start, end=timer()),
        len(filtered_sample_files),
    )

    return filtered_sample_files
