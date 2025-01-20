"""Test dataframes for calculate.py function testing"""

import polars as pl
import pytest

from athena.utils import constants


class TotalPctCoverage:
    @pytest.fixture
    def input_coverage_bed_file(tmp_path):
        target_output = os.path.join(tmp_path, "coverage.bed")

        with open(target_output, "w+") as fh:
            fh.write(
                "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t604\t2556664",
                "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t606\t2556666",
                "chr1\t2556664\t2556733\tTNFRSF14\tNM_003820.4\t0\t610\t2556665",
            )

        return target_output

    required_cols = [
        "chrom",
        "region_start",
        "region_end",
        "gene",
        "transcript",
        "region",
        "position",
        "depth",
    ]

    # ensure we use the same predefined dtypes
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}

    coverage_bed_df = pl.DataFrame(
        [
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                604,
                2556664,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                606,
                2556666,
            ],
            [
                "chr1",
                2556664,
                2556733,
                "TNFRSF14",
                "NM_003820.4",
                0,
                610,
                2556665,
            ],
        ],
        schema=schema,
        orient="row",
    )
