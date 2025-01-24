"""Test data for test_calculate.TestPctThresholds"""

import polars as pl

from athena.utils import constants


def input_df() -> pl.DataFrame:
    """Minimal example of expected min, mean and max columns"""
    required_cols = ["chrom", "position", "gene", "region", "depth"]
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}

    return pl.DataFrame(
        [
            # gene with 2 regions and some coverage
            ["chr1", 1, "BRCA1", "1", 20],
            ["chr1", 2, "BRCA1", "1", 22],
            ["chr1", 3, "BRCA1", "1", 30],
            ["chr1", 4, "BRCA1", "1", 31],
            ["chr1", 1, "BRCA1", "2", 10],
            ["chr1", 2, "BRCA1", "2", 15],
            ["chr1", 3, "BRCA1", "2", 14],
            # gene with a single region with no coverage
            ["chr1", 1, "PALB2", "1", 30],
            ["chr1", 2, "PALB2", "1", 20],
            ["chr1", 1, "PALB2", "2", 0],
            ["chr1", 2, "PALB2", "2", 0],
            # whole gene with no coverage
            ["chr2", 1, "EGFR", "1", 0],
            ["chr2", 2, "EGFR", "1", 0],
            ["chr2", 3, "EGFR", "1", 0],
        ],
        schema=schema,
        orient="row",
    )


def expected_pct_coverage_per_gene() -> pl.DataFrame:
    """
    Below DataFrame construction gives the following with the expected
    per gene coverage values for 10x, 20x and 30x based off of the above
    input_dataframe():

    ┌───────┬──────────┬───────┬────────┬───────┬───────┬──────────┬──────────┐
    │ chrom ┆ position ┆ gene  ┆ region ┆ depth ┆ 10x   ┆ 20x      ┆ 30x      │
    │ ---   ┆ ---      ┆ ---   ┆ ---    ┆ ---   ┆ ---   ┆ ---      ┆ ---      │
    │ cat   ┆ u32      ┆ cat   ┆ cat    ┆ u32   ┆ f32   ┆ f32      ┆ f32      │
    ╞═══════╪══════════╪═══════╪════════╪═══════╪═══════╪══════════╪══════════╡
    │ chr1  ┆ 1        ┆ BRCA1 ┆ 1      ┆ 20    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 2        ┆ BRCA1 ┆ 1      ┆ 22    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 3        ┆ BRCA1 ┆ 1      ┆ 30    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 4        ┆ BRCA1 ┆ 1      ┆ 31    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 1        ┆ BRCA1 ┆ 2      ┆ 10    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 2        ┆ BRCA1 ┆ 2      ┆ 15    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 3        ┆ BRCA1 ┆ 2      ┆ 14    ┆ 100.0 ┆ 57.14286 ┆ 28.57143 │
    │ chr1  ┆ 1        ┆ PALB2 ┆ 1      ┆ 30    ┆ 50.0  ┆ 50.0     ┆ 25.0     │
    │ chr1  ┆ 2        ┆ PALB2 ┆ 1      ┆ 20    ┆ 50.0  ┆ 50.0     ┆ 25.0     │
    │ chr1  ┆ 1        ┆ PALB2 ┆ 2      ┆ 0     ┆ 50.0  ┆ 50.0     ┆ 25.0     │
    │ chr1  ┆ 2        ┆ PALB2 ┆ 2      ┆ 0     ┆ 50.0  ┆ 50.0     ┆ 25.0     │
    │ chr2  ┆ 1        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0      ┆ 0.0      │
    │ chr2  ┆ 2        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0      ┆ 0.0      │
    │ chr2  ┆ 3        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0      ┆ 0.0      │
    └───────┴──────────┴───────┴────────┴───────┴───────┴──────────┴──────────┘
    """
    required_cols = [
        "chrom",
        "position",
        "gene",
        "region",
        "depth",
    ]
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}
    schema = {
        **schema,
        **{"10x": pl.Float32, "20x": pl.Float32, "30x": pl.Float32},
    }

    return pl.DataFrame(
        [
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "BRCA1",
                "region": "1",
                "depth": 20,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "BRCA1",
                "region": "1",
                "depth": 22,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 3,
                "gene": "BRCA1",
                "region": "1",
                "depth": 30,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 4,
                "gene": "BRCA1",
                "region": "1",
                "depth": 31,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "BRCA1",
                "region": "2",
                "depth": 10,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "BRCA1",
                "region": "2",
                "depth": 15,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 3,
                "gene": "BRCA1",
                "region": "2",
                "depth": 14,
                "10x": 100.0,
                "20x": 57.142860412597656,
                "30x": 28.571430206298828,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "PALB2",
                "region": "1",
                "depth": 30,
                "10x": 50.0,
                "20x": 50.0,
                "30x": 25.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "PALB2",
                "region": "1",
                "depth": 20,
                "10x": 50.0,
                "20x": 50.0,
                "30x": 25.0,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "PALB2",
                "region": "2",
                "depth": 0,
                "10x": 50.0,
                "20x": 50.0,
                "30x": 25.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "PALB2",
                "region": "2",
                "depth": 0,
                "10x": 50.0,
                "20x": 50.0,
                "30x": 25.0,
            },
            {
                "chrom": "chr2",
                "position": 1,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr2",
                "position": 2,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr2",
                "position": 3,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
        ],
        schema=schema,
    )


def expected_pct_coverage_per_region() -> pl.DataFrame:
    """
    Below DataFrame construction gives the following with the expected
    per region coverage values for 10x, 20x and 30x based off of the above
    input_dataframe():

    ┌───────┬──────────┬───────┬────────┬───────┬───────┬───────┬──────┐
    │ chrom ┆ position ┆ gene  ┆ region ┆ depth ┆ 10x   ┆ 20x   ┆ 30x  │
    │ ---   ┆ ---      ┆ ---   ┆ ---    ┆ ---   ┆ ---   ┆ ---   ┆ ---  │
    │ cat   ┆ u32      ┆ cat   ┆ cat    ┆ u32   ┆ f32   ┆ f32   ┆ f32  │
    ╞═══════╪══════════╪═══════╪════════╪═══════╪═══════╪═══════╪══════╡
    │ chr1  ┆ 1        ┆ BRCA1 ┆ 1      ┆ 20    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 2        ┆ BRCA1 ┆ 1      ┆ 22    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 3        ┆ BRCA1 ┆ 1      ┆ 30    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 4        ┆ BRCA1 ┆ 1      ┆ 31    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 1        ┆ BRCA1 ┆ 2      ┆ 10    ┆ 100.0 ┆ 0.0   ┆ 0.0  │
    │ chr1  ┆ 2        ┆ BRCA1 ┆ 2      ┆ 15    ┆ 100.0 ┆ 0.0   ┆ 0.0  │
    │ chr1  ┆ 3        ┆ BRCA1 ┆ 2      ┆ 14    ┆ 100.0 ┆ 0.0   ┆ 0.0  │
    │ chr1  ┆ 1        ┆ PALB2 ┆ 1      ┆ 30    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 2        ┆ PALB2 ┆ 1      ┆ 20    ┆ 100.0 ┆ 100.0 ┆ 50.0 │
    │ chr1  ┆ 1        ┆ PALB2 ┆ 2      ┆ 0     ┆ 0.0   ┆ 0.0   ┆ 0.0  │
    │ chr1  ┆ 2        ┆ PALB2 ┆ 2      ┆ 0     ┆ 0.0   ┆ 0.0   ┆ 0.0  │
    │ chr2  ┆ 1        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0   ┆ 0.0  │
    │ chr2  ┆ 2        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0   ┆ 0.0  │
    │ chr2  ┆ 3        ┆ EGFR  ┆ 1      ┆ 0     ┆ 0.0   ┆ 0.0   ┆ 0.0  │
    └───────┴──────────┴───────┴────────┴───────┴───────┴───────┴──────┘
    """
    required_cols = [
        "chrom",
        "position",
        "gene",
        "region",
        "depth",
    ]
    schema = {k: constants.DATAFRAME_TYPES.get(k) for k in required_cols}
    schema = {
        **schema,
        **{"10x": pl.Float32, "20x": pl.Float32, "30x": pl.Float32},
    }

    return pl.DataFrame(
        [
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "BRCA1",
                "region": "1",
                "depth": 20,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "BRCA1",
                "region": "1",
                "depth": 22,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 3,
                "gene": "BRCA1",
                "region": "1",
                "depth": 30,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 4,
                "gene": "BRCA1",
                "region": "1",
                "depth": 31,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "BRCA1",
                "region": "2",
                "depth": 10,
                "10x": 100.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "BRCA1",
                "region": "2",
                "depth": 15,
                "10x": 100.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr1",
                "position": 3,
                "gene": "BRCA1",
                "region": "2",
                "depth": 14,
                "10x": 100.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "PALB2",
                "region": "1",
                "depth": 30,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "PALB2",
                "region": "1",
                "depth": 20,
                "10x": 100.0,
                "20x": 100.0,
                "30x": 50.0,
            },
            {
                "chrom": "chr1",
                "position": 1,
                "gene": "PALB2",
                "region": "2",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr1",
                "position": 2,
                "gene": "PALB2",
                "region": "2",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr2",
                "position": 1,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr2",
                "position": 2,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
            {
                "chrom": "chr2",
                "position": 3,
                "gene": "EGFR",
                "region": "1",
                "depth": 0,
                "10x": 0.0,
                "20x": 0.0,
                "30x": 0.0,
            },
        ],
        schema=schema,
    )
