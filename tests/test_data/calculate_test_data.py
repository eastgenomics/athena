"""Test dataframes for calculate.py function testing"""

import polars as pl


class TotalPctCoverage:
    dataframe = pl.DataFrame(
        [["chr1", 1, 25], ["chr1", 2, 28], ["chr1", 2, 32], ["chr1", 2, 35]],
        schema={
            "chrom": pl.Categorical,
            "position": pl.Int32,
            "depth": pl.UInt32,
        },
        orient="row",
    )

    # DataFrame with 1/10,000 base under 20 depth => 99.99% covered
    dataframe_99_99_pct = pl.DataFrame(
        [["chr1", 1, 19], *[["chr1", x, 25] for x in range(2, 10_001)]],
        schema={
            "chrom": pl.Categorical,
            "position": pl.Int32,
            "depth": pl.UInt32,
        },
        orient="row",
    )


class MinMeanMax:
    """Test data for calculate.min_mean_max"""

    input_calculated_columns_df = pl.DataFrame(
        [
            ["chr1", 1, "BRCA1", "1", 20],
            ["chr1", 2, "BRCA1", "1", 30],
            ["chr1", 3, "BRCA1", "1", 40],
            ["chr1", 1, "BRCA1", "2", 10],
            ["chr1", 2, "BRCA1", "2", 15],
            ["chr1", 3, "BRCA1", "2", 20],
            ["chr1", 1, "PALB2", "1", 30],
            ["chr1", 2, "PALB2", "1", 30],
            ["chr1", 1, "PALB2", "2", 0],
        ],
        schema={
            "chrom": pl.Categorical,
            "position": pl.UInt32,
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "depth": pl.UInt32,
        },
        orient="row",
    )

    expected_grouped_by_gene_df = pl.DataFrame(
        [["BRCA1", 10, 22.5, 40], ["PALB2", 0, 20.0, 30]],
        schema={
            "gene": pl.Categorical,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )

    expected_grouped_by_gene_and_region_df = pl.DataFrame(
        [
            ["BRCA1", "1", 20, 30, 40],
            ["BRCA1", "2", 10, 15.0, 20],
            ["PALB2", "1", 30, 30.0, 30],
            ["PALB2", "2", 0, 0.0, 0],
        ],
        schema={
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )

    expected_joined_to_input_df = pl.DataFrame(
        [
            ["chr1", 1, "BRCA1", "1", 20, 20, 30, 40],
            ["chr1", 2, "BRCA1", "1", 30, 20, 30, 40],
            ["chr1", 3, "BRCA1", "1", 40, 20, 30, 40],
            ["chr1", 1, "BRCA1", "2", 10, 10, 15.0, 20],
            ["chr1", 2, "BRCA1", "2", 15, 10, 15.0, 20],
            ["chr1", 3, "BRCA1", "2", 20, 10, 15.0, 20],
            ["chr1", 1, "PALB2", "1", 30, 30, 30.0, 30],
            ["chr1", 2, "PALB2", "1", 30, 30, 30.0, 30],
            ["chr1", 1, "PALB2", "2", 0, 0, 0.0, 0],
        ],
        schema={
            "chrom": pl.Categorical,
            "position": pl.UInt32,
            "gene": pl.Categorical,
            "region": pl.Categorical,
            "depth": pl.UInt32,
            "min": pl.UInt32,
            "mean": pl.Float32,
            "max": pl.UInt32,
        },
        orient="row",
    )
