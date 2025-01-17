"""Test dataframes for calculate.py function testing"""

import polars as pl


total_pct_coverage_df = pl.DataFrame(
    [["chr1", 1, 25], ["chr1", 2, 28], ["chr1", 2, 32], ["chr1", 2, 35]],
    schema={"chrom": pl.Categorical, "position": pl.Int32, "depth": pl.UInt32},
    orient="row",
)

# DataFrame with 1/10,000 base under 20 depth => 99.99% covered
total_pct_coverage_df_99_99 = pl.DataFrame(
    [["chr1", 1, 19], *[["chr1", x, 25] for x in range(2, 10_001)]],
    schema={"chrom": pl.Categorical, "position": pl.Int32, "depth": pl.UInt32},
    orient="row",
)
