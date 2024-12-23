import polars as pl


def unbin(coverage_data: pl.DataFrame) -> pl.DataFrame:
    """
    Unbin binned coverage data to per base records in DataFrame, dropping
    those outside of region boundaries where the bin originally spanned
    the boundary. This will return a DataFrame with one row per position
    in the given regions.

    ┌───────┬──────────────┬────────────┬───────┬───┬────────┬─────────────────┬───────────────┬───────┐
    │ chrom ┆ region_start ┆ region_end ┆ gene  ┆ … ┆ region ┆ depth_bin_start ┆ depth_bin_end ┆ depth │
    │ ---   ┆ ---          ┆ ---        ┆ ---   ┆   ┆ ---    ┆ ---             ┆ ---           ┆ ---   │
    │ cat   ┆ u32          ┆ u32        ┆ cat   ┆   ┆ cat    ┆ u32             ┆ u32           ┆ u32   │
    ╞═══════╪══════════════╪════════════╪═══════╪═══╪════════╪═════════════════╪═══════════════╪═══════╡
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ … ┆ foo    ┆ 108098346       ┆ 108098347     ┆ 783   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ … ┆ foo    ┆ 108098347       ┆ 108098348     ┆ 771   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ … ┆ foo    ┆ 108098348       ┆ 108098349     ┆ 772   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ … ┆ foo    ┆ 108098349       ┆ 108098350     ┆ 777   │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ … ┆ foo    ┆ 108098350       ┆ 108098351     ┆ 792   │
    │ …     ┆ …            ┆ …          ┆ …     ┆ … ┆ …      ┆ …               ┆ …             ┆ …     │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ … ┆ 2      ┆ 29130709        ┆ 29130710      ┆ 1001  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ … ┆ 2      ┆ 29130710        ┆ 29130711      ┆ 999   │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ … ┆ 2      ┆ 29130711        ┆ 29130712      ┆ 1003  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ … ┆ 2      ┆ 29130712        ┆ 29130713      ┆ 990   │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ … ┆ 2      ┆ 29130713        ┆ 29130714      ┆ 993   │
    └───────┴──────────────┴────────────┴───────┴───┴────────┴─────────────────┴───────────────┴───────┘

                                                    ↓

    ┌───────┬──────────────┬────────────┬───────┬─────────────┬────────┬───────┬───────────┐
    │ chrom ┆ region_start ┆ region_end ┆ gene  ┆ transcript  ┆ region ┆ depth ┆ position  │
    │ ---   ┆ ---          ┆ ---        ┆ ---   ┆ ---         ┆ ---    ┆ ---   ┆ ---       │
    │ cat   ┆ u32          ┆ u32        ┆ cat   ┆ cat         ┆ cat    ┆ u32   ┆ u32       │
    ╞═══════╪══════════════╪════════════╪═══════╪═════════════╪════════╪═══════╪═══════════╡
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ NM_000051.4 ┆ foo    ┆ 783   ┆ 108098346 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ NM_000051.4 ┆ foo    ┆ 771   ┆ 108098347 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ NM_000051.4 ┆ foo    ┆ 772   ┆ 108098348 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ NM_000051.4 ┆ foo    ┆ 777   ┆ 108098349 │
    │ 11    ┆ 108098346    ┆ 108098428  ┆ ATM   ┆ NM_000051.4 ┆ foo    ┆ 792   ┆ 108098350 │
    │ …     ┆ …            ┆ …          ┆ …     ┆ …           ┆ …      ┆ …     ┆ …         │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ NM_007194.4 ┆ 2      ┆ 1001  ┆ 29130709  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ NM_007194.4 ┆ 2      ┆ 999   ┆ 29130710  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ NM_007194.4 ┆ 2      ┆ 1003  ┆ 29130711  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ NM_007194.4 ┆ 2      ┆ 990   ┆ 29130712  │
    │ 22    ┆ 29130385     ┆ 29130714   ┆ CHEK2 ┆ NM_007194.4 ┆ 2      ┆ 993   ┆ 29130713  │
    └───────┴──────────────┴────────────┴───────┴─────────────┴────────┴───────┴───────────┘


    Parameters
    ----------
    coverage_data : pl.DataFrame
        DataFrame on which to calculate values

    Returns
    -------
    pd.DataFrame
        unbinned coverage data
    """

    coverage_data = coverage_data.with_columns(
        position=pl.int_ranges(
            start="depth_bin_start", end="depth_bin_end", dtype=pl.UInt32
        )
    )
    coverage_data = coverage_data.drop(["depth_bin_start", "depth_bin_end"])
    coverage_data = coverage_data.explode("position")
    coverage_data = coverage_data.filter(
        (pl.col("region_start") <= pl.col("position"))
        & (pl.col("position") < pl.col("region_end"))
    )

    return coverage_data
