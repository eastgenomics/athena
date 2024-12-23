import polars as pl


DATAFRAME_TYPES = {
    "chrom": pl.Categorical,
    "region_start": pl.UInt32,
    "region_end": pl.UInt32,
    "gene": pl.Categorical,
    "transcript": pl.Categorical,
    "region": pl.Categorical,
    "region_length": pl.UInt16,
    "depth_bin_start": pl.UInt32,
    "depth_bin_end": pl.UInt32,
    "depth": pl.UInt16,
}
