""" "Data for report.get_total_fully_covered_genes tests"""

import polars as pl

dtypes = {
    "chrom": pl.Categorical,
    "pos": pl.UInt32,
    "gene": pl.Categorical,
    "10x": pl.Float32,
    "30x": pl.Float32,
}


def empty_df() -> pl.DataFrame:
    return pl.DataFrame(
        {"chrom": [], "pos": [], "gene": [], "10x": [], "30x": []},
        schema=dtypes,
    )


def regions_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "pos": [1, 2, 1],
            "gene": ["gene_1", "gene_1", "gene_2"],
            "10x": [100, 100, 100],
            "30x": [100, 94.9, 100],
        },
        schema=dtypes,
    )
