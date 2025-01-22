""" "Data for report.get_total_sub_threshold_genes_and_regions tests"""

import polars as pl

dtypes = {
    "gene": pl.Categorical,
    "region": pl.Categorical,
    "10x": pl.Float32,
    "30x": pl.Float32,
}


def empty_df() -> pl.DataFrame:
    return pl.DataFrame(
        {"gene": [], "region": [], "10x": [], "30x": []},
        schema=dtypes,
    )


def regions_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "gene": ["gene_1", "gene_1", "gene_2"],
            "region": ["exon_1", "exon_2", "exon_1"],
            "10x": [100, 100, 100],
            "30x": [100, 94.9, 100],
        },
        schema=dtypes,
    )
