"""Test data for report.generate_summary_text"""

import polars as pl


def dataframe() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "gene": ["BRCA1", "BRCA2", "PALB2", "EGFR"],
            "transcript": [
                "transcript_1",
                "transcript_2",
                "transcript_3",
                "transcript_4",
            ],
            "10x": [100.0, 100.0, 100.0, 100.0],
            "30x": [100.0, 100.0, 99.99, 89.12],
        }
    )
