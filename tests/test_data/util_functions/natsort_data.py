"""Test data for util_functions.natsort"""

import polars as pl


class MinimalExample:
    @staticmethod
    def unsorted() -> pl.DataFrame:
        return pl.DataFrame(
            {
                "chrom": ["chr1", "chr2", "chr1", "chr2"],
                "pos": [101, 300, 100, 301],
                "depth": [33, 50, 35, 52],
            },
            schema={
                "chrom": pl.Categorical,
                "pos": pl.UInt32,
                "depth": pl.UInt32,
            },
        )

    @staticmethod
    def sorted_by_int_column() -> pl.DataFrame:
        return pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2"],
                "pos": [100, 101, 300, 301],
                "depth": [35, 33, 50, 52],
            },
            schema={
                "chrom": pl.Categorical,
                "pos": pl.UInt32,
                "depth": pl.UInt32,
            },
        )

    @staticmethod
    def sorted_by_categorical_and_int_column() -> pl.DataFrame():
        return pl.DataFrame(
            {
                "chrom": ["chr1", "chr1", "chr2", "chr2"],
                "pos": [101, 100, 300, 301],
                "depth": [33, 35, 50, 52],
            },
            schema={
                "chrom": pl.Categorical,
                "pos": pl.UInt32,
                "depth": pl.UInt32,
            },
        )


class MessyStringColumnData:
    @staticmethod
    def unsorted() -> pl.DataFrame:
        return pl.DataFrame(
            {"col_1": ["foo", "bar", "100.10", "100", "a1"]},
            schema={"col_1": pl.String},
        )

    @staticmethod
    def sorted() -> pl.DataFrame:
        return pl.DataFrame(
            {"col_1": ["100", "100.10", "a1", "bar", "foo"]},
            schema={"col_1": pl.String},
        )


class FloatAndIntColumns:
    @staticmethod
    def unsorted() -> pl.DataFrame:
        return pl.DataFrame(
            {
                "float_col": [10.10, 10.10, 23.23, 0.555],
                "int_col": [5, 3, 1, 25],
            },
            schema={"float_col": pl.Float32, "int_col": pl.UInt32},
        )

    @staticmethod
    def sorted_by_float_and_int() -> pl.DataFrame:
        return pl.DataFrame(
            {
                "float_col": [0.555, 10.10, 10.10, 23.23],
                "int_col": [25, 3, 5, 1],
            },
            schema={"float_col": pl.Float32, "int_col": pl.UInt32},
        )
