"""
Functions relating to loading of various required data files
"""
import numpy as np
import pandas as pd


class LoadData():
    def __init__(self):
        self.dtypes = {
            "chrom": str,
            "exon_start": np.uint32,
            "exon_end": np.uint32,
            "start": np.uint32,
            "end": np.uint32,
            "gene": str,
            "tx": str,
            "transcript": str,
            "exon": str,
            "exon_len": np.uint16,
            "min": np.uint16,
            "mean": float,
            "max": np.uint16,
            "cov_start": np.uint32,
            "cov_end": np.uint32,
            "cov": np.uint16,
            r'[0-9]*x': float
        }


    def filter_dtypes(self, dataframe) -> dict:
        """
        Filter list of dtypes to columns present in df
        Parameters
        ----------
        dataframe : pd.DataFrame
            dataframe to get column names from to select required dtypes

        Returns
        -------
        dict
            dict of column names and dtypes
        """
        return {k: v for k, v in self.dtypes.items() if k in dataframe.columns}


    def read_annotated_bed(self, annotated_bed):
        """
        Read in annotated bed file with per base coverage information for
        the target regions

        Parameters
        ----------
        annotated_bed : str
            filename of annotated bed file

        Returns
        -------
        pd.DataFrame
            dataframe of annotated bed file
        """
        column = [
            "chrom", "exon_start", "exon_end", "gene",
            "transcript", "exon", "cov_start", "cov_end", "cov"
        ]

        annotated_bed = pd.read_csv(
            annotated_bed, sep="\t", names=column, dtype=self.dtypes
        )

        return annotated_bed


    @staticmethod
    def read_hsmetrics(hs_file):
        """
        Read in hsmetrics file to get values for normalisation when
        generating run level stats

        Parameters
        ----------
        hs_file : str
            path to hsmetrics file to read in

        Returns
        -------
        pd.DataFrame
            dataframe of hsmetrics values
        """
        with open(hs_file) as file:
            contents = file.read().splitlines()

        metrics = []

        for idx, line in enumerate(contents):
            if line.startswith('## METRICS CLASS'):
                metrics.extend(contents[idx + 1: idx + 3])

        metrics = [x.split('\t') for x in metrics]

        assert metrics, "METRICS CLASS could not be parsed from hsmetrics file"

        return pd.DataFrame(metrics)
