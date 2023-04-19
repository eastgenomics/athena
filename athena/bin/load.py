"""
Functions relating to loading of various required data files
"""
from pathlib import Path

import numpy as np
import pandas as pd


class loadData():
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
        Read in raw coverage data (annotated bed file) from single stats

        Parameters
        ----------
        raw_coverage : str
            filename of annotated bed file

        Returns
        -------
        pd.DataFrame
            dataframe of annotated bed file
        """
        annotated_bed = pd.read_csv(
            annotated_bed, sep="\t", header=0, dtype=self.dtypes
        )

        return annotated_bed


    def read_exon_stats(self, exon_stats) -> pd.DataFrame:
        """
        Read exon stats file output from stats.Sample() into dataframe

        Parameters
        ----------
        exon_stats : str
            filename of file with per exon coverage stats

        Returns
        -------
        pd.DataFrame
            dataframe of exon coverage stats
        pd.DataFrame | None
            dataframe of hsmetrics values, will return None if not present
        """
        sample = Path(exon_stats).name.replace('_exon_stats.tsv', '')

        with open(exon_stats) as exon_file:
            # check if hsmetrics have been written to first 2 lines of the file
            header = exon_file.readline()
            if header.startswith('BAIT_SET'):
                # hsmetrics present
                header = header.split('\t')
                metrics = [exon_file.readline().split('\t')]
                metrics = pd.DataFrame(metrics, columns=header)

                # set appropriate dtypes on columns we require for
                # calculating run stats normalisation from later\
                metrics = metrics.astype({
                    'ON_TARGET_BASES': int,
                    'PCT_USABLE_BASES_ON_TARGET': float
                })
            else:
                # no metrics present, test the file looks like normal
                # exon stats and read in
                assert header.startswith('chrom'), (
                    f'Exon stats file does not seem valid: {exon_file}'
                )
                metrics = None
                exon_file.seek(0)  # set file pointer back to start of file

            cov_stats = pd.read_csv(
                exon_file, sep="\t", dtype=self.dtypes
            )

        return cov_stats, metrics


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
