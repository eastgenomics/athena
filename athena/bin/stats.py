"""
Functions relating to generating coverage stats from the per-base data
at given thresholds
"""
# from functools import partial
# import os
# import re
# import sys
# import math
# import multiprocessing
# import numpy as np
import pandas as pd
from typing import Union

# from natsort import natsorted
# from pathlib import Path

from .load import loadData
from .utils import unbin


class stats():
    """
    Generates coverage values for min, mean, max and given thresholds
    """
    def __init__(self) -> None:
        self.thresholds = []
        self.coverage_data = None
        self.gene_coverage = None
        self.exon_coverage = None


    def calculate(self) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        Calls all methods to generate coverage values for every region
        (i.e. exon) and summary values for every gene / transcript

        Returns
        -------
        pd.DataFrame
            dataframe of exon level coverage values
        pd.DataFrame
            dataframe of transcript level coverage values
        """
        self.generate_empty_dfs()
        self.coverage_data = unbin(self.coverage_data)


    def generate_empty_dfs(self):
        """
        Generate exon and gene level dataframes for adding values to,
        based off columns in given per base dataframe.
        """
        columns = [
            "chrom", "exon_start", "exon_end", "gene", "tx"
        ]

        self.gene_coverage = pd.DataFrame(columns=columns)
        self.exon_coverage = pd.DataFrame(columns=columns + ['exon'])


    
    def calculate_minimum(self, data, index):
        """
        Calculates the minimum coverage for each region from self.coverage_data

        Parameters
        ----------
        data : pd.DataFrame
            dataframe to which to add the calculated values
        index : list
            column names to use as index for 
        """
        return self.coverage_data.groupby(['chrom', 'exon_start', 'exon_end', 'gene', 'tx', 'exon'])['cov'].mean()
