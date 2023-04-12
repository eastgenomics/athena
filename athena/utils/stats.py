"""
Functions relating to generating coverage stats from the per-base data
at given thresholds
"""

import argparse
from functools import partial
import os
import re
import sys
import math
import multiprocessing
import numpy as np
import pandas as pd

from natsort import natsorted
from pathlib import Path

from load import loadData


class stats():

    def single_sample_coverage_regions(sample_coverage) -> pd.DataFrame:
        """
        _summary_

        Parameters
        ----------
        sample_coverage : _type_
            _description_

        Returns
        -------
        pd.DataFrame
            _description_
        """