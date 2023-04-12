import numpy as np
import pandas as pd


def unbin(binned_data) -> pd.DataFrame:
    """
    Unbin binned coverage data to per base records in dataframe, dropping
    those outside of exon boundaries where the bin originally spanned
    the boundary.

    Parameters
    ----------
    binned_data : pd.DataFrame
        binned coverage data

    Returns
    -------
    pd.DataFrame
        unbinned coverage data
    """
    # generate list of single base values in range of coverage bins
    binned_data['position'] = list(map(list, list(map(
        range, binned_data['cov_start'].values,
        binned_data['cov_end'].values
    ))))

    # split out rows from range so there is one row per position
    unbinned_data = binned_data.explode('position')

    # drop rows where position falls outside of exon boundaries, this
    # occurs from the mosdepth coverage bins spanning the boundaries
    unbinned_data = unbinned_data.iloc[np.where((
        unbinned_data['exon_end'] > unbinned_data['position']
    ) & (
        unbinned_data['exon_start'] <= unbinned_data['position']
    ))]

    # drop columns that are no longer needed
    binned_data.drop(columns=['cov_start', 'cov_end'], inplace=True)

    return binned_data

