import numpy as np
import pandas as pd


def unbin(binned_data) -> pd.DataFrame:
    """
    Unbin binned coverage data to per base records in dataframe, dropping
    those outside of exon boundaries where the bin originally spanned
    the boundary.

    chr  exon_start exon_end  gene      tx         exon cov_start  cov_end  cov
    ---------------------------------------------------------------------------
    1    2488098    2488177  TNFRSF14  NM_003820.3  1    2488096   2488099  233
    1    2488098    2488177  TNFRSF14  NM_003820.3  1    2488099   2488100  236
    1    2488098    2488177  TNFRSF14  NM_003820.3  1    2488100   2488101  237
    1    2488098    2488177  TNFRSF14  NM_003820.3  1    2488101   2488104  235
    1    2488098    2488177  TNFRSF14  NM_003820.3  1    2488104   2488106  238

                                |
                                ▼

    chr   exon_start exon_end  gene      tx          exon   position   cov
    ----------------------------------------------------------------------
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488098    233
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488099    236
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488100    237
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488101    235
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488102    235
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488103    235
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488104    238
    1     2488098    2488177  TNFRSF14  NM_003820.3  1      2488105    238


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
        range, binned_data['cov_start'].values, binned_data['cov_end'].values
    ))))

    # split out rows from range so there is one row per position
    unbinned_data = binned_data.explode('position')
    unbinned_data.astype({'position': np.uint32})

    # drop rows where position falls outside of exon boundaries, this
    # occurs from the mosdepth coverage bins spanning the boundaries
    unbinned_data = unbinned_data.iloc[np.where((
        unbinned_data['exon_end'] > unbinned_data['position']
    ) & (
        unbinned_data['exon_start'] <= unbinned_data['position']
    ))].reset_index(drop=True)

    # drop columns that are no longer needed
    unbinned_data.drop(columns=['cov_start', 'cov_end'], axis=1, inplace=True)

    return unbinned_data
