"""
Functions relating to generating coverage stats from the per-base data
at given thresholds
"""
from functools import partial
import multiprocessing
import sys

from natsort import natsorted
import numpy as np
import pandas as pd

# from load import loadData
from utils import unbin


class stats():
    """
    Generates coverage values for min, mean, max and given thresholds
    """
    def calculate_exon_stats(self, coverage_data, thresholds) -> pd.DataFrame:
        """
        Calls all methods to generate coverage values for every region
        (i.e. exon) and summary values for every gene / transcript

        Parameters
        ----------
        coverage_data : pd.DataFrame
            Dataframe of per base coverage data
        thresholds : list
            list of thresholds at which to calculate coverage

        Returns
        -------
        pd.DataFrame
            dataframe of exon level coverage values
        pd.DataFrame
            dataframe of transcript level coverage values
        """
        gene_coverage, exon_coverage = self.generate_empty_dfs(coverage_data)
        coverage_data = unbin(coverage_data)

        # calculate required exon level stats
        exon_coverage = self.calculate_minimum(
            data=coverage_data,
            output=exon_coverage,
            index=['tx', 'exon'],
            column='cov'
        )

        exon_coverage = self.calculate_mean(
            data=coverage_data,
            output=exon_coverage,
            index=['tx', 'exon'],
            column='cov'
        )

        exon_coverage = self.calculate_maximum(
            data=coverage_data,
            output=exon_coverage,
            index=['tx', 'exon'],
            column='cov'
        )

        exon_coverage = self.calculate_exon_thresholds(
            data=coverage_data,
            output=exon_coverage,
            index=['tx', 'exon'],
            thresholds=thresholds
        )

        return exon_coverage

           # get list of transcripts in data


    def calculate_exon_stats_parallel(self, coverage_data, thresholds) -> pd.DataFrame:
        """
        Wrapper method to call calculate_exon_stats in parallel using
        multiprocessing.Pool() and maximum number of CPU cores available

        Parameters
        ----------
        coverage_data : pd.DataFrame
            Dataframe of per base coverage data
        thresholds : list
            list of thresholds at which to calculate coverage

        Returns
        -------
        pd.DataFrame
            dataframe of exon level coverage values
        """
        # split per base dataframe into an array of dataframes by transcript
        # for passing to calculate_exon_stats() in parrallel
        transcripts = sorted(coverage_data.tx.unique().tolist())
        tx_array = np.array_split(np.array(transcripts), multiprocessing.cpu_count())
        split_dfs = np.asanyarray(
            [coverage_data[coverage_data["tx"].isin(x)]
            for x in tx_array], dtype=object
        )

        with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            try:
                # map threshold arg to function since same is passed to all, then
                # use imap to pass each df to worker to generate stats
                stats_w_arg = partial(
                    self.calculate_exon_stats, thresholds=thresholds
                )
                pool_output = pool.imap_unordered(stats_w_arg, split_dfs)
            except Exception as err:
                print(f"Error in calculating exon stats: {err}")
                pool.close()
                pool.terminate()
            else:
                exon_stats = pd.concat(pool_output, ignore_index=True)

        # correctly sort by gene and exon
        exon_stats.exon = exon_stats.exon.astype('category')
        exon_stats.exon.cat.reorder_categories(
            natsorted(set(exon_stats.exon)), inplace=True, ordered=True)
        exon_stats.sort_values(by=['gene', 'tx', 'exon'], inplace=True)

        return exon_stats


    def generate_empty_dfs(self, per_base_data):
        """
        Generate exon and gene level dataframes for adding values to,
        based off columns in given per base dataframe, keeping just
        single rows for each unique region from the per base rows.

        Parameters
        ----------
        per_base_data : pd.DataFrame
            dataframe of per base coverage data to extract unique exon
            and gene rows from
        
        Returns
        -------
        gene_coverage : pd.DataFrame
            dataframe of unique gene rows
        exon_coverage : pd.DataFrame
            dataframe of unique exon rows
        """
        gene_coverage = per_base_data[[
            'chrom', 'gene', 'tx'
        ]].drop_duplicates().reset_index(drop=True)
        exon_coverage = per_base_data[[
            'chrom', 'exon_start', 'exon_end', 'gene', 'tx', 'exon'
        ]].drop_duplicates().reset_index(drop=True)

        return gene_coverage, exon_coverage


    def calculate_minimum(self, data, output, index, column):
        """
        Calculates the minimum coverage for each region

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of from which calculate values
        output : pd.DataFrame
            dataframe to which to add the calculated values to
        index : list
            column names to use as index for calculating against
        column : str
            column from which to calculate

        Returns
        -------
        pd.DataFrame
            output dataframe with min column appended
        """
        minimum = data.groupby(index)[column].min().reset_index()
        minimum.rename(columns={column: 'min'}, inplace=True)

        output = pd.merge(output, minimum, on=index, validate='1:1')

        return output


    def calculate_maximum(self, data, output, index, column) -> pd.DataFrame:
        """
        Calculates the maximum coverage for each region

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of from which calculate values
        output : pd.DataFrame
            dataframe to which to add the calculated values to
        index : list
            column names to use as index for calculating against
        column : str
            column from which to calculate
        
        Returns
        -------
        pd.DataFrame
            output dataframe with max column appended
        """
        maximum = data.groupby(index)[column].max().reset_index()
        maximum.rename(columns={column: 'max'}, inplace=True)

        return pd.merge(output, maximum, on=index, validate='1:1')


    def calculate_mean(self, data, output, index, column) -> pd.DataFrame:
        """
        Calculates the mean coverage for each region

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of from which calculate values
        output : pd.DataFrame
            dataframe to which to add the calculated values to
        index : list
            column names to use as index for calculating against
        column : str
            column from which to calculate
        
        Returns
        -------
        pd.DataFrame
            output dataframe with mean column appended
        """
        mean = data.groupby(index)[column].mean().reset_index()
        mean.rename(columns={column: 'mean'}, inplace=True)

        return pd.merge(output, mean, on=index, validate='1:1')   


    def calculate_exon_thresholds(self, data, output, index, thresholds
        ) -> pd.DataFrame:
        """
        Calculates the percent coverage at each given threshold for each region

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of from which calculate values
        output : pd.DataFrame
            dataframe to which to add the calculated values to
        index : list
            column names to use as index for calculating against
        thresholds : list
            list of thresholds at which to calculate coverage
        
        Returns
        -------
        pd.DataFrame
            output dataframe with threshold columns appended
        """
        def calculate_pct(data, threshold) -> int:
            """
            Calculates the percent coverage at the threshold for the
            given subset of data (i.e. the exon of each transcript)

            Parameters
            ----------
            data : pd.DataFrame
                subset of dataframe to calculate percent coverage at threshold
            threshold : int
                threshold to calculte coverage at

            Returns
            -------
            int
                percent coverage for the region
            """
            if isinstance(data, pd.core.frame.DataFrame):
                exon_length = int(data.iloc[0].exon_end) - \
                    int(data.iloc[0].exon_start)
                above_threshold = len(
                    data.iloc[np.where(data['cov'] > threshold)].index)

                return (above_threshold / exon_length) * 100

        for threshold in thresholds:
            # numpy vectorised call of the above calculate_pct() function
            # for much faster calculation on large dataframes vs iterating
            threshold_coverage = np.vectorize(calculate_pct, otypes=[object])(
                data.groupby(index)[['cov', 'exon_start', 'exon_end']],
                threshold,
            )

            # join numpy series back to output dataframe
            output = output.join(pd.DataFrame(
                threshold_coverage, columns=['_', f"{threshold}x"]
            )[f"{threshold}x"])
        
            output[f"{threshold}x"] = output[f"{threshold}x"].astype(float)

        return output



if __name__=="__main__":
    data = pd.read_csv(sys.argv[1], sep='\t',
                     names=["chrom", "exon_start", "exon_end", "gene", "tx", "exon",
            "cov_start", "cov_end", "cov"],
            dtype={'chrom': str, 'exon_start': int, 'exon_end': int, 'gene': str, 'tx': str, 'exon': str, 'cov_start': int, 'cov_end': int, 'cov': int})


    stats().calculate_exon_stats_parallel(data, thresholds=[100, 500, 1000, 1500, 2000])

