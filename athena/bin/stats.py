"""
Functions relating to generating coverage stats from the per-base data
at given thresholds
"""
from functools import partial
import multiprocessing
import sys
from time import time

from natsort import natsorted
import numpy as np
import pandas as pd

# from load import loadData
from .utils import unbin


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
        exon_coverage = self._generate_empty_df(
            per_base_data=coverage_data,
            columns=['chrom', 'exon_start', 'exon_end', 'gene', 'transcript', 'exon']
        )
        coverage_data = unbin(coverage_data)

        # calculate required exon level stats
        exon_coverage = self._calculate_minimum(
            data=coverage_data,
            output=exon_coverage,
            index=['transcript', 'exon'],
            column='cov'
        )

        exon_coverage = self._calculate_mean(
            data=coverage_data,
            output=exon_coverage,
            index=['transcript', 'exon']
        )

        exon_coverage = self._calculate_maximum(
            data=coverage_data,
            output=exon_coverage,
            index=['transcript', 'exon'],
            column='cov'
        )

        exon_coverage = self._calculate_exon_thresholds(
            data=coverage_data,
            output=exon_coverage,
            thresholds=thresholds
        )

        return exon_coverage


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
        print(f"Calculating exon stats using {multiprocessing.cpu_count()} cores")
        start = time()

        # split per base dataframe into an array of dataframes by transcript
        # for passing to calculate_exon_stats() in parrallel
        transcripts = sorted(coverage_data.transcript.unique().tolist())
        transcript_array = np.array_split(np.array(transcripts), multiprocessing.cpu_count())
        split_dfs = np.asanyarray(
            [coverage_data[coverage_data["transcript"].isin(x)]
            for x in transcript_array], dtype=object
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
        exon_stats.sort_values(by=['gene', 'transcript', 'exon'], inplace=True)
        exon_stats.reset_index(drop=True, inplace=True)

        total_time = round((time() - start), 2)
        print(f"\nFinished calculating exon stats in {total_time}s")

        return exon_stats


    def calculate_gene_stats(self, coverage_data, exon_data, thresholds) -> pd.DataFrame:
        """
        Calculates per gene coverage stats from the per exon stats

        Parameters
        ----------
        coverage_data : pd.DataFrame
            dataframe of per base coverage data
        exon_data : pd.DataFrame
            dataframe of per exon coverage stats
        thresholds : list
            list of thresholds at which to calculate coverage

        Returns
        -------
        pd.DataFrame
            dataframe of per gene coverage stats
        """
        print('\nCalculating gene stats')
        start = time()

        gene_coverage = self._generate_empty_df(
            per_base_data=exon_data,
            columns=['chrom', 'gene', 'transcript']
        )

        gene_coverage = self._calculate_minimum(
            data=exon_data,
            output=gene_coverage,
            index=['transcript'],
            column='min'
        )

        gene_coverage = self._calculate_mean(
            data=coverage_data,
            output=gene_coverage,
            index=['transcript']
        )

        gene_coverage = self._calculate_maximum(
            data=exon_data,
            output=gene_coverage,
            index=['transcript'],
            column='max'
        )

        gene_coverage = self._calculate_gene_thresholds(
            data=exon_data,
            output=gene_coverage,
            thresholds=thresholds
        )

        total_time = round((time() - start), 2)
        print(f"\nFinished calculating gene stats in {total_time}s")

        return gene_coverage


    def _generate_empty_df(self, per_base_data, columns) -> pd.DataFrame:
        """
        Generate empty dataframes for adding values to, based off columns
        in given per base dataframe, keeping just one row per unique region
        from the per base rows.

        Parameters
        ----------
        per_base_data : pd.DataFrame
            dataframe of per base coverage data to extract unique exon
            and gene rows from
        columns : list
            columns to select from per_base_data
        
        Returns
        -------
        pd.DataFrame
            dataframe of unique rows for the given columns
        """
        return per_base_data[columns].drop_duplicates().reset_index(drop=True)


    def _calculate_minimum(self, data, output, index, column):
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
        minimum = data.groupby(index, as_index=False)[column].min()
        minimum.rename(columns={column: 'min'}, inplace=True)

        output = pd.merge(output, minimum, on=index, validate='1:1')

        return output


    def _calculate_maximum(self, data, output, index, column) -> pd.DataFrame:
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
        maximum = data.groupby(index, as_index=False)[column].max()
        maximum.rename(columns={column: 'max'}, inplace=True)

        return pd.merge(output, maximum, on=index, validate='1:1')


    def _calculate_mean(self, data, output, index) -> pd.DataFrame:
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
        
        Returns
        -------
        pd.DataFrame
            output dataframe with mean column appended
        """
        mean = data.groupby(index, as_index=False)['cov'].mean()
        mean.rename(columns={'cov': 'mean'}, inplace=True)

        return pd.merge(output, mean, on=index, validate='1:1')   


    def _calculate_exon_thresholds(self, data, output, thresholds
        ) -> pd.DataFrame:
        """
        Calculates the percent coverage at each given threshold for each exon.

        This is calculated as:

                        total bases > threshold
                        -----------------------  X 100 
                            exon length

        Parameters
        ----------
        data : pd.DataFrame
            dataframe of per base coverage values from which to calculate
        output : pd.DataFrame
            dataframe to which to add the calculated values to
        thresholds : list
            list of thresholds at which to calculate coverage
        
        Returns
        -------
        pd.DataFrame
            output dataframe with threshold columns appended
        """
        data['exon_length'] = data['exon_end'] - data['exon_start']

        for threshold in thresholds:
            # calculate the total bases above the given threshold
            # force exon length to be returned in output dataframe by
            # selecting the max from the column (where all will be the
            # same as this is being called on each exon)
            over = data.groupby(['transcript', 'exon']).agg(**{
                f"{threshold}x":
                pd.NamedAgg(column="cov", aggfunc=lambda x: sum(x > threshold)),
                'exon_length': pd.NamedAgg(column="exon_length", aggfunc=lambda x: max(x))
            })

            # calculate % at current threshold
            over[f"{threshold}x"] = over[f"{threshold}x"] /\
                 over['exon_length'] * 100

            # add threshold values to output dataframe
            output = pd.merge(
                output, over[[f'{threshold}x']],
                on=['transcript', 'exon'],
                validate='1:1'
            )

        return output


    def _calculate_gene_thresholds(self, data, output, thresholds) -> pd.DataFrame:
        """
        Calculate threshold values at gene level from previously calculated
        exon values, normalising against length of each exon.

        Parameters
        ----------
        data : pd.DataFrame
            per exon data to calculate from
        outut : pd.DataFrame
            per gene dataframe to add gene threshold values to
        thresholds : list
            list of thresholds at which to calculate coverage
            
        Returns
        -------
        pd.DataFrame
            dataframe of per gene values with thresholds appended
        """
        data['exon_length'] = data['exon_end'] - data['exon_start']

        thresholds = [f"{threshold}x" for threshold in thresholds]

        # calculate gene length as the sum of exon lengths for each
        # transcript, then map back onto original dataframe
        gene_length = data.groupby(['transcript'], as_index=False)['exon_length'].sum()
        gene_length.rename(columns={'exon_length': 'gene_length'}, inplace=True)
        merged = pd.merge(data, gene_length, on='transcript', validate='m:1')

        for threshold in thresholds:
            # calculate the fraction each exon contributes to
            # overall gene coverage
            merged[threshold] = merged[threshold] * merged['exon_length'] / merged['gene_length']

        # sum each of the exon coverage fractions for each transcript to
        # get total percent coverage
        summed_thresholds = merged.groupby(
            ['transcript']).agg({x: "sum" for x in thresholds})

        data.drop(['exon_length'], axis=1, inplace=True)

        return pd.merge(output, summed_thresholds, on='transcript', validate='1:1')


if __name__=="__main__":
    data = pd.read_csv(
        sys.argv[1], sep='\t',
        names=[
            "chrom", "exon_start", "exon_end", "gene", "transcript",
            "exon", "cov_start", "cov_end", "cov"
        ],
        dtype={
        'chrom': str, 'exon_start': int, 'exon_end': int, 'gene': str,
        'transcript': str, 'exon': str, 'cov_start': int,
        'cov_end': int, 'cov': int
    })


    thresholds = [10, 50, 100, 150, 200]

    exon_stats = stats().calculate_exon_stats_parallel(
        data, thresholds=thresholds)

    gene_stats = stats().calculate_gene_stats(data, exon_stats, thresholds)

    print(f"Exon stats:\n\n{exon_stats}\n")
    print(f"Gene stats:\n\n{gene_stats}\n")

