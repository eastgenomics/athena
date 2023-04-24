"""
Functions relating to generating coverage stats from the per-base data
at given thresholds
"""
from functools import partial
import multiprocessing
import re
import sys
from time import time

from natsort import natsorted
import numpy as np
import pandas as pd

from .utils import unbin


class Sample():
    """
    Generates per sample coverage values for min, mean, max and given thresholds
    """
    def calculate_exon_stats(self, coverage_data, thresholds) -> pd.DataFrame:
        """
        Calls all methods to generate coverage values for every region (i.e. exon)

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

        # unbin coverage data so each base is an individual dataframe row
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
            natsorted(set(exon_stats.exon)), ordered=True)
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
            dataframe of previously calculated per exon coverage stats
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
                columns=['gene', 'transcript']
        )

        # unbin per base coverage data for calculating mean of each transcript
        coverage_data = unbin(coverage_data)

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

                        total bases >= threshold
                        ------------------------  X 100 
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
                f"{threshold}x":pd.NamedAgg(
                    column="cov",
                    aggfunc=lambda x: sum(x >= int(threshold))
                ),
                'exon_length': pd.NamedAgg(
                    column="exon_length",
                    aggfunc=max
                )
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

        summed_thresholds = summed_thresholds.apply(lambda x: round(x, ndigits=12))

        data.drop(['exon_length'], axis=1, inplace=True)

        return pd.merge(output, summed_thresholds, on='transcript', validate='1:1')


class Run():
    """
    Generates per run coverage values with normalised per sample values.

    Normalisation is performed from hsmetrics values as the
    ON_TARGET_BASES * PCT_USABLE_BASES_ON_TARGET to get the fraction of
    aligned, du-duplicated on target bases. Each set of per exon stats is
    then multiplied against the respective hsmetric values and divided
    by the total sum of the whole runs aligned and de-duplicated on
    target bases.
    """
    def calculate_exon_stats(self, all_exon_stats) -> pd.DataFrame:
        """
        Calculate normalised run level coverage stats for each exon.

        Normalisation is done using the total usable bases of a sample
        over the total usable bases for all samples.

        Parameters
        ----------
        all_exon_stats : list(tuples)
            list of tuples of exon stats dataframe and hsmetrics dataframe,
            one tuple per sample

        Returns
        -------
        pd.DataFrame
            dataframe of normalised run level exon stats
        """
        # parse out just hsmetrics files from each pair of exon_stats -
        # hsmetrics files and calculate run level normalisation value
        all_hsmetrics = [x[1] for x in all_exon_stats]
        normalisation_value = self._calculate_total_bases(
            all_hsmetrics=all_hsmetrics
        )

        # normalise all per sample values and combine into a single df
        run_exon_stats = []

        for sample in all_exon_stats:
            run_exon_stats.append(self._normalise_sample(
                exon_stats=sample[0],
                hsmetrics=sample[1],
                normalisation_value=normalisation_value
            ))

        run_exon_stats = pd.concat(run_exon_stats).reset_index(drop=True)

        # output dataframe of each exon to add our run level stats to
        output = run_exon_stats[[
            'chrom', 'exon_start', 'exon_end', 'gene', 'transcript', 'exon'
        ]].drop_duplicates()

        # get columns we normalised per sample to aggregate per run
        calculate_columns = ['min', 'mean', 'max']
        calculate_columns.extend([
            x for x in run_exon_stats.columns if re.fullmatch(r'\d+x', x)])       

        # calculate total and std deviation for each column
        for column in calculate_columns:
            column_stats = run_exon_stats.groupby(
                ['chrom', 'exon_start', 'exon_end', 'gene', 'transcript', 'exon'],
                observed=True).agg(**{
                    column: pd.NamedAgg(
                        column=column,
                        aggfunc=sum
                    ),
                    f"{column}_std":pd.NamedAgg(
                        column=column,
                        aggfunc=np.std
                    )
                }
            ).reset_index()

            output = pd.merge(
                output, column_stats[[
                    'transcript', 'exon', column, f'{column}_std'
                ]],
                on=['transcript', 'exon'],
                validate='1:1'
            )

        # correctly sort by gene and exon
        output.exon = output.exon.cat.reorder_categories(
            natsorted(set(output.exon)), ordered=True)
        output.sort_values(by=['gene', 'transcript', 'exon'], inplace=True)

        return output


    def calculate_gene_stats(self, run_exon_stats) -> pd.DataFrame:
        """
        Calculate run level gene coverage stats from normalised run level
        per exon coverage stats.

        Parameters
        ----------
        run_exon_stats : pd.DataFrame
            dataframe of normalised run level exon stats

        Returns
        -------
        pd.DataFrame
            dataframe of normalised run level gene stats
        """
        # output dataframe of each exon to add our run level stats to
        output = run_exon_stats[['gene', 'transcript']].drop_duplicates()

        # get columns to aggregate per transcript
        calculate_columns = ['min', 'mean', 'max']
        calculate_columns.extend([
            x for x in run_exon_stats.columns if re.fullmatch(r'\d+x', x)])

        # calculate total and std deviation for each column
        for column in calculate_columns:
            column_stats = run_exon_stats.groupby(
                ['gene', 'transcript'], observed=True).agg(**{
                    column: pd.NamedAgg(
                        column=column,
                        aggfunc=np.mean
                    ),
                    f"{column}_std":pd.NamedAgg(
                        column=column,
                        aggfunc=np.std
                    )
                }).reset_index()

            output = pd.merge(
                output, column_stats[[
                    'transcript', column, f'{column}_std'
                ]],
                on=['transcript'],
                validate='1:1'
            )


        return output


    def _normalise_sample(
            self, exon_stats, hsmetrics, normalisation_value) -> pd.DataFrame:
        """
        Normalise all calculated coverage values for given sample
        against run level bases

        Parameters
        ----------
        exon_stats : pd.DataFrame
            dataframe of sample exon coverage stats
        hsmetrics : pd.DataFrame
            dataframe of sample hemetrics values
        normalisation_value : int
            run level total usable bases to normalise against

        Returns
        -------
        pd.DataFrame
            dataframe of normalised exon stats
        """
        hsmetrics = hsmetrics.astype({
            'ON_TARGET_BASES': int,
            'PCT_USABLE_BASES_ON_TARGET': float
        })

        sample_bases = (
            hsmetrics['ON_TARGET_BASES'] * hsmetrics['PCT_USABLE_BASES_ON_TARGET']
        ).iloc[0] / normalisation_value

        # columns to normalise (min, mean, max and thresholds (100x, 150x...))
        normalise_columns = ['min', 'mean', 'max']
        normalise_columns.extend([
            x for x in exon_stats.columns if re.fullmatch(r'\d+x', x)
        ])

        for column in normalise_columns:
            exon_stats[column] = exon_stats[column] * sample_bases

        return exon_stats


    def _calculate_total_bases(self, all_hsmetrics) -> int:
        """
        Calculate total usable bases from whole run of samples to use
        for normalisation

        Parameters
        ----------
        all_hsmetrics : list(pd.DataFrame)
            list of dataframes of hsmetrics files

        Returns
        -------
        int
            normalisation value to use for all coverage metrics
        """
        all_hsmetrics = pd.concat(all_hsmetrics)

        return (all_hsmetrics['ON_TARGET_BASES'] *\
                 all_hsmetrics['PCT_USABLE_BASES_ON_TARGET']).sum()
