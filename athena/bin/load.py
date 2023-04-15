import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


from .version import VERSION


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


    def filter_dtypes(self, df):
        """
        Filter list of dtypes to columns present in df
        Args:
            - df (df): dataframe to get col names from
        Returns:
            - dict: dict of col names and dtypes
        """
        return {k: v for k, v in self.dtypes.items() if k in df.columns}


    def read_panel_bed(self, bed_file):
        """
        Read in panel bed file to dataframe
        Args: bed_file (file): file handler of bed file
        Returns: panel_bed_df (df): df of panel bed file
        """
        print("Reading panel bed file")

        # first test file has expected columns
        with open(bed_file, 'r') as f:
            assert f.readline().count('\t') == 3, (
                "Provided bed file does not appear to have required 4 columns"
            )

        panel_bed = pd.read_csv(
            bed_file, sep="\t", dtype=self.dtypes, names=[
                "chrom", "start", "end", "transcript"
            ]
        )

        # strip chr from chrom if present
        panel_bed["chrom"] = panel_bed["chrom"].apply(
            lambda x: str(x).replace("chr", "")
        )

        return panel_bed


    def read_transcript_info(self, transcript_file):
        """
        Read in file that contains transcript -> gene symbol, exon no.
        Args: transcript_file (file): file handler
        Returns: transcript_info_df (df): df of transcript info
        """
        print("reading transcript information file")

        # first test file has expected columns
        with open(transcript_file, 'r') as f:
            assert f.readline().count('\t') == 5, (
                "Provided bed file does not appear to have required 4 columns"
            )

        transcript_info_df = pd.read_csv(
            transcript_file, sep="\t", dtype=self.dtypes, names=[
                "chrom", "start", "end", "gene", "transcript", "exon"
            ]
        )

        # strip chr from chrom if present
        transcript_info_df["chrom"] = transcript_info_df["chrom"].apply(
            lambda x: str(x).replace("chr", "")
        )

        return transcript_info_df


    def read_coverage_data(self, coverage_file, chunk_size=None):
        """
        Read in per-base coverage file (i.e. output of mosdepth)
        Args:
            - coverage_file (file): file handler
            - chunk_size (int): this can be 50 million + line file, use chunks
                to read in df and return an iterable to stop bedtools breaking
        Returns: pb_coverage_df (df / list): df of coverage data / list of dfs
                if chunk_size value passed
        """
        print("reading coverage file, this might take a while...")

        if chunk_size:
            # build list of dataframes split by given chunk size
            pb_coverage_df = []

            for df in pd.read_csv(
                coverage_file, sep="\t", compression="infer",
                dtype=self.dtypes,
                chunksize=chunk_size, names=["chrom", "start", "end", "cov"]
            ):
                # strip chr prefix if present
                df["chrom"] = df["chrom"].apply(
                    lambda x: str(x).replace("chr", "")
                )
                # add to list of chunk dfs
                pb_coverage_df.append(df)

            print(
                f"per-base coverage data read in to {len(pb_coverage_df)} "
                f"of {chunk_size} rows"
            )
        else:
            # read file into one df
            pb_coverage_df = pd.read_csv(
                coverage_file, sep="\t", compression="infer",
                dtype=self.dtypes,
                names=["chrom", "start", "end", "cov"]
            )

            # strip chr from chrom if present
            pb_coverage_df["chrom"] = pb_coverage_df["chrom"].apply(
                lambda x: str(x).replace("chr", "")
            )

        return pb_coverage_df


    def read_exon_stats(self, exon_stats):
        """
        Read exon stats file from coverage_single_stats into df

        Args:
            - exon_stats (file): file with per exon coverage stats
        Returns:
            - cov_stats (df): df of exon_stats file
        """
        with open(exon_stats.name) as exon_file:
            cov_stats = pd.read_csv(
                exon_file, sep="\t", comment='#', dtype=self.dtypes
            )

            # strip chr from chrom in cases of diff. formatted bed
            cov_stats["chrom"] = cov_stats["chrom"].apply(
                lambda x: str(x).replace("chr", "")
            )

        return cov_stats


    def read_gene_stats(self, gene_stats):
        """
        Read gene stats file from coverage_single_stats into df

        Args:
            - gene_stats (file): file with per gene coverage stats
        Returns:
            - cov_summary (df): df of gene_stats file
        """
        with open(gene_stats) as gene_file:
            cov_summary = pd.read_csv(
                gene_file, sep="\t", comment='#', dtype=self.dtypes
            )

        return cov_summary


    def read_raw_coverage(self, raw_coverage):
        """
        Read in raw coverage data (annotated bed file) from single stats

        Args:
            - raw_coverage (file): file of annotated bed file
        Returns:
            - raw_coverage (df): df of raw coverage
        """
        column = [
            "chrom", "exon_start", "exon_end", "gene", "transcript", "exon",
            "cov_start", "cov_end", "cov"
        ]

        # read in raw coverage stats file
        raw_coverage = pd.read_csv(
            raw_coverage, sep="\t", names=column, dtype=self.dtypes
        )
        # strip chr from chrom in cases of diff. formatted bed
        raw_coverage["chrom"] = raw_coverage["chrom"].apply(
            lambda x: str(x).replace("chr", "")
        )

        return raw_coverage


    @staticmethod
    def read_template():
        """
        Read in HTML template for report

        Args: None

        Returns:
            - html_template (str): report template
        """
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(bin_dir, "../data/templates/")
        single_template = os.path.join(template_dir, "single_template.html")

        with open(single_template, 'r') as template:
            html_template = template.read()

        return html_template


    def get_low_coverage_regions(self, cov_stats, raw_coverage, threshold):
        """
        Get regions where coverage at given threshold is <100% for
        generating low coverage plots

        Args:
            - cov_stats (df): df of coverage stats for each exon
            - raw_coverage (df): raw bp coverage for each exon
            - threshold (int): defined threshold level (default: 20)

        Returns:
            - low_raw_cov (df): df of raw bp values for each region with
                                coverage less than 100% at threshold
        """
        # get threshold columns and add to column names
        threshold_cols = list(cov_stats.filter(regex='[0-9]+x', axis=1))

        columns = [
            "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
            "min", "mean", "max"
        ]

        columns.extend(threshold_cols)

        # empty df
        # low_stats = pd.DataFrame(columns=columns)

        # get all exons with <100% coverage at given threshold
        # for i, row in cov_stats.iterrows():
        #     if int(row[threshold]) < 100:
        #         low_stats = low_stats.append(row, ignore_index=True)
        low_stats = cov_stats[cov_stats[threshold] < 100].reset_index(drop=True)

        # get list of tuples of genes and exons with low coverage to
        # select out raw coverage
        low_exon_list = low_stats[['gene', 'exon']].values.tolist()
        low_exon_list = [tuple(exon) for exon in low_exon_list]

        # get per base coverage for low coverage regions to plot
        low_raw_cov = raw_coverage[raw_coverage[['gene', 'exon']].apply(
            tuple, axis=1).isin(low_exon_list)].reset_index()

        return low_raw_cov


    @staticmethod
    def get_build_and_stats(gene_stats):
        """
        Get flagstats (if present) and reference build number used for
        alignment from header of gene_stats file.
        Args:
            - gene_stats (file): gene stats file from single_stats
        Returns:
            - flagstat (dict): flagstat values
            - build (str): HTML formatted string of reference build
        """
        flagstat = {}
        # read in flagstat and build from header of gene stats file
        with open(gene_stats) as gene_file:
            for ln in gene_file:
                if ln.startswith("#"):
                    if "build" in ln:
                        # get build number
                        reference = ln.split(":")[1]
                        # add build to string to display
                        if "37" in reference:
                            build = "<li>Reference build used for aligment<b>\
                                {}</b></li>".format(reference)
                        if "38" in build:
                            build = "<li>Reference build used for aligment<b>\
                                {}</b></li>".format(reference)
                    else:
                        # read in flagstat from header
                        key = ln.split(":")[0].strip("#")
                        val = ln.split(":")[1]
                        flagstat[key] = val

        if "build" not in locals():
            # build no. not included in gene_stats file
            build = ""

        return flagstat, build


    @staticmethod
    def get_panel_name(panel):
        """
        Get panel name from panel bed file for displaying in the report
        Args:
            - panel (file): panel bed file
        Returns:
            - panel (str): HTML formatted str for displaying in report
        """
        if panel is not None:
            panel_name = Path(panel).stem

            # format according to output of
            # https://github.com/eastgenomics/eggd_generate_bed
            panel_name = [x.strip("_") for x in panel_name.split("&&") if x]
            panel_name = [
                x.strip("_b37").strip("_b38") for x in panel_name if x
            ]
            panel_name = [x.replace("_", " ") for x in panel_name if x]
            panel_name = ",&nbsp".join(panel_name)
            panel = "<li>Panel(s) / gene(s) included in report: <b>{}</b>\
                </li>".format(panel_name)
        else:
            panel = ""

        return panel


    @staticmethod
    def get_snp_vcfs(snp_vcfs):
        """
        Get names of SNP VCFs (if used) to display in report summary

        Args:
            - snp_vcfs (array): array of file names of SNP VCFs

        Returns:
            - vcfs (str): HTML formatted string of SNP VCFs
        """
        if snp_vcfs:
            # get names of SNP vcfs used to display in report
            vcfs = ", ".join([Path(x).stem for x in snp_vcfs])
            vcfs = f"VCF(s) of known variants included in report: <b>{vcfs}</b>"
        else:
            vcfs = ""

        return vcfs


    @staticmethod
    def get_athena_ver():
        """
        Return version of Athena to display in report footer

        Args: None
        Returns:
            - version (str): version of Athena
        """
        return VERSION


    @staticmethod
    def check_threshold(threshold, cov_stats, cov_summary):
        """
        Check the given low coverage threshold is one of the threshold
        columns in both stats files

        Args:
            - threshold (int): given low threshold value
            - cov_stats (df): df of exon_stats file
            - cov_summary (df): df of gene_stats file
        Returns:
            - threshold (str): given low threshold value
        """
        if "x" not in str(threshold):
            threshold = str(threshold) + "x"

        if threshold not in list(cov_stats) and\
                threshold not in list(cov_summary):
            print("""--threshold must be one of the gene and exon
                stats coverage thresholds. Exiting now.""")
            sys.exit()

        return str(threshold)
