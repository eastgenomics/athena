import argparse
import os
import sys
import pandas as pd
import plotly as plt
import numpy as np

from jinja2 import Environment, FileSystemLoader


class singleReport():

    def load_files(self, stats, raw_coverage):
        """
        Load in raw coverage data, coverage stats file and template.
        """
        # load template
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(bin_dir, "../data/templates/")
        env = Environment(loader = FileSystemLoader(template_dir))
        template = env.get_template('single_template.html')

        # read in coverage stats file
        with open(stats) as stats_file:
            cov_stats = pd.read_csv(stats_file, sep="\t")
        
        column = [
                "chrom", "exon_start", "exon_end",
                "gene", "tx", "exon", "cov_start",
                "cov_end", "cov"
                ]

        # read in raw coverage stats file
        with open(raw_coverage) as raw_file:
            raw_coverage = pd.read_csv(raw_file, sep="\t", names=column)
        
        return cov_stats, raw_coverage


    def report_template(self, sub_20_stats):
        """
        HTML template for report
        """

        html_string = '''
        <html>
            <head>
                <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
                <style>body{ margin:0 100; background:whitesmoke; }</style>
            </head>
            <body>
                <h1>Coverage report for Twist Sample 1 (NA12878)</h1>
                <br></br>
                <h2>Exons with sub-optimal coverage</h2>
                ''' + sub_20_stats + '''
            </body>
        </html>'''

        return html_string
    

    def low_coverage_regions(self, cov_stats, raw_coverage, threshold):
        """
        Get regions where coverage at given threshold is <100%
        """
        # threshold column to check at
        threshold = str(threshold)+"x"

        column = [
                "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
                "min", "mean", "max",
                "10x", "20x", "30x", "50x", "100x"
                ]

        # empty df  
        low_stats = pd.DataFrame(columns=column)
        
        # get all exons with <100% coverage at given threshold
        for i, row in cov_stats.iterrows():
                if int(row[threshold]) < 100:
                    low_stats = low_stats.append(row, ignore_index=True)

        # pandas is terrible and forces floats, change back to int
        dtypes = {
          'chrom': int,
          'exon': int,
          'exon_start': int,
          'exon_end': int,
          'min': int,
          'max': int
        }

        low_stats = low_stats.astype(dtypes)
        
        #for idx, row in low_stats.iterrows():

        # get list of tuples of genes and exons with low coverage
        low_exon_list = low_stats.reset_index()[['gene', 'exon']].values.tolist()
        low_exon_list = [tuple(l) for l in low_exon_list]

        # get raw coverage for low coverage regions to plot
        low_raw_cov = raw_coverage[raw_coverage[['gene', 'exon']].apply(tuple, axis = 1
            ).isin(low_exon_list)].reset_index()

        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        #     print(low_raw_cov)

        sys.exit()

    
    def exon_plot(self, ):
        """
        Plot bp coverage of exon, used for those where coverage is <20x

        Args:
            -
        
        Returns:
            - 
        """





    def generate_report(self, template, cov_stats):
        """
        Generate single sample report from coverage stats

        Args:
            - template (file): template file to make report from
            - cov_stats (df): df of coverage stats for exons
        
        Returns: 

        """
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        report = os.path.join(bin_dir, "../output/", "single_report.html")

        column = [
                "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
                "min", "mean", "max",
                "10x", "20x", "30x", "50x", "100x"
                ]
          
        sub_20x = pd.DataFrame(columns=column)
        
        # get all exons with <100% coverage at 20x
        for i, row in cov_stats.iterrows():
                if int(row["20x"]) < 100:
                    sub_20x = sub_20x.append(row, ignore_index=True)

        # pandas is terrible and forces floats, change back to int
        dtypes = {
          'chrom': int,
          'exon': int,
          'exon_start': int,
          'exon_end': int,
          'min': int,
          'max': int
        }

        sub_20x = sub_20x.astype(dtypes)
        

        columns = ["min", "mean", "max", "10x", "20x", "30x", "50x", "100x"]

        stats = pd.pivot_table(cov_stats, index=["gene", "tx", "chrom", "exon", "exon_start", "exon_end"], 
                        values=["min", "mean", "max", "10x", "20x", "30x", "50x", "100x"])

        sub_20_stats = pd.pivot_table(sub_20x, index=["gene", "tx", "chrom", "exon", "exon_start", "exon_end"], 
                        values=["min", "mean", "max", "10x", "20x", "30x", "50x", "100x"])

        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        #     print(stats)


        stats = stats.reindex(columns, axis=1)
        sub_20_stats = sub_20_stats.reindex(columns, axis=1)

        stats_html = sub_20_stats.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')

        html_string = self.report_template(stats_html)

        file = open("report.html", 'w')
        file.write(html_string)
        file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate coverage report for a single sample.'
        )
    parser.add_argument('stats', help='stats file on which to generate report from')
    parser.add_argument('raw_coverage', help='raw coverage file that stats were generated from')
    parser.add_argument('--threshold', nargs='?', default=20, help="threshold to define low coverage, if not given 20 will be used as default")
    parser.add_argument('--output', help='Output file name')
    parser.add_argument('--plots', help='', nargs='?')
    args = parser.parse_args()

    # generate report
    report = singleReport()
    print(args) 

    cov_stats, raw_coverage = report.load_files(args.stats, args.raw_coverage)
    report.low_coverage_regions(cov_stats, raw_coverage, args.threshold)
    #report.generate_report(template, cov_stats)