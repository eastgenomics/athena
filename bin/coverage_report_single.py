import argparse
import os
import sys
import pandas as pd
import numpy as np

from jinja2 import Environment, FileSystemLoader


class singleReport():

    def load_files(self, file):
        """
        Load in coverage data and template.
        """
        # load template
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(bin_dir, "../data/templates/")
        env = Environment(loader = FileSystemLoader(template_dir))
        template = env.get_template('single_template.html')


        # read in coverage stats file
        with open(file) as file:
            cov_stats = pd.read_csv(file, sep="\t")
        
        # will need to read in plots too
    
        return template, cov_stats


    def generate_report(self, template, cov_stats):
        """
        Generate single sample report from coverage stats

        Args:
            - template (file): template file to make report from
            - cov_stats (df): df of coverage stats for exons
        
        Returns: 

        """
        sample = "Twist Sample 1 (NA12878)"

        bin_dir = os.path.dirname(os.path.abspath(__file__))
        report = os.path.join(bin_dir, "../output/", "single_report.html")

        column = [
                "gene", "tx", "chrom", "exon", "exon_start", "exon_end",
                "min", "mean", "max",
                "10x", "20x", "30x", "50x", "100x"
                ]
        dtypes = {
          'chrom': int,
          'exon': int,
          'exon_start': int,
          'exon_end': int,
          'min': int,
          'max': int
        }
          
        sub_20x = pd.DataFrame(columns=column)
        

        for i, row in cov_stats.iterrows():
                if int(row["20x"]) < 100:
                    sub_20x = sub_20x.append(row, ignore_index=True)

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


        with open(report, 'w') as file:
            file.write(template.render({
                "sample": sample,
                "stats_table": stats.to_html(classes="table table-striped", justify="left"),
                "sub_20x": sub_20_stats.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')
            }))

f = open('/home/jack/report.html','w')
f.write(html_string)
f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate coverage report for a single sample.'
        )
    parser.add_argument('--stats', help='stats file on which to generate report from')
    parser.add_argument('--output', help='Output file name')
    parser.add_argument('--plots', help='', nargs='?')
    args = parser.parse_args()

    # generate report
    report = singleReport()
    template, cov_stats = report.load_files(args.stats)
    report.generate_report(template, cov_stats)