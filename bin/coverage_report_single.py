import argparse
import os
import sys
import pandas as pd

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
        
    
        return template, cov_stats


    def generate_report(self, template, cov_stats):
        """
        Generate single sample report from coverage stats
        """





if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate coverage report for a single sample.'
        )
    parser.add_argument('--file', help='bed file on which to generate report from')
    parser.add_argument('--output', help='Output file name')
    parser.add_argument('--plots', help='', nargs='?')
    args = parser.parse_args()

    # generate report
    report = singleReport()
    template, cov_stats = report.load_files(args.file)
    report.generate_report(template, cov_stats)