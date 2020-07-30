"""
Script to generate single sample coverage report.
Takes single sample coverage stats as input, along with the raw
coverage input file and an optional "low" coverage threshold (default 20).

Jethro Rainford 200722
"""

import argparse
import os
import sys
import tempfile
import pandas as pd
import plotly.tools as plotly_tools
import plotly
import plotly.express as px
import matplotlib.pyplot as plt
from plotly.offline import plot
import plotly.graph_objs as go
import numpy as np
import math

from jinja2 import Environment, FileSystemLoader
from plotly.graph_objs import *


class singleReport():

    def load_files(self, exon_stats, gene_stats, raw_coverage):
        """
        Load in raw coverage data, coverage stats file and template.
        """
        # load template
        bin_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(bin_dir, "../data/templates/")
        env = Environment(loader = FileSystemLoader(template_dir))
        template = env.get_template('single_template.html')

        # read in exon stats file
        with open(exon_stats) as exon_file:
            cov_stats = pd.read_csv(exon_file, sep="\t")
        
        # read in gene stats file
        with open(gene_stats) as gene_file:
            cov_summary = pd.read_csv(gene_file, sep="\t")


        column = [
                "chrom", "exon_start", "exon_end",
                "gene", "tx", "exon", "cov_start",
                "cov_end", "cov"
                ]
        
        # read in raw coverage stats file
        with open(raw_coverage) as raw_file:
            raw_coverage = pd.read_csv(raw_file, sep="\t", names=column)
        
        return cov_stats, cov_summary, raw_coverage


    def report_template(self, total_stats, gene_stats, sub_20_stats, fig, report_vals):
        """
        HTML template for report
        """
        print(type(total_stats))
        html_string = '''
        <html>
            <head>

                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <style>
                    body{
                        margin-top: 100;
                        margin-bottom: 150;
                        width: 75%;
                        margin-left: auto;
                        margin-right: auto
                        }
                    
                    tr:hover {background-color:#ecebfc !important}

                    # .collapsible {
                    #     background-color: #777;
                    #     color: white;
                    #     cursor: pointer;
                    #     padding: 18px;
                    #     width: 100%;
                    #     border: none;
                    #     text-align: left;
                    #     outline: none;
                    #     font-size: 15px;
                    #     }

                    .active, .collapsible:hover {
                        background-color: #555;
                        }

                    .collapsible:after {
                        content: '+';
                        color: white;
                        font-weight: bold;
                        float: right;
                        margin-left: 5px;
                        }

                    .active:after {
                        content: "-";
                        }
                    .content {
                        padding: 0 18px;
                        max-height: 0;
                        overflow: hidden;
                        transition: max-height 0.5s ease-out;
                        background-color: #f1f1f1;
                        }              
                </style>

            <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.4.4/jquery.js"></script>
            <script type="text/javascript">


            </script>

            </head>
            <body>
            <div class="well">
                <h1>Coverage report for: Twist Sample 1 (NA12878)</h1>
                <br></br>
                <h2>Summary</h2>
                <br>
                <p style="font-size:18px">
                The following report provides an assessment of coverage for the given sample.<br></br>
                It contains the following sections: <br></br>
                <ul>
                    <li> A table of exons with less than 100% coverage at '''+report_vals["threshold"]+'''x.</li><br>
                    <li> A series of interactive plots of exons with sub-optimal coverage.</li><br>
                    <li> A summary table of coverage across all genes.</li><br>
                    <li> A full table of per exon coverage across all genes.</li>
                </ul>    
                <br>
                </p>
                <p style="font-size:18px">

                Of the '''+report_vals["total_genes"]+''' genes included in the panel, 
                '''+report_vals["exon_issues"]+''' exons in 
                '''+report_vals["gene_issues"]+''' genes had sub optimal
                coverage. To assess what portion of the exon(s) have sub-optimal coverage,
                the plots display an interactive window on hovering, showing the coverage
                at that current base.<br>

                </p>
                <br></br>

                <h2>Exons with sub-optimal coverage</h2>
                
                <table id=sub>
                
                data-search="true"
                data-show-refresh="true"
                data-show-toggle="true"
                data-show-fullscreen="true"
                data-show-columns="true"
                data-show-columns-toggle-all="true"
                data-detail-view="true"
                data-show-export="true"

                    ''' + sub_20_stats + '''
                </table>


                <br></br>
                '''+ fig +'''
                <br></br><br></br>

                <h2> Per gene coverage summary </h2>
                
                <button type="button" class="collapsible">Show / hide table</button>
                <div class="content">
                    <table>
                       ''' + gene_stats + '''
                    </table>
                </div>

                <br></br><br></br>

                <h2> Coverage for all regions of all genes </h2>

                <button type="button" class="collapsible">Show / hide table</button>
                <div class="content">
                    <table>
                        ''' + total_stats + '''
                    </table>
                </div>

            <script>
            var coll = document.getElementsByClassName("collapsible");
            var i;

            for (i = 0; i < coll.length; i++) {
            coll[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var content = this.nextElementSibling;
                if (content.style.maxHeight){
                content.style.maxHeight = null;
                } else {
                content.style.maxHeight = content.scrollHeight + "px";
                } 
            });
            }
            </script>

            </body>
            </div>
        </html>'''

        return html_string
    

    def colour_row(self, row):
        
        color = 'background-color: {}'.format('green' if row["20x"] > 90 else 'red')
        
        #return (color, color, color, color, color, color, color, color, color, color, color, color, color, color)
        
        #print(row["20x"])

        # if row["20x"] >= 95:
        #     bar = "bar(color='#2E7F18', vmin=0,vmax=100)"

        # if  95 > row["20x"] >= 90:
        #     bar = "bar(color='#45731E', vmin=0,vmax=100)"

        # if  90 > row["20x"] >= 70:
        #     bar = "bar(color='#675E24', vmin=0,vmax=100)"

        # if  70 > row["20x"] >= 50:
        #     bar = "bar(color='#8D472B', vmin=0,vmax=100)"

        # if  50 > row["20x"] >= 30:
        #     bar = "bar(color='#B13433', vmin=0,vmax=100)"

        # if  30 > row["20x"] >= 10:
        #     bar = "bar(color='#C82538', vmin=0,vmax=100)"

        # if  10 > row["20x"]:
        #     bar = "bar(color='grey', vmin=0,vmax=100)"




        return bar



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
        
        # get list of tuples of genes and exons with low coverage to select out raw coverage
        low_exon_list = low_stats.reset_index()[['gene', 'exon']].values.tolist()
        low_exon_list = [tuple(l) for l in low_exon_list]

        # get raw coverage for low coverage regions to plot
        low_raw_cov = raw_coverage[raw_coverage[['gene', 'exon']].apply(tuple, axis = 1
            ).isin(low_exon_list)].reset_index()

        print(low_raw_cov)

        return low_raw_cov
    

    def low_exon_plot(self, low_raw_cov, threshold):
        """
        Plot bp coverage of exon, used for those where coverage is <20x

        Args:
            - low_raw_cov (df): df of raw coverage for exons with low coverage
        
        Returns:
            - 
        """
        # get list of tuples of genes and exons to define plots
        genes = low_raw_cov.drop_duplicates(["gene", "exon"])[["gene", "exon"]].values.tolist()
        genes = [tuple(l) for l in genes]

        # sort list of genes/exons by gene and exon
        genes = sorted(genes, key=lambda element: (element[0], element[1]))

        plot_titles = [str(x[0])+" exon: "+str(x[1]) for x in genes]

        print(plot_titles)

        low_raw_cov["exon_len"] = low_raw_cov["exon_end"] - low_raw_cov["exon_start"]
        #low_raw_cov["label_name"] = low_raw_cov["gene"]+" exon: "+(low_raw_cov["exon"].astype(str))

        low_raw_cov["relative_position"] = low_raw_cov["exon_end"] - round(((low_raw_cov["cov_end"] + low_raw_cov["cov_start"])/2))
        
        # highest coverage value to set y axis for all plots
        max_y = max(low_raw_cov["cov"].tolist())

        # set no. rows to number of plots / number of columns to define grid
        columns = 4
        rows = math.ceil(len(genes)/4)

        # define grid to add plots to
        fig = plotly_tools.make_subplots(
                            rows=rows, cols=columns, print_grid=True, 
                            horizontal_spacing= 0.06, vertical_spacing= 0.06, 
                            subplot_titles=plot_titles)

        plots = []
        plot_titles = []
        
        # counter for grid
        row_no = 1
        col_no = 1

        for gene in genes:
            # make plot for each gene / exon
            print(gene)

            title = low_raw_cov.loc[(low_raw_cov["gene"] == gene[0]) & (low_raw_cov["exon"] == gene[1])]         
            
            # counter for grid, by gets to 5th entry starts new row
            if row_no // 5 == 1:
                col_no += 1
                row_no = 1

            exon_cov = low_raw_cov.loc[(low_raw_cov["gene"] == gene[0]) & (low_raw_cov["exon"] == gene[1])]

            # built list of threshold points to plot line
            yval = [threshold]*max_y

            # generate plot and threshold line to display

            # if any plots have no coverage, just display empty plot            
            if sum(exon_cov["cov"]) != 0:
                plot = go.Scatter(
                            x=exon_cov["cov_start"], y=exon_cov["cov"],
                            mode="lines",
                            hovertemplate = '<i>position: </i>%{x}'+ '<br>coverage: %{y}<br>',
                            )   
            else:
                plot = go.Scatter(
                                x=exon_cov["cov_start"], y=exon_cov["cov"],
                                mode="markers", marker={"opacity":0}
                                )


            threshold_line = go.Scatter(x=exon_cov["cov_start"], y=yval, hoverinfo='skip', 
                            mode="lines", line = dict(color = 'rgb(205, 12, 24)', 
                            width = 1))
                        
            plots.append(plot)            

            # add to subplot grid
            print(col_no, row_no)

            fig.add_trace(plot, col_no, row_no)
            fig.add_trace(threshold_line, col_no, row_no)


            row_no = row_no + 1

        fig["layout"].update(height=1750, showlegend=False)
              
                
        fig.update_xaxes(nticks=3, ticks="", showgrid=True, tickformat=',d')
        fig.update_yaxes(title='coverage')    
        fig.update_xaxes(title='exon position', color='#FFFFFF')    

        plotly.io.write_html(fig, "plots.html")

        fig = fig.to_html(full_html=False)

        return fig


    def generate_report(self, cov_stats, cov_summary, fig, threshold):
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
        
        # get all exons with <100% coverage at threshold 
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

        total_stats = pd.pivot_table(cov_stats, index=["gene", "tx", "chrom", "exon", "exon_start", "exon_end"], 
                        values=["min", "mean", "max", "10x", "20x", "30x", "50x", "100x"])

        sub_20_stats = pd.pivot_table(sub_20x, index=["gene", "tx", "chrom", "exon", "exon_start", "exon_end"], 
                        values=["min", "mean", "max", "10x", "20x", "30x", "50x", "100x"])

        # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        #     print(stats)


        total_stats = total_stats.reindex(columns, axis=1)
        sub_20_stats = sub_20_stats.reindex(columns, axis=1)

        total_stats.reset_index(inplace=True)
        sub_20_stats.reset_index(inplace=True)

        # get values for report
        total_genes = len(cov_summary["gene"])
        gene_issues = len(list(set(sub_20_stats["gene"].tolist())))
        exon_issues = len(sub_20_stats["exon"])

        # empty dict to add values for displaying in report text
        report_vals = {}

        report_vals["total_genes"] = str(total_genes)
        report_vals["gene_issues"] = str(gene_issues)
        report_vals["threshold"] = str(threshold)
        report_vals["exon_issues"] = str(exon_issues)

        sub_20_stats['20x'] = sub_20_stats['20x'].apply(lambda x: int(x))

        # set ranges for colouring cells
        x0 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 10) & (sub_20_stats['20x'] > 0)].index, '20x']
        x10 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 30) & (sub_20_stats['20x'] >= 10)].index, '20x']
        x30 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 50) & (sub_20_stats['20x'] >= 30)].index, '20x']
        x50 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 70) & (sub_20_stats['20x'] >= 50)].index, '20x']
        x70 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 90) & (sub_20_stats['20x'] >= 70)].index, '20x']
        x90 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] < 95) & (sub_20_stats['20x'] >= 90)].index, '20x']
        x95 = pd.IndexSlice[sub_20_stats.loc[(sub_20_stats['20x'] >= 95)].index, '20x']
        

        # apply colours to coverage cell based on value, 0 is given solid red
        s = sub_20_stats.style.apply(
            lambda x: ["background-color: #d70000" if x["20x"] == 0 and idx==10 else "" for idx,v in enumerate(x)], axis=1)\
            .bar(subset=x0, color='red', vmin=0, vmax=100)\
            .bar(subset=x10, color='#990000', vmin=0, vmax=100)\
            .bar(subset=x30, color='#C82538', vmin=0, vmax=100)\
            .bar(subset=x50, color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=x70, color='#FF4500', vmin=0, vmax=100)\
            .bar(subset=x90, color='#45731E', vmin=0, vmax=100)\
            .bar(subset=x95, color='#007600', vmin=0, vmax=100)\
            .set_table_attributes('table border="1" class="dataframe table table-hover table-bordered"')\

        # generate html string from table objects
        gene_stats = cov_summary.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')
        total_stats = total_stats.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')
        sub_20_stats = s.render()

        #sub_20_stats = sub_20_stats.to_html().replace('<table border="1" class="dataframe">','<table class="table table-striped">')


        # add tables & plots to template
        html_string = self.report_template(total_stats, gene_stats, sub_20_stats, fig, report_vals)

        file = open("coverage_report.html", 'w')
        file.write(html_string)
        file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate coverage report for a single sample.'
        )
    parser.add_argument(
        'exon_stats', help='exon stats file (from coverage_stats_single.py)')
    parser.add_argument(
        'gene_stats', help='gene stats file (from coverage_stats_single.py)')
    parser.add_argument(
        'raw_coverage', help='raw coverage file that stats were generated from')
    parser.add_argument(
        '--threshold', nargs='?', default=20, help="threshold to define low coverage (int), if not given 20 will be used as default. Must be one of the thresholds in the input file.")
    parser.add_argument(
        '--output', help='Output file name')

    args = parser.parse_args()

    # initialise
    report = singleReport()

    # read in files
    cov_stats, cov_summary, raw_coverage = report.load_files(
                                                        args.exon_stats, 
                                                        args.gene_stats, 
                                                        args.raw_coverage
                                                        )
    
    # get regions with low coverage
    low_raw_cov = report.low_coverage_regions(cov_stats, raw_coverage, args.threshold)
    
    # generate plot of sub optimal regions
    fig = report.low_exon_plot(low_raw_cov, args.threshold)
     
    # generate report
    report.generate_report(cov_stats, cov_summary, fig, args.threshold)