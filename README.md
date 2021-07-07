<p align="center">
    <img height="250" width="250" src="data/static/images/logo.png">
</p>


# Athena [![GitHub release][release-image]][release-url] [![made-with-python][python-image]][python-url]


Athena is a tool to generate coverage statistics for NGS data, and combine these into an interactive HTML report. This gives both summary level and in depth information as to the coverage of the data, including various tables and plots to visualise the data. Examples of the output statistics files and [report][report-link] may be found in `data/example`.<br>


## Installation

Dependencies may be installed from the requirements.txt file using ```pip install -r requirements.txt```.
This should contains all the required python packages required to generate coverage statistics and reports.
In addition, [BEDtools][bedtools-url] is also required to be installed and on path.

Tested on Ubuntu 18.04.4 and macOS 10.15.4

## Usage

It is written to take in per base coverage data (as output from tools such as [mosdepth][mosdepth-url]) as input to calculate coverage for target regions defined in a bed file. <br></br>

The general workflow for generating the statistics and report is as follows: <br>
- Annotate each region of the bed file with the gene, exon and per base coverage data using `annotate_bed.py`
- Generate per exon and per gene statistics using `coverage_stats_single.py`
- Generate HTML coverage report with `coverage_report_single.py`

For DNAnexus cloud platform users, an Athena [dx applet][dx-url] has also been built.


### Expected file formats

As a minimum, Athena requires 3 input files. These are a bed file for the gene panel, a file of transcript information and the output of your coverage tool (mosdepth, samtools etc.). These files MUST have the following columns:

- panel bed file: `chromosome  start  end  transcript`
- transcript file: `chromosome  start  end  gene  transcript  exon`
- coverage file: `chromosome  start  end  coverage`

n.b. the process for creating the transcript file may be found [here][transcript-file-url].

### Annotating BED file
The BED file containing regions of interest is first required to be annotated with gene, exon and coverage information prior to analysis. This may be done using [BEDtools intersect][bedtools-intersect-url], with a file containing transcript to gene and exon information, and then the per base coverage data. Currently, 100% overlap is required between coordinates in the panel bed file and the transcript annotation file, therefore you must ensure any added flank regions etc. are the same.<br>

Included is a Python script (`annotate_bed.py`) to perform the required BED file annotation.

Expected inputs:

```
-p / --panel_bed : Input panel bed file; must have ONLY the following 4 columns chromosome, start position, end position, gene/transcript

-t / --transcript_file : Transcript annotation file, contains required gene and exon information. Must have ONLY the following 6 columns:
chromosome, start, end, gene, transcript, exon

-c / --coverage_file : Per base coverage file (output from mosdepth or similar)

-s / -chunk_size : (optional) nrows to split per-base coverage file for intersecting, with <16GB RAM can lead to bedtools intersect failing. Reccomended values: 16GB RAM -> 20000000; 8GB RAM -> 10000000

-n / --output_name : (optional) Prefix for naming output file, if not given will use name from per base coverage file

Example usage:

$ annotate_bed.py -p panel_bed_file.bed -t transcript_file.tsv -c {input_file}.per_base.bed
```
<br>
This wraps the bedtools intersect commands below. These commands are given as an example, the output file column ordering must match that given in /data/example example_annotated_bed for calculating coverage statistics:
<br>

```
$ bedtools intersect -a beds/sorted_bed_file.bed -b beds/exons_nirvana2010_no_PAR_Y_noflank.bed -wa -wb -f 1.0 -r | awk 'OFS="\t" {if ($4 == $9) print}' | cut -f 1,2,3,8,9,10 > sample1_genes_exons.bed

    - sorted_bed_file.bed -- bed file defining regions of interest (columns: chromosome, start, end, transcript)
    - exons_nirvana2010_no_PAR_Y.bed -- a bed file containing transcript -> exon and gene information
    - sample1_genes_exons.bed -- bed file with added genes and exons (expected columns: chromosome, start, end, transcript, gene, exon)


$ bedtools intersect -wa -wb -a sample1_genes_exons.bed -b data/sample1.per-base.bed | cut -f 1,2,3,4,5,6,8,9,10 > sample1_gene_exon_coverage.bed

    - sample1_genes_exons.bed -- file output from above command
    - sample1_per_base.bed -- per base coverage file output from mosdepth or similar
    - sample1_gene_exon_coverage.bed -- annotated bed file ready for analysis
```


### Generating coverage statistics
The `coverage_stats_single.py` script generates both a tsv of per per exon and per gene coverage statistics. This gives a minimum, mean and maxmimum coverage for each region, along with coverage at defined thresholds.

Expected inputs:

```
--file: annotated bed file on which to generate report from
--build: text file with build number used for alignment, output from mosdepth (optional)
--outfile: output file name prefix, if not given the input file name will be used as the name prefix
--thresholds: threshold values to calculate coverage for as comma seperated integers (default: 10, 20, 30, 50, 100)
--flagstat: flagstat file for sample, required for generating run statistics (in development)
--cores: Number of CPU cores to utilise, for larger numbers of genes this will drastically reduce run time. If not given will use maximum available

Example usage:

$ python3 bin/coverage_stats_single.py  --file annotated_bed_file --build {sample}_reference_build.txt --thresholds 20, 40, 60, 80 --outfile example_sample
```

Example output files are given in `/data/example/`


### Generating coverage reports
The `coverage_report_single.py` script generates the full HTML report. It requires several files as input (some optional):

```
-e / --exon_stats: per exon statistics file (from `coverage_stats_single.py`)
-g / --gene_stats: per gene statistics file (from `coverage_stats_single.py`)
-r / --raw_coverage: annotated bed file with coverage data (generated from annotate_bed.sh / bedtools intersect)
-s / --snps: VCF(s) of known SNPs to check coverage of (optional; i.e. HGMD, ClinVar)
-t / --threshold: threshold value defining sub-optimal coverage (optional; default if not given: 20)
-n / --sample_name: name for title of report (optional; gene_stats file name will be used if not given)
-o / --output: name for output report (optional; sample name will be used if not given)
-p / --panel: panel bed file used for initial annotation, name will be displayed in summary of report (optional)
-l / --limit: number of genes at which to limit including full gene plots, large numbers of genes may take a long time to generate the plots (optional)
-m / --summary: boolean flag to add clinical report summary text in summary section, includes list of all genes with transcripts (optional; default False)
--cores: Number of CPU cores to utilise, for larger numbers of genes this will drastically reduce run time. If not given will use maximum available

Example usage:

$ python3 bin/coverage_report_single.py --gene_stats output/sample1-exon-coverage_gene_stats.tsv --exon_stats output/sample1-exon-coverage_exon_stats.tsv --raw_coverage sample1_gene_exon_coverage.bed -t 30 -n sample1
```


Any bugs or suggestions for improvements please raise an issue.


[release-image]: https://img.shields.io/github/v/release/eastgenomics/athena
[release-url]: https://github.com/eastgenomics/athena/releases
[python-image]: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
[python-url]: https://www.python.org/

[report-link]: https://htmlpreview.github.io/?https://github.com/eastgenomics/athena/blob/master/data/example/Example_coverage_report.html

[bedtools-url]: https://bedtools.readthedocs.io/en/latest/content/installation.html
[bedtools-intersect-url]: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
[mosdepth-url]: https://github.com/brentp/mosdepth

[dx-url]: https://github.com/eastgenomics/eggd_athena
[transcript-file-url]: https://cuhbioinformatics.atlassian.net/wiki/spaces/P/pages/2241101840/Generating+transcripts+file+for+Athena
