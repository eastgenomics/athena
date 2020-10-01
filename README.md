<p align="center">
    <img height="235" width="244" src="data/static/images/logo.png">
</p>


# Athena

Athena is a tool to generate coverage statistics for NGS data, and combine these into an interactive HTML report. This gives both summary level and in depth information as to the coverage of the data, including various tables and plots to visualise the data. Examples of the output statistics files and report may be found in `data/example`.<br>


## Installation

Dependencies may be installed from the requirements.txt file using ```pip install -r requirements.txt```.
This should contain everything required to generate coverage statistics and reports. 

Tested on Ubuntu 18.04.4 and macOS 10.15.4

## Usage

It is written to take in per base coverage data (as output from tools such as mosdepth and samtools mpileup) as input to calculate coverage for target regions defined in a bed file. <br></br>

The general workflow for generating the statistics and report is as follows: <br>
- Annotate each region of the bed file with the gene, exon and per base coverage data using `annotate_bed.sh`
- Generate per exon and per gene statistics using `coverage_stats_single.py`
- Generate HTML coverage report with `coverage_report_single.py`


### Annotating BED file
The BED file containing regions of interest is first required to be annotated with gene, exon and coverage information prior to analysis. This may be done using bedtools intersect (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html), with a file containing transcript to gene and exon information, and then the per base coverage data. <br>

Included is a Bash script (`annotate_bed.sh`) to perform the required BED file annotation.

This expects the following as input:

```
-i : Input panel bed file; must have columns chromosome, start position, end position, transcript.
-g : Exons nirvana file, contains required gene and exon information.
-b : Per base coverage file (output from mosdepth or similar).
-o : Output file name prefix, will have the .bed suffix.

Example usage:

$ annotate_bed.sh -i panel_bed_file.bed -g exons_nirvana -b {input_file}.per_base.bed -o output_file
```
<br>
This wraps the bedtools intersect commands below. These commands are given as an example, the output file column ordering must match that given in `/data/example/example_annotated_bed` for calculating coverage statistics: <br>

```
$ bedtools intersect -a beds/sorted_bed_file.bed -b beds/exons_nirvana2010_no_PAR_Y_noflank.bed -wa -wb | awk 'OFS="\t" {if ($4 == $9) print}' | cut -f 1,2,3,8,9,10 > sample1_genes_exons.bed

    - sorted_bed_file.bed -- bed file defining regions of interest (columns: chromosome, start, end, transcript)
    - exons_nirvana2010_no_PAR_Y.bed -- a bed file containing transcript -> exon and gene information
    - sample1_genes_exons.bed -- bed file with added genes and exons (expected columns: chromosome, start, end, transcript, gene, exon)


$ bedtools intersect -wa -wb -a sample1_genes_exons.bed -b data/sample1.per-base.bed | cut -f 1,2,3,4,5,6,8,9,10 > sample1_gene_exon_coverage.bed

    - sample1_genes_exons.bed -- file output from above command
    - sample1_per_base.bed -- per base coverage file output from mosdepth or similar
    - sample1_gene_exon_coverage.bed -- annotated bed file ready for analysis
```


### Generating coverage statistics
The `coverage_stats_single.py` script generates both a tsv of per per exon and per gene coverage statistics. This gives a minimum, mean and maxmimum coverage for each region, along with coverage at defined thresholds. Inputs include:

```
--file: annotated bed file on which to generate report from
--build: text file with build number used for alignment, output from mosdepth (optional)
--outfile: output file name prefix, if not given the input file name will be used as the name prefix
--thresholds: threshold values to calculate coverage for as comma seperated integers (default: 10, 20, 30, 50, 100)
--flagstat: flagstat file for sample, required for generating run statistics (in development)

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
-s / --snps: VCF(s) of known SNPs to check coverage of (i.e. HGMD, ClinVar; optional)
-t / --threshold: threshold value defining sub-optimal coverage (default if not given: 20)
-n / --sample_name: optional name for title of report (gene_stats file name will be used if not given)
-o / --output: optional name for output report (sample name will be used if not given)

Example usage:

$ python3 bin/coverage_report_single.py --gene_stats output/sample1-exon-coverage_gene_stats.tsv --exon_stats output/sample1-exon-coverage_exon_stats.tsv --raw_coverage sample1_gene_exon_coverage.bed -t 30 -n sample1
```


### For development

Features to be developed:
- Generate run level statistics from multiple samples
- Generate run level report from multiple samples
- Add interactive elements to tables to increase useability (i.e sorting, filtering, searching)

Any bugs or suggestions for improvements please raise an issue.
