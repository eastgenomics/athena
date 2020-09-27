# eggd_coverage_report

This is a tool to generate coverage statistics for NGS data, and combine these into an interactive HTML report. This gives both summary level and in depth information as to the coverage of the data, including various tables and plots to visualise the data. Examples of the output statistics files and report may be found in `data/example`.<br>

It is written to take in per base coverage data (as output from tools such as mosdepth and samtools mpileup) as input to calculate coverage for target regions defined in a bed file. <br></br>

The general workflow for generating the statistics and report is as follows: <br>
- Annotate bed file with gene and exon information for each region
- Annotate bed file with per base coverage data
- Generate per exon and per gene statistics using coverage_stats_single.py
- Generate HTML report with coverage_report_single.py

## Installation

Dependencies may be installed from the requirements.txt file using ```pip install -r requirements.txt```.
This should contain everything required to generate coverage statistics and reports. 
Installation on macOS may have issues importing packages that have been installed, in this case use ```pip install -m``` for those packages with issues.

Tested on Ubuntu 18.04.4 and macOS 10.15.4

## Usage

### Annotating BED file
The BED file containing regions of interest is first required to be annotated with gene, exon and coverage information prior to analysis. This may be done using bedtools intersect (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html), with a file containing transcript to gene and exon information, and the per base coverage data. <br>

The following commands are given as an example, column ordering must match that given in `/data/example/example_input_coverage.txt`: <br>

- ```$ bedtools intersect -a beds/sorted_bed_file.bed -b beds/exons_nirvana2010_no_PAR_Y_noflank.bed -wa -wb | awk 'OFS="\t" {if ($4 == $9) print}' | cut -f 1,2,3,8,9,10 > sample1_genes_exons.bed```
    - sorted_bed_file.bed -- bed file defining regions of interest (expected columns: chromosome, start, end, transcript)
    - exons_nirvana2010_no_PAR_Y.bed -- a bed file containing transcript -> exon and gene information
    - sample1_genes_exons.bed -- bed file with added genes and exons (expected columns: chromosome, start, end, transcript, gene, exon)

- ```$ bedtools intersect -wa -wb -a sample1_genes_exons.bed -b data/sample1.per-base.bed | cut -f 1,2,3,4,5,6,8,9,10 > sample1_gene_exon_coverage.bed```
    - sample1_genes_exons.bed -- file output from above command
    - sample1_per_base.bed -- per base coverage file output from mosdepth or similar
    - sample1_gene_exon_coverage.bed -- annotated bed file ready for analysis

There is a bash script included that wraps the above commands - ```annotate_bed.sh```. This takes the bed file, exons_nirvana.tsv, 
per base coverage bed and an output file prefix name as input, and outputs the required file for performaing coverage calculations. <br>
Unless specified this will output the file in to ```/output```.


### Generating coverage statistics
The `coverage_stats_single.py` script generates both a tsv of per gene and per exon coverage statistics. This gives a minimum, mean and maxmimum coverage for each region, along with coverage at defined thresholds (10x, 20x, 30x, 50x, 100x). As input, this requires just the annotated bed file generated from above. It optionally also takes the reference build .txt file output from mosdepth to display in the report, and also a flagstat output file (in development). Both the per gene and per exon tsv file will be written to the output directory.

### Generating coverage reports
The `coverage_report_single.py` script generates the full HTML report. It requires several files as input (some optional):

- `-e / --exon_stats`: per exon statistics file (from `coverage_stats_single.py`)
- `-g / --gene_stats`: per gene statistics file (from `coverage_stats_single.py`)
- `-r / --raw_coverage`: annotated bed file with coverage data (generated from bedtools intersect)
- `-s / --snps`: VCF(s) of known SNPs to check coverage of (i.e. HGMD, ClinVar; optional)
- `-t / --threshold`: threshold value defining sub-optimal coverage (default if not given: 20)
- `-n / --sample_name`: optional name for title of report
- `-o / --output`: optional name for output report (sample name will be used if not given)

```$ python3 bin/coverage_report_single.py --gene_stats output/sample1-exon-coverage_gene_stats.tsv --exon_stats output/sample1-exon-coverage_exon_stats.tsv --raw_coverage sample1_gene_exon_coverage.bed -t 30 -n sample1```

### For development
Features to be developed:
- Generate run level statistics from multiple samples
- Generate run level report from multiple samples
- Add interactive elements to tables to increase useability (i.e sorting, filtering, searching)
- Add better styling to report

Any bugs or suggestions for improvements please raise an issue.
