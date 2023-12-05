<!-- dx-header -->
# Athena (DNAnexus Platform App)

DNAnexus app Athena.

<!-- /dx-header -->

## What does this app do?
Generates coverage reports to assess quality of NGS data.
<br>

## What are typical use cases for this app?

Used to generate coverage reports to assess coverage of genomic regions defined in a bed file.
<br>

## Building the app

The source code for the app is in the `dnanexus/` directory, and may be built from the included makefile. This is will set up the app structure locally in a temporary directory with required dependencies, build the app to the platform and then clear the temp directory.

The app may be built with:
```
cd dnanexus/
make
```

## What data are required for this app to run?

### Required inputs:

- Panel BED file
- Per base coverage BED file (output from mosdepth)
- Exons nirvana (exon annotation file; originally generated from Illumina Nirvana RefSeq gff file, analogous file can be created by following instructions [here][transcript-file-url])

n.b. check [Athena][athena-url] readme for required formats of input files, specific formatting is required and if not correct will result in errors / empty files being generated.


### Optional inputs:

- thresholds: thresholds at which to calculate coverage (default: 10, 20, 30, 50, 100)
- build: text file containing reference build used for alignment (as output from [eggd_mosdepth][eggd_mosdepth-url])
- name: sample name, used to name output files and within report title. If not given this will be parsed from the per base coverage bed.
- cutoff threshold: threshold at which to define sub-optimal coverage (must be one of the threshold values; default: 20)
- snps: VCF(s) of SNPs for which to calculate coverage for (i.e. HGMD, ClinVar)
- limit: number of genes in panel at which to not generate full gene plots, for large panels this may take a long time and will result in larger output report file size.
- panel: boolean option to display panel used in report (default: True)
- summary: boolean option to include summary of genes / transcripts used in report (default: False)
- per_chromosome_coverage: boolean option to include global per-chromosome coverage plots in report (default: False)
<br>

## What does this app output?

- {sample_name}_coverage_report.html: coverage report for sample.
- {sample_name}_exon_stats.tsv: contains per exon coverage metrics, used for generating report.
- {sample_name}_gene_stats.tsv: contains per gene coverage metrics, used for generating report.
- {sample_name}_annotated.bed: raw annotated bed file, contains per base coverage data.


#### This app was made by EMEE GLH

[athena-url]: https://github.com/eastgenomics/athena
[eggd_mosdepth-url]: https://github.com/eastgenomics/eggd_mosdepth
[transcript-file-url]: https://cuhbioinformatics.atlassian.net/wiki/spaces/P/pages/2241101840/Generating+transcripts+file+for+Athena
