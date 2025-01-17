<p align="center">
    <img height="250" width="250" src="athena/data/images/logo.png">
</p>


# Athena [![GitHub release][release-image]][release-url] [![made-with-python][python-image]][python-url]

Athena is a tool to generate coverage statistics for NGS data, and combine these into an interactive HTML report. This gives both summary level and in depth information as to the coverage of the data, including various tables and plots to visualise the data. Examples of the output statistics files and [report][report-link] may be found in `data/example`.


## :desktop_computer: Usage

As a minimum, a panel bed file and binned coverage data must be provided as input. The panel bed file is expected to have the following columns (with no header):

* chrom (`str`)
* start (`int`)
* end (`int`)
* gene (`str`)
* transcript (`str`)
* region (`str`)

The coverage data is expected to be binned data as output from a tool such as (mosdepth)[mosdepth-url], with the columns:

* chrom (`str`)
* start (`int`)
* end (`int`)
* depth (`int`)

The coverage data is intersected against the provided bed file using [bedtools intersect][bedtools-intersect-url] to identify the coverage in all provided regions, and then the binned coverage data unbinned to get per base data.

This data is then summarised into the coverage per transcript and per region (i.e exons / introns) against the provided `--thresholds` and `--minimum` cut off to define as low coverage. Optional plots of the low covered regions, full genes and whole chromosomes may also be included.

## :page_facing_up: Inputs




## <img src="athena/data/images/moby.png" width="34"/> Docker

A Docker image is provided for running Athena within a container. This may be built and run as follows:

```
version=$(grep -Po "([\d]\.){2}[\d]" athena/version.py)
docker build . -t athena:$version
```

```
$ docker run athena:2.0.0 athena report --help
usage: athena.py report [-h] [-r REGIONS] [-c COVERAGE]
                        [--normal_coverage NORMAL_COVERAGE]
                        [--hsmetrics HSMETRICS] [-a ANNOTATED_BED]
                        [-t THRESHOLDS [THRESHOLDS ...]] [-m MINIMUM]
                        [--panel PANEL]
                        [--clinical_indication CLINICAL_INDICATION] [-b BUILD]
                        [-o OUTPUT]
                        [--panel_filters PANEL_FILTERS [PANEL_FILTERS ...]]
                        [--summary] [--summary_file] [--limit LIMIT]
                        [--plot_sub_threshold] [--plot_chromosomes]
                        [--write_data] [--force] [--verbose]

optional arguments:
  -h, --help            show this help message and exit
  -r REGIONS, --regions REGIONS
                        Bed file of target regions to provide coverage data
                        for
  -c COVERAGE, --coverage COVERAGE
                        Bed file of raw coverage data output from samtools /
                        mosdepth
  --normal_coverage NORMAL_COVERAGE
                        tsv of previously calculated normal coverage values
                        for n samples. Requires hsmetrics file for sample
                        providing to --hsmetrics.
  --hsmetrics HSMETRICS
                        hsmetrics file for current sample, required for
                        --normal_coverage
  -a ANNOTATED_BED, --annotated_bed ANNOTATED_BED
  -t THRESHOLDS [THRESHOLDS ...], --thresholds THRESHOLDS [THRESHOLDS ...]
                        Thresholds at which to calculate percent coverage
  -m MINIMUM, --minimum MINIMUM
                        Minimum threshold value to use as cut off for defining
                        as low coverage region. Must be one of --threshold
                        values.
  --panel PANEL         Name of sequencing panel the report is for
  --clinical_indication CLINICAL_INDICATION
                        Clinical indication the report is for
  -b BUILD, --build BUILD
                        Reference build of sample data
  -o OUTPUT, --output OUTPUT
                        Prefix for naming output files. Defaults to prefix of
                        coverage bed file.
  --panel_filters PANEL_FILTERS [PANEL_FILTERS ...]
                        Preset filters of genes / transcripts to set for the
                        full gene plots, these will be presented in a drop
                        down menu for filtering the plots. These should be
                        passed as key:value pairs of panel name to display in
                        the drop down and a comma separated list of gene
                        symbols to filter with. Example: 'Cancer:BRCA1,BRCA2'
                        'Cardiac:MYH7,TNNT2'
  --summary             Display summary of genes / transcripts in report in
                        summary section
  --summary_file        Output text in summary section to a text file
  --limit LIMIT         Number of genes at which to skip full gene plot
                        generation. For large panels this significantly
                        increases the report file size.
  --plot_sub_threshold  Generates interactive plots of regions where coverage
                        is below that defined with --minimum
  --plot_chromosomes    Generates full chromosome plots of each chromosome
  --write_data          Controls if to write out the gene, region and per base
                        data to files
  --force               Force overwriting of existing files with same output
                        name
  --verbose             Increase logging verbosity to DEBUG level
```

## DNAnexus

A DNAnexus app is provided for running Athena within the DNAnexus platform. Please see the [app readme][dnanexus-readme] for details on building and running of the app.


[release-image]: https://img.shields.io/github/v/release/eastgenomics/athena
[release-url]: https://github.com/eastgenomics/athena/releases
[python-image]: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
[python-url]: https://www.python.org/

[report-link]: https://htmlpreview.github.io/?https://github.com/eastgenomics/athena/blob/master/data/example/Example_coverage_report.html

[dnanexus-readme]: https://github.com/eastgenomics/athena/blob/master/dnanexus/readme.md

[bedtools-url]: https://bedtools.readthedocs.io/en/latest/content/installation.html
[bedtools-intersect-url]: https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
[mosdepth-url]: https://github.com/brentp/mosdepth
