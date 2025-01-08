#!/bin/bash

set -exo pipefail

# prefixes all lines of commands written to stdout with datetime
PS4='\000[$(date)]\011'
export TZ=Europe/London
set -exo pipefail

# set frequency of instance usage in logs to 10 seconds
kill $(ps aux | grep pcp-dstat | head -n1 | awk '{print $2}')
/usr/bin/dx-dstat 10

_set_bool_inputs() {
    : '''
    Set bool inputs to be passable as Python bool args.

    If not specified will be unset, this allows them to be implicitly skipped.
    '''
    [ "$summary" == 'true' ] && summary='--summary ' || unset summary
    [ "$summary_file" == 'true' ] && summary_file='--summary_file ' || unset summary_file
}

_set_string_inputs() {
    : '''
    Format string inputs correctly to be passable as inputs.

    If not specified will be unset, this allows them to be implicitly skipped.
    '''
    [ -n "$panel" ] && panel="--panel '${panel}' " || unset panel
    [ -n "$panel_filters" ] && panel_filters="--panel_filters ${panel_filters} " || unset panel_filters
    [ -n "$indication" ] && indication="--indication '${indication}' " || unset indication
}

_upload_outputs() {
    report=$(find . -type f -maxdepth 1 -name "*_coverage_report.html")
    gene_coverage=$(find . -type f -maxdepth 1 -name "*.gene_coverage.tsv")
    region_coverage=$(find . -type f -maxdepth 1 -name "*.region_coverage.tsv")
    annotated_bed=$(find . -type f -maxdepth 1 -name "*.coverage.bed")
    summary_text=$(find . -type f -maxdepth 1 -name "*_summary.txt")

    dx-jobutil-add-output report "$(dx upload "$report" --brief)" --class=file
    dx-jobutil-add-output gene_coverage "$(dx upload "$gene_coverage" --brief)" --class=file
    dx-jobutil-add-output region_coverage "$(dx upload "$region_coverage" --brief)" --class=file
    dx-jobutil-add-output annotated_bed "$(dx upload "$annotated_bed" --brief)" --class=file

    if [[ -n "$summary_text" ]]; then
        dx-jobutil-add-output summary_text "$(dx upload "$summary_text" --brief)" --class=file
    fi

    echo "Uploaded all output files"
}

main() {

    dx-download-all-inputs --parallel

    per_base_coverage=$(find /home/dnanexus/in/coverage_files -type f -name "*.per-base.bed.gz")
    build_file=$(find /home/dnanexus/in/coverage_files -name "*build.txt")
    [ -n "$build_file" ] && build="--build $(cat $build_file)"

    chmod a+x bedtools
    sudo mv bedtools /usr/local/bin

    echo "Installing python packages"
    time sudo -H python3 -m pip install --no-index --no-deps packages/*

    _set_bool_inputs
    _set_string_inputs

    python3 athena/athena.py \
        --regions "$regions_path" \
        --coverage "$per_base_coverage" \
        --thresholds $thresholds \
        --minimum $minimum \
        --debug \
        --force \
        $panel \
        $indication \
        $build \
        $panel_filters \
        $summary \
        $summary_file

    _upload_outputs
}
