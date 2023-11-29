#!/bin/bash
# eggd_athena

set -exo pipefail

main() {

    # download inputs
    dx download "$panel_bed"
    dx download "$exons_file"

    # download mosdepth files
    for i in "${!mosdepth_files[@]}"
    do
        dx download "${mosdepth_files[$i]}"
    done

    # set per base and build files from downloaded mosdepth file array
    # if this is run with just per-base bed file in the array, build will
    # evaluate to an empty string and not be passed into the report script
    pb_bed=$(find . -name "*.per-base.bed.gz")
    build=$(find . -name "*.reference_build.txt")

    # download SNPs if given to SNPs dir
    mkdir snps
    for i in "${!snps[@]}"
    do
        dx download "${snps[$i]}" -o snps/
    done

    # set up bedtools
    gunzip bedtools.static.binary.gz
    mv bedtools.static.binary bedtools
    chmod a+x bedtools
    sudo mv bedtools /usr/local/bin

    # install required python packages from local packages dir
    echo "Installing python packages"
    time sudo -H python3 -m pip install --no-index --no-deps packages/*

    echo "Finished setup. Beginning analysis."

    # if sample naming given replace spaces with "_" and "/" with "-"
    if [ "$name" ]; then name=${name// /_}; fi
    if [ "$name" ]; then name=${name//\//-}; fi

    # build string of args and annotate bed file
    annotate_args="--chunk_size 20000000 --panel_bed $panel_bed_name --transcript_file $exons_file_name --coverage_file $pb_bed"
    if [ "$name" ]; then annotate_args+=" --output_name $name"; fi
    echo "Performing bed file annotation with following arguments: " $annotate_args

    time python3 athena/bin/annotate_bed.py $annotate_args
    annotated_bed=$(find athena/output/ -name "*_annotated.bed")

    # build string of inputs to pass to stats script
    stats_args=""

    if [ "$thresholds" ]; then stats_args+=" --thresholds $thresholds"; fi
    if [ "$build_name" ]; then stats_args+=" --build $build"; fi
    if [ "$name" ]; then stats_args+=" --outfile ${name}"; fi

    stats_cmd="--file $annotated_bed"
    stats_cmd+=$stats_args
    echo "Generating coverage stats with: " $stats_cmd

    # generate single sample stats
    time python3 athena/bin/coverage_stats_single.py $stats_cmd

    exon_stats=$(find athena/output/ -name "*exon_stats.tsv")
    gene_stats=$(find athena/output/ -name "*gene_stats.tsv")

    # build string of inputs for report script
    report_args=""

    if [ "$cutoff_threshold" ]; then report_args+=" --threshold $cutoff_threshold"; fi
    if [ "$name" ]; then report_args+=" --sample_name $name"; fi
    if [ "$panel" = true ]; then report_args+=" --panel $panel_bed_name"; fi
    if [ "$panel_filters" ]; then report_args+=" --panel_filters ${panel_filters} "; fi
    if [ "$summary" = true ]; then report_args+=" --summary"; fi
    if [ "${!snps[@]}" ]; then
        snp_vcfs=$(find ~/snps/ -name "*.vcf*")
        echo $snp_vcfs
        report_args+=" --snps $snp_vcfs";
    fi
    shopt -s nocasematch
    if [[ "$per_chromosome_coverage" == "true" ]]; then report_args+=" --per_base_coverage $pb_bed"; fi

    report_cmd="athena/bin/coverage_report_single.py --exon_stats $exon_stats --gene_stats $gene_stats --raw_coverage $annotated_bed --limit $limit"
    report_cmd+=$report_args
    echo "Generating report with: " $report_cmd

    # generate report
    time python3 $report_cmd

    report=$(find athena/output/ -name "*coverage_report.html")

    # compress annotated bed since it can be large
    gzip "$annotated_bed"
    annotated_bed_gz="${annotated_bed}.gz"

    echo "Completed. Uploading files"

    exon_stats=$(dx upload $exon_stats --brief)
    gene_stats=$(dx upload $gene_stats --brief)
    report=$(dx upload $report --brief)
    annotated_bed=$(dx upload $annotated_bed_gz --brief)

    dx-jobutil-add-output exon_stats "$exon_stats" --class=file
    dx-jobutil-add-output gene_stats "$gene_stats" --class=file
    dx-jobutil-add-output report "$report" --class=file
    dx-jobutil-add-output annotated_bed "$annotated_bed" --class=file
}
