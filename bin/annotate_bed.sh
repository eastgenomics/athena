#!/bin/bash

input_bed=""
gene_file=""
bp_coverage=""

Help()
{
   # Display Help
    echo "
    This script may be used to perform bedtools intersect commands to generate the required annotated bed file
    for generating single sample coverage statistics."
    echo ""
    echo "Usage:"
    echo ""
    echo "-i    Input panel bed file; must have columns chromosome, start position, end position, transcript."
    echo "-g    Exons nirvana file, contains required gene and exon information."
    echo "-b    Per base coverage file (output from mosdepth or similar)."
    echo "-h    Print this Help."
    echo ""
}

# display help message on -h
while getopts ":i:g:b:h" option; do
   case $option in
        i) input_bed="$OPTARG"
        ;;
        g) gene_file="$OPTARG"
        ;;
        b) bp_coverage="$OPTARG"
        ;;
        h) # display Help
            Help
            exit 1
        ;;
        \?) # incorrect option
            echo "Error: Invalid option, please see usage below."
            Help
            exit
        ;;
        \*) # incorrect option
        ;;
   esac
done

# check for missing args
if  [ -z $input_bed ] ||
    [ -z $gene_file  ] ||
    [ -z $bp_coverage ]; then

    echo "Error: Missing arguments, please see usage below."
        Help
        exit 0
fi

# check for empty and non existent input files
for file in $input_bed $gene_file $bp_coverage; do
    [ ! -s $file ] && echo "$file does not exist or is empty. Exiting." && exit;
done

# name outfile from mosdepth coverage file
outfile=$(basename $bp_coverage)
outfile=${outfile/.per-base.bed.gz/_annotated.bed}

# add gene and exon annotation to panel bed file from exons nirvana tsv
bedtools intersect -a $input_bed -b $gene_file -wa -wb | awk 'OFS="\t" {if ($4 == $9) print}' | cut -f 1,2,3,8,9,10 > ${tmp}.txt

# add coverage annotation from per base coverage bed file
bedtools intersect -wa -wb -a $tmp.txt -b $bp_coverage | cut -f 1,2,3,4,5,6,8,9,10 > $outfile

echo "Done. Output file: " $outfile.bed
