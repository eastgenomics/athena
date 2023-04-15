    parser = argparse.ArgumentParser(
        description='Generate coverage stats and HTML report.'
    )

    parser.add_argument(
        '--panel_bed', '-p',
        help='panel bed file'
    )
    parser.add_argument(
        '--transcript_file', '-t',
        help='file with gene and exon information'
    )
    parser.add_argument(
        '--coverage_file', '-c',
        help='per base coverage data file'
    )
    parser.add_argument(
        '--chunk_size', '-s', type=int,
        help='number lines to read per-base coverage file in one go'
    )
    parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )
    parser.add_argument(
        '-l', '--limit', nargs='?',
        help="Number of genes at which to limit including full gene plots,\
        large numbers of genes takes a long time to generate the plots.",
        default=-1,
        required=False
    )
    parser.add_argument(
        '-m', '--summary',
        help="If passed, a short paragraph will be included in the\
        summary section. This includes details on the sequencing and the\
        genes/transcripts used in the panel.",
        default=False, action='store_true'
    )
    parser.add_argument(
        '-n', '--sample_name', nargs='?',
        help="Name of sample to display in report, if not\
            specified this will be the prefix of the\
            gene_stats input file.",
        required=False
    )
    parser.add_argument(
        '-o', '--output', nargs='?',
        help='name preifx for output file, if none will use coverage file',
        required=False
    )
    parser.add_argument(
        '--panel', nargs='?',
        help='(Optional) Panel bed file used from annotation, if passed\
        name of file will be displayed in report to show what\
        panel(s) / gene(s) were included.',
        required=False
    )
    parser.add_argument(
        '--thresholds', nargs='*',
        default=[10, 20, 30, 50, 100],
        help='threshold values to calculate coverage for as comma\
            seperated integers (default: 10, 20, 30, 50, 100).'
    )
    parser.add_argument(
        '--sub_threshold', nargs='?',
        default=20, type=int,
        help="threshold to define low coverage (int). Must be one defined \
            by --thresholds",
        required=False
    )

    # sub commands for running seperate parts
    subparsers = parser.add_subparsers(
        dest='development', help='sub-commands for development and testing'
    )

    # parser to annotate panel bed file
    annotate_parser = subparsers.add_parser(
        'annotate', help='annotate bed file with transcript & coverage info, \
            requires panel bed, transcript and coverage args'
    )
    annotate_parser.add_argument(
        '--panel_bed', '-p', required=True,
        help='panel bed file'
    )
    annotate_parser.add_argument(
        '--transcript_file', '-t', required=True,
        help='file with gene and exon information'
    )
    annotate_parser.add_argument(
        '--coverage_file', '-c', required=True,
        help='per base coverage data file'
    )
    annotate_parser.add_argument(
        '--chunk_size', '-s', type=int,
        help='number lines to read per-base coverage file in one go'
    )
    annotate_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )

    # parser to generate coverage stats
    stats_parser = subparsers.add_parser(
        'stats', help='generate per gene and per exon coverage stats'
    )
    stats_parser.add_argument(
        'annotated_bed',
        help='annotated bed file from annotated_bed'
    )
    stats_parser.add_argument(
        '--flagstat', nargs='?',
        help='Optional flagstat file, needed for generating run stats.'
    )
    stats_parser.add_argument(
        '--build', nargs='?',
        help='Optional text file with build number used for alignment.'
    )
    stats_parser.add_argument(
        '--outfile', nargs='?', help='Output file name prefix, if not\
        given the input file name will be used as the name prefix.',
        type=str
    )
    stats_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )


    # parser to generate report
    report_parser = subparsers.add_parser(
        'report', help='generate report from gene and exon stats'
    )
    report_parser.add_argument(
        'annotated_bed',
        help='annotated bed file from annotated_bed'
    )
    report_parser.add_argument(
        'exon_stats',
        help='exon stats file (output/{sample_id}_exon_stats.tsv)'
    )
    report_parser.add_argument(
        'gene_stats',
        help='gene stats file (output/{sample_id}_gene_stats.tsv)'
    )
    report_parser.add_argument(
        '-s', '--snps', nargs='*',
        help='Optional; check coverage of VCF(s) of SNPs.'
    )
    report_parser.add_argument(
        '-t', '--threshold', nargs='?',
        default=20, type=int,
        help="threshold to define low coverage (int), if not\
            given 20 will be used as default. Must be one of\
            the thresholds in the input file.",
        required=False
    )
    report_parser.add_argument(
        '-n', '--sample_name', nargs='?',
        help="Name of sample to display in report, if not\
            specified this will be the prefix of the\
            gene_stats input file.",
        required=False
    )
    report_parser.add_argument(
        '-o', '--output', nargs='?',
        help='Output report name, if not specified the sample\
        name from the report will be used.',
        required=False
    )
    report_parser.add_argument(
        '-p', '--panel', nargs='?',
        help='(Optional) Panel bed file used from annotation, if passed\
        name of file will be displayed in report to show what\
        panel(s) / gene(s) were included.',
        required=False
    )
    report_parser.add_argument(
        '-l', '--limit', nargs='?',
        help="Number of genes at which to limit including full gene plots,\
        large numbers of genes takes a long time to generate the plots.",
        default=-1,
        required=False
    )
    report_parser.add_argument(
        '-m', '--summary',
        help="If passed, a short paragraph will be included in the\
        summary section. This includes details on the sequencing and the\
        genes/transcripts used in the panel.",
        default=False, action='store_true'
    )
    report_parser.add_argument(
        '--cores', nargs='?', default=None,
        help='Number of cores to utilise, for larger numbers of genes this\
        will drastically reduce run time. If not given will use maximum\
        available'
    )
