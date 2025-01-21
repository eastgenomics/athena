"""Test data for test_io.TestReadHsmetrics"""

import os

import polars as pl
import pytest


@pytest.fixture
def input_hsmetrics_file(tmp_path):
    """Example hsmetrics file contents"""
    target_output = os.path.join(tmp_path, "hsmetrics.tsv")

    with open(target_output, "w+") as fh:
        fh.write(
            "## htsjdk.samtools.metrics.StringHeader$\n"
            "# CollectHsMetrics BAIT_INTERVALS=[targets.picard]"
            " TARGET_INTERVALS=[targets.picard]"
            " INPUT=/home/dnanexus/in/sorted_bam/133830882-24326R0034-24NGCEN83-9527-F-99347387_markdup.bam"
            " OUTPUT=/home/dnanexus/out/eggd_picard_stats/QC/133830882-24326R0034-24NGCEN83-9527-F-99347387_markdup.hsmetrics.tsv"
            " PER_TARGET_COVERAGE=/home/dnanexus/out/eggd_picard_stats/QC/133830882-24326R0034-24NGCEN83-9527-F-99347387_markdup.pertarget_coverage.tsv"
            " COVERAGE_CAP=100000 REFERENCE_SEQUENCE=genome.fa   "
            " METRIC_ACCUMULATION_LEVEL=[ALL_READS] NEAR_DISTANCE=250"
            " MINIMUM_MAPPING_QUALITY=20 MINIMUM_BASE_QUALITY=20"
            " CLIP_OVERLAPPING_READS=true INCLUDE_INDELS=false"
            " SAMPLE_SIZE=10000 ALLELE_FRACTION=[0.001, 0.005, 0.01, 0.02,"
            " 0.05, 0.1, 0.2, 0.3, 0.5] VERBOSITY=INFO QUIET=false"
            " VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5"
            " MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false"
            " CREATE_MD5_FILE=false GA4GH_CLIENT_SECRETS=client_secrets.json"
            " USE_JDK_DEFLATER=false USE_JDK_INFLATER=false$\n"
            "## htsjdk.samtools.metrics.StringHeader$\n"
            "# Started on: Fri Dec 20 18:52:28 UTC 2024$\n"
            "$\n"
            "## METRICS CLASS\tpicard.analysis.directed.HsMetrics$\n"
            "BAIT_SET\tBAIT_TERRITORY\tBAIT_DESIGN_EFFICIENCY\tON_BAIT_BASES\tNEAR_BAIT_BASES\tOFF_BAIT_BASES\tPCT_SELECTED_BASES\tPCT_OFF_BAIT\tON_BAIT_VS_SELECTED\tMEAN_BAIT_COVERAGE\tPCT_USABLE_BASES_ON_BAIT\tPCT_USABLE_BASES_ON_TARGET\tFOLD_ENRICHMENT\tHS_LIBRARY_SIZE\tHS_PENALTY_10X\tHS_PENALTY_20X\tHS_PENALTY_30X\tHS_PENALTY_40X\tHS_PENALTY_50X\tHS_PENALTY_100X\tTARGET_TERRITORY\tGENOME_SIZE\tTOTAL_READS\tPF_READS\tPF_BASES\tPF_UNIQUE_READS\tPF_UQ_READS_ALIGNED\tPF_BASES_ALIGNED\tPF_UQ_BASES_ALIGNED\tON_TARGET_BASES\tPCT_PF_READS\tPCT_PF_UQ_READS\tPCT_PF_UQ_READS_ALIGNED\tMEAN_TARGET_COVERAGE\tMEDIAN_TARGET_COVERAGE\tMAX_TARGET_COVERAGE\tMIN_TARGET_COVERAGE\tZERO_CVG_TARGETS_PCT\tPCT_EXC_DUPE\tPCT_EXC_ADAPTER\tPCT_EXC_MAPQ\tPCT_EXC_BASEQ\tPCT_EXC_OVERLAP\tPCT_EXC_OFF_TARGET\tFOLD_80_BASE_PENALTY\tPCT_TARGET_BASES_1X\tPCT_TARGET_BASES_2X\tPCT_TARGET_BASES_10X\tPCT_TARGET_BASES_20X\tPCT_TARGET_BASES_30X\tPCT_TARGET_BASES_40X\tPCT_TARGET_BASES_50X\tPCT_TARGET_BASES_100X\tAT_DROPOUT\tGC_DROPOUT\tHET_SNP_SENSITIVITY\tHET_SNP_Q\tSAMPLE\tLIBRARY\tREAD_GROUP$\n"
            "targets\t685969\t1\t2328224315\t943367002\t2826386472\t0.536504\t0.463496\t0.711649\t3394.066372\t0.376355\t0.140761\t1746.272157\t9529653\t5.762479\t5.785637\t5.806222\t5.824234\t5.841216\t5.952375\t685969\t3137454505\t43799950\t43799950\t6186249337\t28824517\t28763941\t6097977789\t4013847317\t870783730\t1\t0.658095\t0.997898\t1269.421402\t1266\t3208\t0\t0.000384\t0.341774\t0\t0.052351\t0.011294\t0.169711\t0.282131\t1.245752\t0.999781\t0.999755\t0.999577\t0.999213\t0.999071\t0.998914\t0.998759\t0.998003\t6.637439\t0.31588\t0.166542\t1\t\t\t$\n"
            "$\n"
            "## HISTOGRAM\tjava.lang.Integer$\n"
            "coverage_or_base_quality\thigh_quality_coverage_count\tunfiltered_baseq_count$\n"
            "0\t150\t0$\n"
            "1\t18\t0$\n"
            "2\t17\t0$\n"
            "3\t9\t0$\n"
            "4\t10\t0$\n"
            "5\t5\t0$\n"
            "6\t13\t0$\n"
            "7\t17\t0$\n"
            "8\t25\t0$\n"
            "9\t26\t0$\n"
            "10\t47\t0$\n"
        )

    return target_output


def expected_hsmetrics_content_df():
    """Example required contents from hsmetrics file with actual metrics"""
    return pl.DataFrame(
        [
            {
                "BAIT_SET": "targets",
                "BAIT_TERRITORY": "685969",
                "BAIT_DESIGN_EFFICIENCY": "1",
                "ON_BAIT_BASES": "2328224315",
                "NEAR_BAIT_BASES": "943367002",
                "OFF_BAIT_BASES": "2826386472",
                "PCT_SELECTED_BASES": "0.536504",
                "PCT_OFF_BAIT": "0.463496",
                "ON_BAIT_VS_SELECTED": "0.711649",
                "MEAN_BAIT_COVERAGE": "3394.066372",
                "PCT_USABLE_BASES_ON_BAIT": "0.376355",
                "PCT_USABLE_BASES_ON_TARGET": "0.140761",
                "FOLD_ENRICHMENT": "1746.272157",
                "HS_LIBRARY_SIZE": "9529653",
                "HS_PENALTY_10X": "5.762479",
                "HS_PENALTY_20X": "5.785637",
                "HS_PENALTY_30X": "5.806222",
                "HS_PENALTY_40X": "5.824234",
                "HS_PENALTY_50X": "5.841216",
                "HS_PENALTY_100X": "5.952375",
                "TARGET_TERRITORY": "685969",
                "GENOME_SIZE": "3137454505",
                "TOTAL_READS": "43799950",
                "PF_READS": "43799950",
                "PF_BASES": "6186249337",
                "PF_UNIQUE_READS": "28824517",
                "PF_UQ_READS_ALIGNED": "28763941",
                "PF_BASES_ALIGNED": "6097977789",
                "PF_UQ_BASES_ALIGNED": "4013847317",
                "ON_TARGET_BASES": "870783730",
                "PCT_PF_READS": "1",
                "PCT_PF_UQ_READS": "0.658095",
                "PCT_PF_UQ_READS_ALIGNED": "0.997898",
                "MEAN_TARGET_COVERAGE": "1269.421402",
                "MEDIAN_TARGET_COVERAGE": "1266",
                "MAX_TARGET_COVERAGE": "3208",
                "MIN_TARGET_COVERAGE": "0",
                "ZERO_CVG_TARGETS_PCT": "0.000384",
                "PCT_EXC_DUPE": "0.341774",
                "PCT_EXC_ADAPTER": "0",
                "PCT_EXC_MAPQ": "0.052351",
                "PCT_EXC_BASEQ": "0.011294",
                "PCT_EXC_OVERLAP": "0.169711",
                "PCT_EXC_OFF_TARGET": "0.282131",
                "FOLD_80_BASE_PENALTY": "1.245752",
                "PCT_TARGET_BASES_1X": "0.999781",
                "PCT_TARGET_BASES_2X": "0.999755",
                "PCT_TARGET_BASES_10X": "0.999577",
                "PCT_TARGET_BASES_20X": "0.999213",
                "PCT_TARGET_BASES_30X": "0.999071",
                "PCT_TARGET_BASES_40X": "0.998914",
                "PCT_TARGET_BASES_50X": "0.998759",
                "PCT_TARGET_BASES_100X": "0.998003",
                "AT_DROPOUT": "6.637439",
                "GC_DROPOUT": "0.31588",
                "HET_SNP_SENSITIVITY": "0.166542",
                "HET_SNP_Q": "1",
                "SAMPLE": "",
                "LIBRARY": "",
                "READ_GROUP$": "$",
            }
        ]
    )
