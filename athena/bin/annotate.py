"""
Functions to annotate a panel bed file with transcript information
and per base coverage data using bedtools.

Requires: bedtools
"""
from pathlib import Path

import pybedtools as bedtools


class annotateBed():

    def add_transcript_info(self, panel_bed, transcript_info):
        """
        Use pybedtools to annotate panel bed file with coverage data

        Parameters
        ----------
        panel_bed : str
            panel bed file to intersect
        transcript_info : str
            transcript info file to intersect

        Returns
        -------
        bed_w_transcripts : pd.DataFrame
            dataframe of intersected panel bed file with
            transcript information
        """
        print("\nCalling bedtools to add transcript info")

        bed = bedtools.BedTool(panel_bed)
        all_transcripts = bedtools.BedTool(transcript_info)

        # intersecting panel bed file with transcript/gene/exon information
        # requires 100% overlap on panel -> transcript coordinates
        bed_w_transcript = bed.intersect(
            all_transcripts, wa=True, wb=True, F=1.0, sorted=True
        )

        # convert pybedtools object to dataframe
        bed_w_transcript = bed_w_transcript.to_dataframe(names=[
            "p_chrom", "p_start", "p_end", "p_transcript",
            "t_chrom", "t_start", "t_end", "t_gene", "t_transcript", "t_exon"
        ])

        # get list of unique transcripts from panel bed to ensure output
        # contains expected transcripts
        panel_transcripts = bed.to_dataframe(
            names=['chrom', 'start', 'end', 'transcript']
        ).transcript.unique().tolist()


        # check for empty file
        assert len(bed_w_transcript.index) > 0, """Empty file returned from
            intersecting panel bed and transcript file. Check if flanks are
            being used as 100% coordinate overlap is currently required."""

        # panel bed file defines transcript to use, filter transcript file for
        # just those transcripts
        bed_w_transcript = bed_w_transcript[
            bed_w_transcript["p_transcript"] == bed_w_transcript["t_transcript"]
        ]

        intersect_transcripts = bed_w_transcript.t_transcript.unique().tolist()

        # # ensure no transcripts dropped from panel due to missing from
        # # transcripts file
        assert len(panel_transcripts) == len(intersect_transcripts), (
            f"Transcript(s) dropped from panel during intersecting with "
            f"transcript file. Total before {len(panel_transcripts)}. Total "
            f"after {len(intersect_transcripts)}. Dropped transcripts: "
            f"{set(panel_transcripts) - set(intersect_transcripts)}"
        )

        # drop duplicate columns
        bed_w_transcript = bed_w_transcript.drop(columns=[
            't_chrom', 't_start', 't_end', 'p_transcript'
        ])

        return bed_w_transcript


    def add_coverage(self, bed_w_transcript, per_base_coverage, build):
        """
        Use pybedtools to add coverage bin data to selected panel regions

        Parameters
        ----------
        bed_w_transcript : pd.DataFrame
            panel bed file with transcript information
        coverage_dfs : str
            filename of per base coverage data
        build : int
            reference build to use (37 | 38), this determines which
            genome reference file to use when calling bedtools intersect
            from `athena/data/genomes`

        Returns
        -------
        bed_w_coverage : pd.DataFrame
            panel bed with transcript and coverage info
        """
        print("\nCalling bedtools to add coverage info")

        if build == 37:
            genome_file = "human.hg19.genome"
        elif build == 38:
            genome_file = "human.hg38.genome"
        else:
            raise RuntimeError('Only build 37 or 38 currently supported')

        genome_file = Path(__file__).parent.resolve().joinpath(
            f'../data/genomes/{genome_file}')

        # set bedtools objects
        bed_w_transcript = bedtools.BedTool.from_dataframe(bed_w_transcript)
        per_base_coverage = bedtools.BedTool(per_base_coverage)

        # intersect per base coverage data onto panel bed
        # sorted - required to invoke sweeping algorithm to reduce
        #   memory usage from large files
        # nonamecheck - set to allow differences between chromosome
        #   naming (i.e. chr1 vs chr01 vs 1)
        # g - genome file to use, required when using sorted so bedtools
        #   knows ranges of regions relative to chromosome lengths
        bed_w_coverage = bed_w_transcript.intersect(
            per_base_coverage,
            wa=True,
            wb=True,
            sorted=True,
            nonamecheck=True,
            g=genome_file
        )

        bed_w_coverage = bed_w_coverage.to_dataframe(names=[
            "chrom", "exon_start", "exon_end", "gene", "transcript",
            "exon", "cov_chrom", "cov_start", "cov_end", "cov"
        ])

        # check again for empty output of bedtools, can happen due to memory
        # maxing out and doesn't seem to raise an exception...
        assert len(bed_w_coverage.index) > 0, (
            "Error calling bedtools to intersect per base coverage data with "
            "panel bed, this is likely due to running out of available memory"
        )

        # drop duplicate chromosome col and rename
        bed_w_coverage.drop(columns=["cov_chrom"], inplace=True)

        return bed_w_coverage
