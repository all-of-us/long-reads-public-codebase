########################################################################################
# SexCheck (upstream workflow name: SexCheckNaive)
# Phase 1 stage a, QC. Not on the critical path.
#
# Infers sex from one BAM's chrX and chrY coverage relative to chr1 (mosdepth) and
# compares it with the expected sex.
#
# Runs once per sample. This copy is a later version than the Phase 1 runs.
# Inputs:  bam, bai, expected_sex_type, mosdepth_summary_txt (optional)
# Outputs: inferred_sex_info
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/tasks/QC/AlignedMetrics.wdl" as AM
import "repo/wdl/tasks/QC/SexConcordance.wdl" as SC

workflow SexCheckNaive {
    input {
        File bam
        File bai
        String expected_sex_type

        File? mosdepth_summary_txt
    }

    if (!defined(mosdepth_summary_txt)) {
        call AM.MosDepthWGS {input: bam = bam, bai = bai}
    }

    call SC.SummarizeCoverages {input: mosdepth_summary_txt = select_first([mosdepth_summary_txt, MosDepthWGS.summary_txt])}
    call SC.MakeACall {
        input:
            cov_chr1 = SummarizeCoverages.cov_chr1,
            cov_chrX = SummarizeCoverages.cov_chrX,
            cov_chrY = SummarizeCoverages.cov_chrY,
            expected_sex_type = expected_sex_type
    }

    output {
        Map[String, String] inferred_sex_info = MakeACall.inferred_sex_info
    }
}
