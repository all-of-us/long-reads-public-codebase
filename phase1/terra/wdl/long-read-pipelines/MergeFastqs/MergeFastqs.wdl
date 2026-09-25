########################################################################################
# MergeFastqs
# Phase 1 stage a (per-sample processing), optional preprocessing.
#
# Concatenates a sample's gzipped FASTQ files (for example, from several SMRT cells)
# into one file.
#
# Runs once per sample.
# Inputs:  fastq_gzs, sample, gcs_out_root_dir
# Outputs: FASTQ  the merged FASTQ
# Previous: HiFiBamToFastQ. Next: PBAssembleWithHifiasm.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/tasks/Utils.wdl"
import "repo/wdl/tasks/Finalize.wdl" as FF

workflow MergeFastqs {
    input {
        Array[File] fastq_gzs
        String sample

        String gcs_out_root_dir
    }

    if (1<length(fastq_gzs)) {
        call Utils.MergeFastqGzs {input:
            fastq_gzs = fastq_gzs,
            prefix = sample
        }
    }
    File to_save = select_first([MergeFastqGzs.merged_fastq_gz, fastq_gzs[0]])
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/FASTQs/"
    call FF.FinalizeToFile {input:
        file = to_save, outdir = outdir, name = "~{sample}.fastq.gz"
    }

    output {
        File FASTQ = FinalizeToFile.gcs_path
    }
}
