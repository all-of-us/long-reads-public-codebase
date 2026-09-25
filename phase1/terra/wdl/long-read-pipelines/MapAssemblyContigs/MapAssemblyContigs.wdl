########################################################################################
# MapAssemblyContigs
# Phase 1 stage a (per-sample processing), optional.
#
# Aligns assembly contigs to the reference with minimap2.
#
# Runs once per assembly (haplotype).
# Inputs:  contigs_fa, sample, read_group, ref_map_file (table of reference files)
# Outputs: aligned_bam  contigs aligned to the reference
# Previous: PBAssembleWithHifiasm.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/tasks/AlignReads.wdl"

workflow MapAssemblyContigs {
    input {
        File contigs_fa
        String read_group
        String sample
        File ref_map_file
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call AlignReads.Minimap2 { input:
        reads = [contigs_fa],
        ref_fasta = ref_map['fasta'],
        RG = "@RG\\tID:~{read_group}\\tSM:~{sample}",
        map_preset = "asm5"
    }
    output {
        File aligned_bam = Minimap2.aligned_bam
        File aligned_bai = Minimap2.aligned_bai
    }
}