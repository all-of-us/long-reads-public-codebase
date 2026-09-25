########################################################################################
# PhysicalPhasing (the WDL names its workflow PhasicalPhasing)
# Phase 1 stage c (reference panel), launch step 6.
#
# Physically phases one sample's small variants and SVs together with HiPhase, using
# the sample's read alignments. Existing phasing in the input VCFs is removed first.
# Tandem-repeat genotypes (TRGT) can be phased at the same time.
#
# Runs once per sample. There is no example_inputs.json because no Phase 1 run was found.
# Inputs:  all_chr_bam, small_vcfs_directory and sv_vcfs_directory (locations of the
#          sample's SNV/indel and SV VCFs), reference_fasta, sample_id,
#          trgt_vcf (optional), hiphase_extra_args
# Outputs: hiphase_short_vcf, hiphase_sv_vcf  physically phased VCFs
# Previous: LRJointCallGVCFs and MergeRegenotypedIntersampleVcf.
# Next:     ChromosomePhasedPanelCreationFromHiPhase
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
########################################################################################

version 1.0

import "repo/wdl/methods/phasing/Helper.wdl" as H


workflow PhasicalPhasing {

    input {
        File all_chr_bam
        File all_chr_bai
        File reference_fasta
        File reference_fasta_fai
        String small_vcfs_directory
        String sv_vcfs_directory
        File? trgt_vcf
        File? trgt_vcf_tbi
        String prefix
        Int hiphase_memory
        String hiphase_extra_args
        String sample_id
    }

    call ConvertLowerCase {input:
        vcf = sv_vcfs_directory + "/" + sample_id + '.vcf.gz',
        prefix = sample_id + ".uppercased_sv_cleaned"
            
    }

    call UnphaseGenotypes as UnphaseSVGenotypes { input:
        vcf = ConvertLowerCase.subset_vcf,
        vcf_tbi = ConvertLowerCase.subset_tbi,
        prefix = prefix + ".unphased"
    }


    call H.HiphaseAll as HiphaseAll { input:
        bam = all_chr_bam,
        bai = all_chr_bai,
        unphased_snp_vcf = small_vcfs_directory + "/" + sample_id + '.vcf.gz',
        unphased_sv_vcf = UnphaseSVGenotypes.unphased_vcf,
        unphased_trgt_vcf = trgt_vcf,
        unphased_trgt_vcf_tbi = trgt_vcf_tbi,
        ref_fasta = reference_fasta,
        ref_fasta_fai = reference_fasta_fai,
        samplename = sample_id,
        memory = hiphase_memory,
        extra_args = hiphase_extra_args
    }

    output {
        File hiphase_short_vcf = HiphaseAll.phased_snp_vcf
        File hiphase_sv_vcf = HiphaseAll.phased_sv_vcf
        File? hiphase_trgt_vcf = HiphaseAll.phased_trgt_vcf
        File? haplotag_file = HiphaseAll.haplotag_file
    }
}

task ConvertLowerCase {
    input {
        File vcf
        String prefix
    }
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/mnt/disks/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.vcf
        bgzip ~{prefix}.vcf ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{prefix}.vcf.gz"
        File subset_tbi = "~{work_dir}/~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "8 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }
}

task UnphaseGenotypes {

    input {
        File vcf
        File vcf_tbi
        String prefix
    }

    command <<<
        set -euxo pipefail

        # set (a)ll genotypes to (u)nphased and sort by allele (e.g., 1|0 becomes 0/1)
        bcftools +setGT ~{vcf} -Oz -o ~{prefix}.vcf.gz -- --target-gt a --new-gt u
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File unphased_vcf = "~{prefix}.vcf.gz"
        File unphased_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
