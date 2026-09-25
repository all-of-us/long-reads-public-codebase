########################################################################################
# GetRegenotypedVcfKanpigMergedV11
# Phase 1 stage c (reference panel), step 8 of the callset-integration launch order.
#
# Genotypes one sample at every site of the cohort SV VCF with Kanpig, from the
# sample's read alignments and a sex-specific ploidy BED.
#
# Runs once per sample.
# Inputs:  intersample_vcf_gz (from RemoveSamplesV11), sample_id, alignments_bam,
#          reference_fa, ploidy_bed_female / ploidy_bed_male, is_male / sex
# Outputs: regenotyped_kanpig  the sample's genotypes at all cohort sites
# Previous: RemoveSamplesV11. Next: MergeRegenotypedIntersampleVcf.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


# Re-genotypes an inter-sample VCF.
#
workflow GetRegenotypedVcfKanpigMergedV11 {
    input {
        String sample_id
        Boolean is_male = false
        String sex = "F"
        File intersample_vcf_gz
        File intersample_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to have a single (artificial) sample column and all GTs equal to 0/1."
    }

    call GetRegenotypedVcfImpl {
        input:
            sample_id = sample_id,
            is_male = is_male,
            sex = sex,
            intersample_vcf_gz = intersample_vcf_gz,
            intersample_tbi = intersample_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male
    }
    
    output {
        File regenotyped_kanpig = GetRegenotypedVcfImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = GetRegenotypedVcfImpl.regenotyped_kanpig_tbi
    }
}


task GetRegenotypedVcfImpl {
    input {
        String sample_id
        Boolean is_male
        String sex
        File intersample_vcf_gz
        File intersample_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
        intersample_vcf_gz: "Assumed to have a single (artificial) sample column and all GTs equal to 0/1."
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "kanpig_regenotyped"
    String kanpig_params_multisample = "--sizemin 20 --sizemax 10000 --chunksize  500 --gpenalty 0.04 --hapsim 0.97"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 10*ceil(size(intersample_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))
    
    Int n_cpu = 8
    Int mem_gb = 16  # 2*n_cpu suggested by Adam
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        df -h
        
        # Making sure the VCF has the right sample ID
        echo ~{sample_id} > samples.txt
        ${TIME_COMMAND} bcftools reheader --threads ${N_THREADS} --samples samples.txt ~{intersample_vcf_gz} > tmp1.vcf.gz
        rm -f ~{intersample_vcf_gz}
        tabix -f tmp1.vcf.gz
        
        # Re-genotyping
        touch ~{alignments_bai}
        if [ ~{is_male} == "true" -o ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/kanpig --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_multisample} --reference ~{reference_fa} --input tmp1.vcf.gz --bam ~{alignments_bam} --out tmp2.vcf.gz
        rm -f tmp1.vcf.gz
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp2.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_kanpig = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: n_cpu
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
