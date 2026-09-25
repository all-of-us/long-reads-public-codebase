########################################################################################
# PAV2SVs
# Phase 1 stage a (per-sample processing).
#
# Converts a PAV VCF into an SV VCF (records of at least min_sv_length bp) and a
# separate SNV VCF, for use as the PAV input to TruvariIntrasampleV11.
#
# Runs once per sample.
# Inputs:  input_vcf_gz (PAV output), min_sv_length
# Outputs: sv_vcf_gz   SVs;  snp_vcf_gz  SNVs
# Previous: pav. Next: TruvariIntrasampleV11.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


#
workflow PAV2SVs {
    input {
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    call PAV2SVsImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length
    }
    
    output {
        File sv_vcf_gz = PAV2SVsImpl.sv_vcf_gz
        File sv_tbi = PAV2SVsImpl.sv_tbi
        File snp_vcf_gz = PAV2SVsImpl.snp_vcf_gz
        File snp_tbi = PAV2SVsImpl.snp_tbi
    }
}


task PAV2SVsImpl {
    input {
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} gunzip -c ~{input_vcf_gz} > input.vcf
        ${TIME_COMMAND} java -cp ~{docker_dir} PAV2SVs input.vcf ~{min_sv_length} sv.vcf snp.vcf
        ${TIME_COMMAND} bgzip sv.vcf
        tabix -f sv.vcf.gz
        ${TIME_COMMAND} bgzip snp.vcf
        tabix -f snp.vcf.gz
    >>>
    
    output {
        File sv_vcf_gz = work_dir + "/sv.vcf.gz"
        File sv_tbi = work_dir + "/sv.vcf.gz.tbi"
        File snp_vcf_gz = work_dir + "/snp.vcf.gz"
        File snp_tbi = work_dir + "/snp.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
