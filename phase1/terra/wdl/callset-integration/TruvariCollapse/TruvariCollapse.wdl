########################################################################################
# TruvariCollapse
# Helper, not on the Phase 1 launch path.
#
# Runs truvari collapse on one VCF, merging SVs that describe the same event.
#
# Inputs:  intput_vcf_gz (sic) and its index
# Outputs: truvari_collapsed  the collapsed VCF
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


#
workflow TruvariCollapse {
    input {
        File intput_vcf_gz
        File intput_vcf_gz_tbi
    }
    parameter_meta {
    }

    call TruvariCollapseImpl {
        input:
            intput_vcf_gz = intput_vcf_gz,
            intput_vcf_gz_tbi = intput_vcf_gz_tbi
    }
    
    output {
        File truvari_collapsed = TruvariCollapseImpl.truvari_collapsed
    }
}


task TruvariCollapseImpl {
    input {
        File intput_vcf_gz
        File intput_vcf_gz_tbi
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(intput_vcf_gz,"GB"))) + 50
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} truvari collapse -i ~{intput_vcf_gz} > truvari_collapsed.vcf
    >>>
    
    output {
        File truvari_collapsed = "~{work_dir}/truvari_collapsed.vcf"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
