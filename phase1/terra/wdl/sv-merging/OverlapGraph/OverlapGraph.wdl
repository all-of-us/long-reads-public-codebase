########################################################################################
# OverlapGraph
# Helper from the SV-merging comparisons, not on the Phase 1 launch path.
#
# Builds a graph of overlapping SVs from a merged callset (for one chromosome and SV
# type) and reports its degree and connected-component statistics, to compare merging
# and re-genotyping methods.
#
# Inputs:  bucket_dir, only_chr, sv_type, merger, genotyper, n_samples, limits
# Outputs: report
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


#
workflow OverlapGraph {
    input {
        String bucket_dir
        String only_chr
        Int sv_type
        String merger
        String genotyper
        Int n_samples
        Int max_degree
        Int max_component_size
    }
    parameter_meta {
    }
    call OverlapGraphImpl {
        input:
            bucket_dir = bucket_dir,
            only_chr = only_chr,
            sv_type = sv_type,
            merger = merger,
            genotyper = genotyper,
            n_samples = n_samples,
            max_degree = max_degree,
            max_component_size = max_component_size
    }
    output {
        File report = OverlapGraphImpl.report
    }
}


task OverlapGraphImpl {
    input {
        String bucket_dir
        String only_chr
        Int sv_type
        String merger
        String genotyper
        Int n_samples
        Int max_degree
        Int max_component_size
    }
    parameter_meta {
        genotyper: "null: no regenotyping."
    }
    
    String docker_dir = "/sv-merging"
    String work_dir = "/cromwell_root/sv-merging"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if [ ~{merger} == "joint" ]; then
            ID="joint900_pass_chr2122_standardized"
        elif [ ~{genotyper} != "null" ]; then
            ID="regenotyped_~{genotyper}_standardized"
        else
            ID="standardized"
        fi
        VCF_GZ_FILE="${ID}.vcf.gz"
        gsutil cp ~{bucket_dir}/~{merger}/${VCF_GZ_FILE} .
        gunzip ${VCF_GZ_FILE}
        if [ ~{n_samples} -ne -1 ]; then
            seq 0 $(( ~{n_samples}-1 )) > list.txt
            N_ROWS=$(( ~{n_samples}/${N_THREADS} ))
            split -d -l ${N_ROWS} list.txt chunk-
            for SAMPLES_FILE in $(ls chunk-*); do
                java -Xmx4g -cp ~{docker_dir} OverlapGraph ${ID}.vcf ~{only_chr} ~{sv_type} ${SAMPLES_FILE} 1 1000 10000 . &
            done
            wait
        else
            java -Xmx2g -cp ~{docker_dir} OverlapGraph ${ID}.vcf ~{only_chr} ~{sv_type} null 1 1000 10000 .
        fi
        if [ ~{genotyper} != "null" ]; then
            INFIX="_~{genotyper}"
        else
            INFIX=""
        fi
        cat *_chr22_del_degree_histogram.txt > ~{merger}${INFIX}_chr22_del_degree_histogram.matrix
        cat *_chr22_del_componentSize_histogram.txt > ~{merger}${INFIX}_chr22_del_componentSize_histogram.matrix
        rm -f *.txt ${ID}.vcf
        tar -czvf report.tar.gz *.matrix
    >>>
    output {
        File report = work_dir + "/report.tar.gz"
    }
    runtime {
        docker: "fcunial/sv-merging"
        cpu: 16  # Arbitrary
        memory: "32GB"  # Arbitrary
        disks: "local-disk 250 HDD"  # Arbitrary
        preemptible: 0
    }
}
