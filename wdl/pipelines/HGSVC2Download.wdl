version 1.0


#
workflow HGSVC2Download {
    input {
        Array[String] sample_ids
        Array[File] bam_addresses
        Int target_coverage
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
        remote_dir: "Root directory in the remote bucket. Every sample is stored in a subdirectory."
    }
    
    scatter(i in range(length(sample_ids))) {
        call HGSVC2DownloadImpl {
            input:
                sample_id = sample_ids[i],
                bam_addresses = bam_addresses[i],
                target_coverage = target_coverage,
                remote_dir = remote_dir,
                n_cpus = n_cpus,
                ram_size_gb = ram_size_gb
        }
    }
    
    output {
    }
}


task HGSVC2DownloadImpl {
    input {
        String sample_id
        File bam_addresses
        Int target_coverage
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = (3*target_coverage)*2 + 120
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        TARGET_N_BYTES=$(( ~{target_coverage} * 3000000000 * 2 ))
        touch tmp1.fastq
        while read ADDRESS; do
            ${TIME_COMMAND} wget ${ADDRESS}
            FILE_NAME=$(basename ${ADDRESS})
            FILE_NAME=${FILE_NAME%.bam}
            ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n >> tmp1.fastq 
            N_BYTES=$(wc -c tmp1.fastq)
            if [ ${N_BYTES} -gt ${TARGET_N_BYTES} ]; then
                break
            fi
        done < ~{bam_addresses}
        head -c ${TARGET_N_BYTES} tmp1.fastq > tmp2.fastq
        rm -f tmp1.fastq
        N_ROWS=$(wc -l < tmp2.fastq)
        N_ROWS=$(( (${N_ROWS}/4)*4 ))
        head -n ${N_ROWS} tmp2.fastq | gzip > ~{sample_id}.fastq.gz
        rm -f tmp2.fastq
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.fastq.gz ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/hgsvc2"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
