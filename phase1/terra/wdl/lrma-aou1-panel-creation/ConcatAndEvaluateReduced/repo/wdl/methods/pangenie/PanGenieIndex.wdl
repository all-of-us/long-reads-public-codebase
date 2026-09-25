version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow PanGenieIndex {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String pangenie_docker
        File? monitoring_script
        RuntimeAttributes? pangenie_runtime_attributes
    }

    call Index {
        input:
            panel_vcf_gz = panel_vcf_gz,
            panel_vcf_gz_tbi = panel_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            docker = pangenie_docker,
            monitoring_script = monitoring_script,
            runtime_attributes = pangenie_runtime_attributes
    }

    output {
        Array[File] pangenie_index_chromosome_graphs = Index.pangenie_index_chromosome_graphs
        Array[File] pangenie_index_chromosome_kmers = Index.pangenie_index_chromosome_kmers
        File pangenie_index_unique_kmers_map = Index.pangenie_index_unique_kmers_map
        File pangenie_index_path_segments_fasta = Index.pangenie_index_path_segments_fasta
    }
}

task Index {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File reference_fasta
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31
        String? extra_args

        RuntimeAttributes runtime_attributes = {}
    }

    Int num_chromosomes = length(chromosomes)

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # subset reference and panel VCF to chromosomes
        samtools faidx -r <(echo -e "~{sep="\n" chromosomes}") ~{reference_fasta} > reference.subset.fa
        bcftools view ~{panel_vcf_gz} -r ~{sep="," chromosomes} > panel.subset.vcf

        NPROC=$(nproc)
        # NUM_CHROMOSOMES=~{num_chromosomes}
        # NUM_THREADS=$(( NPROC < NUM_CHROMOSOMES ? NPROC : NUM_CHROMOSOMES ))
        NUM_THREADS=$NPROC

        /pangenie/build/src/PanGenie-index \
            -o ~{output_prefix} \
            -r reference.subset.fa \
            -t $NUM_THREADS \
            -v panel.subset.vcf \
            -k ~{kmer_length} \
            ~{extra_args}
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        Array[File] pangenie_index_chromosome_graphs = glob("~{output_prefix}_*_Graph.cereal")
        Array[File] pangenie_index_chromosome_kmers = glob("~{output_prefix}_*_kmers.tsv.gz")
        File pangenie_index_unique_kmers_map = "~{output_prefix}_UniqueKmersMap.cereal"
        File pangenie_index_path_segments_fasta = "~{output_prefix}_path_segments.fasta"
    }
}
