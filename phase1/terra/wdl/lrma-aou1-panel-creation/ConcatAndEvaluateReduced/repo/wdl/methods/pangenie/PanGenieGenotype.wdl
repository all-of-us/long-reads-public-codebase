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

workflow PanGenieGenotype {
    input {
        Array[File] pangenie_index_chromosome_graphs
        Array[File] pangenie_index_chromosome_kmers
        File pangenie_index_unique_kmers_map
        File pangenie_index_path_segments_fasta
        String index_prefix
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        File input_cram
        File? input_crai
        Boolean subset_reads = true
        String sample_name

        String? pangenie_extra_args

        String docker
        String pangenie_docker
        File? monitoring_script
        RuntimeAttributes? pangenie_runtime_attributes
    }

    if (!defined(input_crai)) {
        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                docker = docker,
                monitoring_script = monitoring_script
        }
    }

    if (subset_reads) {
        call PreprocessCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                input_cram_idx = select_first([input_crai, IndexCaseReads.cram_idx]),
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = sample_name,
                chromosomes = chromosomes,
                docker = docker,
                monitoring_script = monitoring_script
        }
    }

    if (!subset_reads) {
        call PreprocessCaseReadsWithoutSubsetting {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                input_cram_idx = select_first([input_crai, IndexCaseReads.cram_idx]),
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = sample_name,
                docker = docker,
                monitoring_script = monitoring_script
        }
    }

    call Genotype {
        input:
            pangenie_index_chromosome_graphs = pangenie_index_chromosome_graphs,
            pangenie_index_chromosome_kmers = pangenie_index_chromosome_kmers,
            pangenie_index_unique_kmers_map = pangenie_index_unique_kmers_map,
            pangenie_index_path_segments_fasta = pangenie_index_path_segments_fasta,
            index_prefix = index_prefix,
            input_fasta = select_first([PreprocessCaseReads.preprocessed_fasta, PreprocessCaseReadsWithoutSubsetting.preprocessed_fasta]),
            sample_name = sample_name,
            output_prefix = sample_name,
            docker = pangenie_docker,
            extra_args = pangenie_extra_args,
            monitoring_script = monitoring_script,
            runtime_attributes = pangenie_runtime_attributes
    }

    output {
        File? cram_idx = IndexCaseReads.cram_idx
        File genotyping_vcf_gz = Genotype.genotyping_vcf_gz
        File genotyping_vcf_gz_tbi = Genotype.genotyping_vcf_gz_tbi
        File genotyping_naively_phased_vcf_gz = Genotype.genotyping_naively_phased_vcf_gz
        File genotyping_naively_phased_vcf_gz_tbi = Genotype.genotyping_naively_phased_vcf_gz_tbi
        File histogram = Genotype.histogram
#        File path_segments_fasta = Genotype.path_segments_fasta
    }
}

# some 1000G CRAM indices have issues due to htsjdk version, so we reindex; see e.g. https://github.com/broadinstitute/gatk/issues/7076
task IndexCaseReads {
    input {
        File input_cram

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    # if we instead simply use File cram_idx = "~{input_cram}.crai" in the output block,
    # Terra tries to localize an index adjacent to the CRAM in PreprocessCaseReads???
    String output_prefix = basename(input_cram, ".cram")

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        samtools index -@ $(nproc) ~{input_cram} ~{output_prefix}.cram.crai
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task PreprocessCaseReads {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # hacky way to get chromosomes into bed file
        grep -P '~{sep="\\t|" chromosomes}\t' ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosomes.bed

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools view --reference ~{reference_fasta} -@ $(nproc) -L chromosomes.bed -u ~{input_cram} | \
            samtools fasta --reference ~{reference_fasta} -@ $(nproc) | \
            sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task PreprocessCaseReadsWithoutSubsetting {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools fasta --reference ~{reference_fasta} -@ $(nproc) ~{input_cram} | sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task Genotype {
    input {
        Array[File] pangenie_index_chromosome_graphs
        Array[File] pangenie_index_chromosome_kmers
        File pangenie_index_unique_kmers_map
        File pangenie_index_path_segments_fasta
        String index_prefix
        File input_fasta
        String sample_name
        String output_prefix

        String docker
        File? monitoring_script
        String? extra_args

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        NPROC=$(nproc)
        NUM_THREADS=$NPROC

        mv ~{sep=" " pangenie_index_chromosome_graphs} .
        mv ~{sep=" " pangenie_index_chromosome_kmers} .
        mv ~{pangenie_index_unique_kmers_map} .
        mv ~{pangenie_index_path_segments_fasta} .

        /pangenie/build/src/PanGenie \
            -i ~{input_fasta} \
            -f ~{index_prefix} \
            -o ~{output_prefix} \
            -s ~{sample_name} \
            -t $NUM_THREADS \
            -j $NUM_THREADS \
            ~{extra_args}

        bgzip -c ~{output_prefix}_genotyping.vcf > ~{output_prefix}_genotyping.vcf.gz
        bcftools index -t ~{output_prefix}_genotyping.vcf.gz

        # naively set all GTs to phased for Vcfdist evaluation
        bcftools +setGT ~{output_prefix}_genotyping.vcf.gz -Oz -o ~{output_prefix}_genotyping_naively_phased.vcf.gz -- -t a -n p
        bcftools index -t ~{output_prefix}_genotyping_naively_phased.vcf.gz
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
        File genotyping_vcf_gz = "~{output_prefix}_genotyping.vcf.gz"
        File genotyping_vcf_gz_tbi = "~{output_prefix}_genotyping.vcf.gz.tbi"
        File genotyping_naively_phased_vcf_gz = "~{output_prefix}_genotyping_naively_phased.vcf.gz"
        File genotyping_naively_phased_vcf_gz_tbi = "~{output_prefix}_genotyping_naively_phased.vcf.gz.tbi"
        File histogram = "~{output_prefix}_histogram.histo"
#        File path_segments_fasta = "~{output_prefix}_path_segments.fasta"
    }
}
