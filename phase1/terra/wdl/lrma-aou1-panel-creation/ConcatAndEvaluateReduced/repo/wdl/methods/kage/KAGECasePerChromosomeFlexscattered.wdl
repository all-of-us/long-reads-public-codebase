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

workflow KAGECasePerChromosome {
    input {
        File input_cram
        File? input_crai
        String sample_name
        File reference_fasta
        File reference_fasta_fai
        File reference_dict

        # per chromosome quantities, grouped into scatter chunks
        Array[Array[String]+] chromosomes
        Array[Array[File]+] panel_index
        Array[Array[File]+] panel_kmer_index_only_variants_with_revcomp
        Array[Array[File]+] panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        Array[Array[File]+] panel_multi_split_vcf_gz_tbi

        Float average_coverage

        String samtools_docker
        String kage_docker
        File? monitoring_script

        RuntimeAttributes kage_runtime_attributes = {}
    }

    if (!defined(input_crai)) {
        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                docker = samtools_docker,
                monitoring_script = monitoring_script
        }
    }

    scatter (scatter_index in range(length(chromosomes))) {
        call KAGE {
            input:
                scatter_index = scatter_index,
                input_cram = input_cram,
                input_cram_idx = select_first([input_crai, IndexCaseReads.cram_idx]),
                panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp[scatter_index],
                panel_index = panel_index[scatter_index],
                panel_multi_split_vcf_gz = panel_multi_split_vcf_gz[scatter_index],
                panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi[scatter_index],
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                chromosomes = chromosomes[scatter_index],
                subset_reads = true,
                average_coverage = average_coverage,
                output_prefix = sample_name,
                sample_name = sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = kage_runtime_attributes
        }
    }

    call SerializeArrays {
        input:
            vcfs_array = flatten(KAGE.chromosome_kage_vcf_gzs),
            tbis_array = flatten(KAGE.chromosome_kage_vcf_gz_tbis),
            prefix = sample_name,
            docker = samtools_docker
    }

    output {
        File cram_idx = select_first([input_crai, IndexCaseReads.cram_idx])
        Array[File] chromosome_kmer_counts = flatten(KAGE.chromosome_kmer_counts)
        Array[File] chromosome_kage_vcf_gzs = flatten(KAGE.chromosome_kage_vcf_gzs)
        Array[File] chromosome_kage_vcf_gz_tbis = flatten(KAGE.chromosome_kage_vcf_gz_tbis)
        # serialized array outputs
        File chromosome_kage_vcf_gzs_txt = SerializeArrays.vcfs_txt
        File chromosome_kage_vcf_gz_tbis_txt = SerializeArrays.tbis_txt
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
        set -euox pipefail

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

task KAGE {
    input {
        Int scatter_index
        File input_cram
        File input_cram_idx
        Array[File] panel_kmer_index_only_variants_with_revcomp
        Array[File] panel_index
        Array[File] panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        Array[File] panel_multi_split_vcf_gz_tbi
        String output_prefix
        String sample_name
        Float average_coverage
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        Array[String]+ chromosomes
        Boolean subset_reads

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 10000000"

        Boolean? ignore_helper_model = true
        String? kage_genotype_extra_args

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    # explicitly specify CPU and memory separately for KAGE commands;
    # using nproc w/ automatic CPU/memory scaling of custom instances on Terra can be suboptimal or cause OOM
    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    Int num_chromosomes = length(chromosomes)

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        CHROMOSOMES=(~{sep=" " chromosomes})
        KMER_INDICES=(~{sep=" " panel_kmer_index_only_variants_with_revcomp})
        INDICES=(~{sep=" " panel_index})
        MULTI_SPLITS=(~{sep=" " panel_multi_split_vcf_gz})

        NUM_CHROMOSOMES=~{num_chromosomes}
        NUM_DIGITS=${#NUM_CHROMOSOMES}

        for (( c = 0; c < ~{num_chromosomes}; c++ ))
        do
            CHROMOSOME=${CHROMOSOMES[$c]}
            KMER_INDEX=${KMER_INDICES[$c]}
            INDEX=${INDICES[$c]}
            MULTI_SPLIT=${MULTI_SPLITS[$c]}

            C_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $c)

            date
            mkdir outputs-$C_WITH_LEADING_ZEROS
            mkfifo outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.preprocessed.fa
            if ~{subset_reads}; then
                # hacky way to get chromosomes into bed file
                grep -P "$CHROMOSOME\t" ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > outputs-$C_WITH_LEADING_ZEROS/chromosome.bed

                echo "Subsetting reads..."
                samtools view --reference ~{reference_fasta} -@ $(nproc) ~{if subset_reads then "--regions-file outputs-$C_WITH_LEADING_ZEROS/chromosome.bed" else ""} -u -X ~{input_cram} ~{input_cram_idx} | \
                    samtools fasta --reference ~{reference_fasta} -@ $(nproc) > outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.preprocessed.fa &
            else
                echo "Not subsetting reads..."
                samtools fasta --reference ~{reference_fasta} -@ $(nproc) ~{input_cram} > outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.preprocessed.fa &
            fi

            kmer_mapper map \
                ~{kmer_mapper_args} \
                -t ~{cpu_resolved} \
                -i $KMER_INDEX \
                -f outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.preprocessed.fa \
                -o ~{output_prefix}.shard-~{scatter_index}-output-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kmer_counts.npy

            kage genotype \
                -i $INDEX \
                -c ~{output_prefix}.shard-~{scatter_index}-output-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kmer_counts.npy \
                --average-coverage ~{average_coverage} \
                -s ~{sample_name} \
                ~{true='-I true' false='-I false' ignore_helper_model} \
                ~{kage_genotype_extra_args} \
                -o outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.kage.bi.vcf
            bcftools view outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.kage.bi.vcf -Ob -o outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.kage.bi.bcf
            bcftools index outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.kage.bi.bcf

            # we need to add split multiallelics to biallelic-only KAGE BCF
            # create single-sample header from LOO panel w/ split multiallelics
            bcftools view --no-version -h -G $MULTI_SPLIT | \
                sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
                sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.header.txt
            # create single-sample missing genotypes from LOO panel w/ split multiallelics
            bcftools view --no-version -H -G $MULTI_SPLIT | \
                sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.GT.txt
            # create single-sample BCF w/ split multiallelics
            bcftools view <(cat outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.header.txt outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.GT.txt) -Ob -o outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.bcf
            bcftools index outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.bcf

            bcftools concat --no-version -a outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.kage.bi.bcf outputs-$C_WITH_LEADING_ZEROS/~{output_prefix}.$CHROMOSOME.multi.split.bcf -Oz -o ~{output_prefix}.shard-~{scatter_index}-output-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.vcf.gz
            bcftools index -t ~{output_prefix}.shard-~{scatter_index}-output-$C_WITH_LEADING_ZEROS.$CHROMOSOME.kage.vcf.gz

            rm -r outputs-$C_WITH_LEADING_ZEROS
        done
    >>>

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 115]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        Array[File] chromosome_kmer_counts = glob("~{output_prefix}.*.kmer_counts.npy")       # must be careful with globbing and lexicographical ordering here!
        Array[File] chromosome_kage_vcf_gzs = glob("~{output_prefix}.*.kage.vcf.gz")
        Array[File] chromosome_kage_vcf_gz_tbis = glob("~{output_prefix}.*.kage.vcf.gz.tbi")
    }
}

task SerializeArrays {
    input {
        Array[String] vcfs_array
        Array[String] tbis_array
        String prefix

        String docker

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        mv ~{write_lines(vcfs_array)} ~{prefix}.chromosome_kage_vcf_gzs.txt
        mv ~{write_lines(tbis_array)} ~{prefix}.chromosome_kage_vcf_gz_tbis.txt
    }

    output {
        File vcfs_txt = "~{prefix}.chromosome_kage_vcf_gzs.txt"
        File tbis_txt = "~{prefix}.chromosome_kage_vcf_gz_tbis.txt"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 3]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}
