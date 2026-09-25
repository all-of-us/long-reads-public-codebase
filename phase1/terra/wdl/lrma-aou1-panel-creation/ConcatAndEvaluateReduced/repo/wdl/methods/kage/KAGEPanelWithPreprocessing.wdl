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

workflow KAGEPanelWithPreprocessing {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        Boolean do_preprocessing = true # if false, input_vcf_gz is assumed to contain only biallelics
        File reference_fasta
        File reference_fasta_fai
        String output_prefix
        Array[String] chromosomes
        Boolean do_single_chromosome = false

        String docker
        File? monitoring_script

        Int? spacing

        RuntimeAttributes? runtime_attributes
        RuntimeAttributes? medium_runtime_attributes
        RuntimeAttributes? large_runtime_attributes
        Int? cpu_make_count_model
    }

    String output_prefix_ = if do_single_chromosome then output_prefix + "." + chromosomes[0] else output_prefix

    if (do_preprocessing) {
        scatter (i in range(length(chromosomes))) {
            call PreprocessPanelVCF {
                input:
                    input_vcf_gz = input_vcf_gz,
                    input_vcf_gz_tbi = input_vcf_gz_tbi,
                    chromosomes = [chromosomes[i]],
                    output_prefix = output_prefix_ + "." + chromosomes[i],
                    docker = docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = runtime_attributes
            }
        }

        call ConcatVcfs as ConcatVcfsPreprocessedPanel {
            input:
                vcf_gzs = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                vcf_gz_tbis = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                output_prefix = output_prefix_ + ".preprocessed",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call ConcatVcfs as ConcatVcfsPreprocessedPanelSplit {
            input:
                vcf_gzs = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                vcf_gz_tbis = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                output_prefix = output_prefix_ + ".preprocessed.split",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call ConcatVcfs as ConcatVcfsPreprocessedPanelBi {
            input:
                vcf_gzs = PreprocessPanelVCF.preprocessed_panel_bi_vcf_gz,
                vcf_gz_tbis = PreprocessPanelVCF.preprocessed_panel_bi_vcf_gz_tbi,
                output_prefix = output_prefix_ + ".preprocessed.bi",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call ConcatVcfs as ConcatVcfsPreprocessedPanelMultiSplit {
            input:
                vcf_gzs = PreprocessPanelVCF.preprocessed_panel_multi_split_vcf_gz,
                vcf_gz_tbis = PreprocessPanelVCF.preprocessed_panel_multi_split_vcf_gz_tbi,
                output_prefix = output_prefix_ + ".preprocessed.multi.split",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }
    }

    File preprocessed_panel_bi_vcf_gz_ = select_first([ConcatVcfsPreprocessedPanelBi.vcf_gz, input_vcf_gz])
    File preprocessed_panel_bi_vcf_gz_tbi_ = select_first([ConcatVcfsPreprocessedPanelBi.vcf_gz_tbi, input_vcf_gz_tbi])

    call MakeSitesOnlyVcfAndNumpyVariants {
        input:
            input_vcf_gz = preprocessed_panel_bi_vcf_gz_,
            input_vcf_gz_tbi = preprocessed_panel_bi_vcf_gz_tbi_,
            chromosomes = chromosomes,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    scatter (i in range(length(chromosomes))) {
        call MakeChromosomeGenotypeMatrix {
            input:
                input_vcf_gz = preprocessed_panel_bi_vcf_gz_,
                input_vcf_gz_tbi = preprocessed_panel_bi_vcf_gz_tbi_,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = medium_runtime_attributes
        }

        call MakeChromosomeGraph {
            input:
                input_vcf_gz = preprocessed_panel_bi_vcf_gz_,
                input_vcf_gz_tbi = preprocessed_panel_bi_vcf_gz_tbi_,
                reference_fasta = reference_fasta,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = medium_runtime_attributes
        }

        call MakeChromosomeVariantToNodes {
            input:
                chromosome_graph = MakeChromosomeGraph.chromosome_graph,
                sites_only_vcf_gz = MakeSitesOnlyVcfAndNumpyVariants.sites_only_vcf_gz,
                sites_only_vcf_gz_tbi = MakeSitesOnlyVcfAndNumpyVariants.sites_only_vcf_gz_tbi,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call MakeChromosomeHaplotypeToNodes {
            input:
                chromosome_variant_to_nodes = MakeChromosomeVariantToNodes.chromosome_variant_to_nodes,
                chromosome_genotype_matrix = MakeChromosomeGenotypeMatrix.chromosome_genotype_matrix,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = medium_runtime_attributes
        }
    }

    call MergeChromosomeGraphs {
        input:
            chromosome_graphs = MakeChromosomeGraph.chromosome_graph,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = large_runtime_attributes
    }

    call MergeChromosomeVariantToNodes {
        input:
            chromosome_variant_to_nodes = MakeChromosomeVariantToNodes.chromosome_variant_to_nodes,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    call MergeChromosomeHaplotypeToNodes {
        input:
            chromosome_haplotype_to_nodes = MakeChromosomeHaplotypeToNodes.chromosome_haplotype_to_nodes,
            chromosome_haplotype_nodes = MakeChromosomeHaplotypeToNodes.chromosome_haplotype_nodes,
            num_nodes = MergeChromosomeVariantToNodes.num_nodes,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    call MakeHelperModel {
        input:
        chromosome_genotype_matrices = MakeChromosomeGenotypeMatrix.chromosome_genotype_matrix,
        variant_to_nodes = MergeChromosomeVariantToNodes.variant_to_nodes,
        output_prefix = output_prefix_,
        docker = docker,
        monitoring_script = monitoring_script,
        runtime_attributes = medium_runtime_attributes
    }

    scatter (i in range(length(chromosomes))) {
        call SampleChromosomeKmersFromLinearReference {
            input:
                reference_fasta = reference_fasta,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                spacing = spacing,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = large_runtime_attributes
        }
    }

    call MergeFlatKmers as MergeChromosomeKmersFromLinearReference {
        input:
            flat_kmers = SampleChromosomeKmersFromLinearReference.chromosome_linear_kmers,
            num_nodes = MergeChromosomeVariantToNodes.num_nodes,
            reference_fasta_fai = reference_fasta_fai,
            output_prefix = "~{output_prefix_}.linear_kmers",
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = large_runtime_attributes
    }

    call MakeLinearReferenceKmerCounter {
        input:
            linear_kmers = MergeChromosomeKmersFromLinearReference.merged_kmers,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = large_runtime_attributes
    }

    scatter (i in range(length(chromosomes))) {
        call GetChromosomeShortVariantKmers {
            input:
                chromosome_graph = MakeChromosomeGraph.chromosome_graph[i],
                chromosome_position_id_index = MakeChromosomeGraph.chromosome_position_id_index[i],
                chromosome_variant_to_nodes = MakeChromosomeVariantToNodes.chromosome_variant_to_nodes[i],
                linear_kmers_counter = MakeLinearReferenceKmerCounter.linear_kmers_counter,
                sites_only_vcf_gz = MakeSitesOnlyVcfAndNumpyVariants.sites_only_vcf_gz,
                sites_only_vcf_gz_tbi = MakeSitesOnlyVcfAndNumpyVariants.sites_only_vcf_gz_tbi,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = large_runtime_attributes
        }

        call SampleChromosomeStructuralVariantKmers {
            input:
                chromosome_graph = MakeChromosomeGraph.chromosome_graph[i],
                chromosome_variant_to_nodes = MakeChromosomeVariantToNodes.chromosome_variant_to_nodes[i],
                linear_kmers_counter = MakeLinearReferenceKmerCounter.linear_kmers_counter,
                chromosome = chromosomes[i],
                output_prefix = output_prefix_,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = large_runtime_attributes
        }

        call MergeFlatKmers as MergeChromosomeShortAndStructuralVariantKmers {
            input:
                flat_kmers = [GetChromosomeShortVariantKmers.chromosome_short_variant_kmers, SampleChromosomeStructuralVariantKmers.chromosome_structural_variant_kmers],
                output_prefix = "~{output_prefix_}.~{chromosomes[i]}.variant_kmers",
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = large_runtime_attributes
        }
    }

    call MergeFlatKmers as MergeChromosomeVariantKmers {
        input:
            flat_kmers = MergeChromosomeShortAndStructuralVariantKmers.merged_kmers,
            num_nodes = MergeChromosomeVariantToNodes.num_nodes,
            reference_fasta_fai = reference_fasta_fai,
            output_prefix = "~{output_prefix_}.variant_kmers",
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = large_runtime_attributes
    }

    call MakeReverseVariantKmerIndex {
        input:
            variant_kmers = MergeChromosomeVariantKmers.merged_kmers,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    call MakeVariantKmerIndexWithReverseComplements {
        input:
            variant_kmers = MergeChromosomeVariantKmers.merged_kmers,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    call MakeCountModel {
        input:
            graph = MergeChromosomeGraphs.graph,
            haplotype_to_nodes = MergeChromosomeHaplotypeToNodes.haplotype_to_nodes,
            haplotype_nodes = MergeChromosomeHaplotypeToNodes.haplotype_nodes,
            kmer_index_only_variants_with_revcomp = MakeVariantKmerIndexWithReverseComplements.kmer_index_only_variants_with_revcomp,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = large_runtime_attributes,
            cpu = cpu_make_count_model
    }

    call RefineCountModel {
        input:
            variant_to_nodes = MergeChromosomeVariantToNodes.variant_to_nodes,
            sampling_count_model = MakeCountModel.sampling_count_model,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    call FindTrickyVariants {
        input:
            variant_to_nodes = MergeChromosomeVariantToNodes.variant_to_nodes,
            sampling_count_model = MakeCountModel.sampling_count_model,
            reverse_variant_kmers = MakeReverseVariantKmerIndex.reverse_variant_kmers,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    call MakeIndexBundle {
        input:
            variant_to_nodes = MergeChromosomeVariantToNodes.variant_to_nodes,
            numpy_variants = MakeSitesOnlyVcfAndNumpyVariants.numpy_variants,
            refined_sampling_count_model = RefineCountModel.refined_sampling_count_model,
            tricky_variants = FindTrickyVariants.tricky_variants,
            helper_model = MakeHelperModel.helper_model,
            helper_model_combo_matrix = MakeHelperModel.helper_model_combo_matrix,
            kmer_index_only_variants_with_revcomp = MakeVariantKmerIndexWithReverseComplements.kmer_index_only_variants_with_revcomp,
            output_prefix = output_prefix_,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = medium_runtime_attributes
    }

    output {
        File index = MakeIndexBundle.index
        File kmer_index_only_variants_with_revcomp = MakeVariantKmerIndexWithReverseComplements.kmer_index_only_variants_with_revcomp

        File? preprocessed_panel_vcf_gz = ConcatVcfsPreprocessedPanel.vcf_gz
        File? preprocessed_panel_vcf_gz_tbi = ConcatVcfsPreprocessedPanel.vcf_gz_tbi
        File? preprocessed_panel_split_vcf_gz = ConcatVcfsPreprocessedPanelSplit.vcf_gz
        File? preprocessed_panel_split_vcf_gz_tbi = ConcatVcfsPreprocessedPanelSplit.vcf_gz_tbi
        File preprocessed_panel_bi_vcf_gz = preprocessed_panel_bi_vcf_gz_
        File preprocessed_panel_bi_vcf_gz_tbi = preprocessed_panel_bi_vcf_gz_tbi_
        File? preprocessed_panel_multi_split_vcf_gz = ConcatVcfsPreprocessedPanelMultiSplit.vcf_gz
        File? preprocessed_panel_multi_split_vcf_gz_tbi = ConcatVcfsPreprocessedPanelMultiSplit.vcf_gz_tbi
    }
}

task PreprocessPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        Array[String] chromosomes
        String output_prefix
        String? extra_view_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # use tee to pipe the output of the first bcftools command to the three subsequent command blocks
        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} ~{extra_view_args} -Ou | \
            bcftools norm --no-version -m+ -N -Ou | \
            bcftools plugin fill-tags --no-version -Ou -- -t AF,AC,AN | tee \
        >(
            # TODO check whether dropping AF=0 alleles here has any effect
            bcftools view --no-version --trim-alt-alleles -Ou | \
                bcftools view --no-version --min-alleles 2 -Ou | \
                bcftools plugin fill-tags --no-version -Ou -- -t AF,AC,AN | tee \
            >(
                bcftools norm --no-version -m- -N -Ou | \
                    bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.split.vcf.gz -- -t AF,AC,AN &&
                bcftools index -t ~{output_prefix}.preprocessed.split.vcf.gz
            ) | \
            bcftools view --no-version -Oz -o ~{output_prefix}.preprocessed.vcf.gz &&
            bcftools index -t ~{output_prefix}.preprocessed.vcf.gz
        ) \
        >(
            # we need to drop multiallelics before trimming, otherwise there may be representation issues in the graph
            bcftools view --no-version --min-alleles 2 --max-alleles 2  -Ou | \
                bcftools view --no-version --trim-alt-alleles -Ou | \
                bcftools view --no-version --min-alleles 2 -Ou | \
                bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.bi.vcf.gz -- -t AF,AC,AN &&
            bcftools index -t ~{output_prefix}.preprocessed.bi.vcf.gz
        ) | \
        (
            bcftools view --no-version --min-alleles 3 -Ou | \
                bcftools norm --no-version -m- -N -Ou | \
                bcftools view --no-version --trim-alt-alleles -Ou | \
                bcftools view --no-version --min-alleles 2 -Ou | \
                bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.multi.split.vcf.gz -- -t AF,AC,AN &&
            bcftools index -t ~{output_prefix}.preprocessed.multi.split.vcf.gz
        )
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_panel_vcf_gz = "~{output_prefix}.preprocessed.vcf.gz"
        File preprocessed_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.vcf.gz.tbi"
        File preprocessed_panel_split_vcf_gz = "~{output_prefix}.preprocessed.split.vcf.gz"
        File preprocessed_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.split.vcf.gz.tbi"
        File preprocessed_panel_bi_vcf_gz = "~{output_prefix}.preprocessed.bi.vcf.gz"
        File preprocessed_panel_bi_vcf_gz_tbi = "~{output_prefix}.preprocessed.bi.vcf.gz.tbi"
        File preprocessed_panel_multi_split_vcf_gz = "~{output_prefix}.preprocessed.multi.split.vcf.gz"
        File preprocessed_panel_multi_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.multi.split.vcf.gz.tbi"
    }
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    Int disk_size_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools concat ~{sep=" " vcf_gzs} --naive -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File vcf_gz = "~{output_prefix}.vcf.gz"
        File vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}

task MakeSitesOnlyVcfAndNumpyVariants {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # need to retain GT header lines
        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} | cut -f 1-9 - | bgzip > ~{output_prefix}.sites.vcf.gz
        bcftools index -t ~{output_prefix}.sites.vcf.gz

        obgraph make_numpy_variants \
            -v <(bgzip -c -d ~{output_prefix}.sites.vcf.gz) \
            -o ~{output_prefix}.numpy_variants.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File sites_only_vcf_gz = "~{output_prefix}.sites.vcf.gz"
        File sites_only_vcf_gz_tbi = "~{output_prefix}.sites.vcf.gz.tbi"
        File numpy_variants = "~{output_prefix}.numpy_variants.pkl"
    }
}

task MakeChromosomeGenotypeMatrix {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        lite_utils make_genotype_matrix \
            -v ~{input_vcf_gz} \
            -c ~{chromosome} \
            -o ~{output_prefix}.~{chromosome}.genotype_matrix.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_genotype_matrix = "~{output_prefix}.~{chromosome}.genotype_matrix.pkl"
    }
}

task MakeChromosomeGraph {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File reference_fasta
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        obgraph make \
            -r ~{reference_fasta} \
            -c ~{chromosome} \
            -v ~{input_vcf_gz} \
            -o ~{output_prefix}.~{chromosome}.obgraph.pkl

        obgraph add_allele_frequencies \
            -c ~{chromosome} \
            -v ~{input_vcf_gz} \
            -g ~{output_prefix}.~{chromosome}.obgraph.pkl

        obgraph make_position_id \
            -g ~{output_prefix}.~{chromosome}.obgraph.pkl \
            -o ~{output_prefix}.~{chromosome}.position_id_index.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_graph = "~{output_prefix}.~{chromosome}.obgraph.pkl"
        File chromosome_position_id_index = "~{output_prefix}.~{chromosome}.position_id_index.pkl"
    }
}

task MergeChromosomeGraphs {
    input {
        Array[File] chromosome_graphs
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        obgraph merge_graphs \
            -g ~{sep=" " chromosome_graphs} \
            -o ~{output_prefix}.obgraph.pkl

        obgraph make_position_id \
            -g ~{output_prefix}.obgraph.pkl \
            -o ~{output_prefix}.position_id_index.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File graph = "~{output_prefix}.obgraph.pkl"
        File position_id_index = "~{output_prefix}.position_id_index.pkl"
    }
}

task MakeChromosomeVariantToNodes {
    input {
        File chromosome_graph
        File sites_only_vcf_gz
        File sites_only_vcf_gz_tbi
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view \
            -Oz \
            -r ~{chromosome} \
            ~{sites_only_vcf_gz} \
            -o ~{chromosome}.sites.vcf.gz
        bcftools index -t ~{chromosome}.sites.vcf.gz

        obgraph make_variant_to_nodes \
            -g ~{chromosome_graph} \
            -v ~{chromosome}.sites.vcf.gz \
            -o ~{output_prefix}.~{chromosome}.variant_to_nodes.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_variant_to_nodes = "~{output_prefix}.~{chromosome}.variant_to_nodes.pkl"
    }
}

task MakeHelperModel {
    input {
        Array[File] chromosome_genotype_matrices
        File variant_to_nodes
        String output_prefix

        String docker
        File? monitoring_script
        Int? num_threads = 2 # single thread has bug
        Int? window_size = 100

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        lite_utils merge_genotype_matrices_and_convert_to_unphased \
            -g ~{sep=" " chromosome_genotype_matrices} \
            -o ~{output_prefix}.genotype_matrix.pkl

        kage create_helper_model \
            -g ~{output_prefix}.genotype_matrix.pkl \
            -v ~{variant_to_nodes} \
            -w ~{window_size} \
            -t ~{num_threads} \
            -o ~{output_prefix}.helper_model
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File helper_model = "~{output_prefix}.helper_model.pkl"
        File helper_model_combo_matrix = "~{output_prefix}.helper_model_combo_matrix.pkl"
    }
}

task SampleChromosomeKmersFromLinearReference {
    input {
        File reference_fasta
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script
        Int? spacing = 1
        Int? kmer_length = 31
        Boolean? include_reverse_complement = true

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # subset reference to chromosome
        samtools faidx ~{reference_fasta} \
            -r <(echo ~{chromosome}) \
            -o chromosome.fa
        samtools faidx chromosome.fa

        # get chromosome length from fasta index
        CHROMOSOME_LENGTH=$(cut -f 2 chromosome.fa.fai)

        NUM_THREADS=$(nproc)

        graph_kmer_index make \
            -t $NUM_THREADS \
            -s ~{spacing} \
            -k ~{kmer_length} \
            --include-reverse-complement ~{include_reverse_complement} \
            -R chromosome.fa \
            -n ~{chromosome} \
            -G "$CHROMOSOME_LENGTH" \
            -o ~{output_prefix}.~{chromosome}.linear_kmers.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 8])
        memory: select_first([runtime_attributes.command_mem_gb, 63]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_linear_kmers = "~{output_prefix}.~{chromosome}.linear_kmers.pkl"
    }
}

task MakeLinearReferenceKmerCounter {
    input {
        File linear_kmers
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        graph_kmer_index count_kmers \
             --subsample-ratio 1 \
            -f ~{linear_kmers} \
            -o ~{output_prefix}.linear_kmers_counter.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File linear_kmers_counter = "~{output_prefix}.linear_kmers_counter.pkl"
    }
}

task GetChromosomeShortVariantKmers {
    input {
        File chromosome_graph
        File chromosome_position_id_index
        File chromosome_variant_to_nodes
        File linear_kmers_counter
        File sites_only_vcf_gz
        File sites_only_vcf_gz_tbi
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31
        Int? chunk_size = 20000
        Int? max_variant_nodes = 3

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # check /dev/shm implementation between Batch and Life Sciences
        df -h /dev/shm

        bcftools view \
            -Oz \
            -r ~{chromosome} \
            ~{sites_only_vcf_gz} \
            -o ~{chromosome}.sites.vcf.gz
        bcftools index -t ~{chromosome}.sites.vcf.gz

        NUM_THREADS=$(nproc)

        graph_kmer_index make_unique_variant_kmers \
            -D true \
            -t $NUM_THREADS \
            -k ~{kmer_length} \
            -c ~{chunk_size} \
            --max-variant-nodes ~{max_variant_nodes} \
            -g ~{chromosome_graph} \
            -p ~{chromosome_position_id_index} \
            -V ~{chromosome_variant_to_nodes} \
            -I ~{linear_kmers_counter} \
            -v ~{chromosome}.sites.vcf.gz \
            -o ~{output_prefix}.~{chromosome}.short_variant_kmers.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 8])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_short_variant_kmers = "~{output_prefix}.~{chromosome}.short_variant_kmers.pkl"
    }
}

task SampleChromosomeStructuralVariantKmers {
    input {
        File chromosome_graph
        File chromosome_variant_to_nodes
        File linear_kmers_counter
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        graph_kmer_index sample_kmers_from_structural_variants \
            -k ~{kmer_length} \
            -g ~{chromosome_graph} \
            -V ~{chromosome_variant_to_nodes} \
            -i ~{linear_kmers_counter} \
            -o ~{output_prefix}.~{chromosome}.structural_variant_kmers.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_structural_variant_kmers = "~{output_prefix}.~{chromosome}.structural_variant_kmers.pkl"
    }
}

task MakeChromosomeHaplotypeToNodes {
    input {
        File chromosome_variant_to_nodes
        File chromosome_genotype_matrix
        String chromosome
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        obgraph make_haplotype_to_nodes_bnp \
            -g ~{chromosome_variant_to_nodes} \
            -v ~{chromosome_genotype_matrix} \
            -o ~{output_prefix}.~{chromosome}.haplotype_to_nodes.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File chromosome_haplotype_to_nodes = "~{output_prefix}.~{chromosome}.haplotype_to_nodes.pkl"
        File chromosome_haplotype_nodes = "~{output_prefix}.~{chromosome}.haplotype_to_nodes.pkl.haplotype_nodes"
    }
}

task MergeFlatKmers {
    input {
        Array[File] flat_kmers
        String output_prefix

        String docker
        File? monitoring_script
        File? num_nodes
        File? reference_fasta_fai

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        lite_utils merge_flat_kmers \
            --flat-kmers ~{sep=" " flat_kmers} \
            ~{"--num-nodes " + num_nodes} \
            ~{"--reference-fasta-fai " + reference_fasta_fai} \
            -o ~{output_prefix}.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File merged_kmers = "~{output_prefix}.pkl"
    }
}

task MergeChromosomeVariantToNodes {
    input {
        Array[File] chromosome_variant_to_nodes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        lite_utils merge_chromosome_variant_to_nodes \
            --variant-to-nodes ~{sep=" " chromosome_variant_to_nodes} \
            --output-prefix ~{output_prefix}
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File variant_to_nodes = "~{output_prefix}.variant_to_nodes.pkl"
        File num_nodes = "~{output_prefix}.num_nodes.pkl"
    }
}

task MergeChromosomeHaplotypeToNodes {
    input {
        Array[File] chromosome_haplotype_to_nodes
        Array[File] chromosome_haplotype_nodes
        File num_nodes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        lite_utils merge_chromosome_haplotype_to_nodes \
            --haplotype-to-nodes ~{sep=" " chromosome_haplotype_to_nodes} \
            --num-nodes ~{num_nodes} \
            --output-prefix ~{output_prefix}
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File haplotype_to_nodes = "~{output_prefix}.haplotype_to_nodes.pkl"
        File haplotype_nodes = "~{output_prefix}.haplotype_to_nodes.pkl.haplotype_nodes"
    }
}

task MakeReverseVariantKmerIndex {
    input {
        File variant_kmers
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        graph_kmer_index make_reverse \
            -f ~{variant_kmers} \
            -o ~{output_prefix}.reverse_variant_kmers.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File reverse_variant_kmers = "~{output_prefix}.reverse_variant_kmers.pkl"
    }
}

task MakeVariantKmerIndexWithReverseComplements {
    input {
        File variant_kmers
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31
        Int? hash_modulo = 200000033

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        graph_kmer_index make_from_flat \
            -f ~{variant_kmers} \
            --kmer-size ~{kmer_length} \
            --hash-modulo ~{hash_modulo} \
            --add-reverse-complements True \
            --skip-frequencies True \
            -o ~{output_prefix}.kmer_index_only_variants_with_revcomp.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kmer_index_only_variants_with_revcomp = "~{output_prefix}.kmer_index_only_variants_with_revcomp.pkl"
    }
}

task MakeCountModel {
    input {
        File graph
        File haplotype_to_nodes
        File haplotype_nodes
        File kmer_index_only_variants_with_revcomp
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kage sample_node_counts_from_population \
            -t ~{cpu_resolved} \
            -g ~{graph} \
            -H ~{haplotype_to_nodes} \
            -i ~{kmer_index_only_variants_with_revcomp} \
            -o ~{output_prefix}.sampling_count_model.pkl
    }

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File sampling_count_model = "~{output_prefix}.sampling_count_model.pkl"
    }
}

task RefineCountModel {
    input {
        File variant_to_nodes
        File sampling_count_model
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kage refine_sampling_model \
            -s ~{sampling_count_model} \
            -v ~{variant_to_nodes} \
            -o ~{output_prefix}.refined_sampling_count_model.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File refined_sampling_count_model = "~{output_prefix}.refined_sampling_count_model.pkl"
    }
}

task FindTrickyVariants {
    input {
        File variant_to_nodes
        File sampling_count_model
        File reverse_variant_kmers
        String output_prefix

        String docker
        File? monitoring_script
        Int? max_counts = 1000

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kage find_tricky_variants \
            -m ~{sampling_count_model} \
            -v ~{variant_to_nodes} \
            -r ~{reverse_variant_kmers} \
            -M ~{max_counts} \
            -o ~{output_prefix}.tricky_variants.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File tricky_variants = "~{output_prefix}.tricky_variants.pkl"
    }
}

task MakeIndexBundle {
    input {
        File variant_to_nodes
        File numpy_variants
        File refined_sampling_count_model
        File tricky_variants
        File helper_model
        File helper_model_combo_matrix
        File kmer_index_only_variants_with_revcomp
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kage make_index_bundle \
            -g ~{variant_to_nodes} \
            -v ~{numpy_variants} \
            -A ~{refined_sampling_count_model} \
            -x ~{tricky_variants} \
            -f ~{helper_model} \
            -F ~{helper_model_combo_matrix} \
            -i ~{kmer_index_only_variants_with_revcomp} \
            -o ~{output_prefix}.index.pkl
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File index = "~{output_prefix}.index.pkl"
    }
}
