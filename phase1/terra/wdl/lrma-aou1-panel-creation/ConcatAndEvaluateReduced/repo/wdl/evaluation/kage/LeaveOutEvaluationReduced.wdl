version 1.0

import "../../methods/kage/KAGEPanelWithPreprocessing.wdl" as KAGEPanelWithPreprocessing
import "../../methods/kage/KAGECasePerChromosomeFlexscattered.wdl" as KAGECasePerChromosome
import "../../methods/kage/GLIMPSEBatchedCasePerChromosome.wdl" as GLIMPSEBatchedCasePerChromosome
import "../../methods/pangenie/PanGenieIndex.wdl" as PanGenieIndex
import "../../methods/pangenie/PanGenieGenotype.wdl" as PanGenieGenotype

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

workflow LeaveOutEvaluation {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File case_reference_fasta
        File case_reference_fasta_fai
        File case_reference_dict
        File reference_fasta
        File reference_fasta_fai
        Array[File] genetic_maps
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        String output_prefix
        Array[String] chromosomes
        Array[String] leave_out_sample_names
        Int case_average_coverage
        Int kage_merge_batch_size
        Int glimpse_batch_size
        String? extra_chunk_args
        String? extra_phase_args
        Boolean do_pangenie
        String? extra_view_args

        # reduced arguments
        Int? num_short_variants_to_retain
        Boolean do_genotype_SVs = true

        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        Map[String, File] leave_out_crams

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        String docker
        String samtools_docker
        String kage_docker
        String pangenie_docker
        File? monitoring_script

        Int? cpu_make_count_model
        RuntimeAttributes? runtime_attributes
        RuntimeAttributes? medium_runtime_attributes
        RuntimeAttributes? large_runtime_attributes
        RuntimeAttributes? pangenie_runtime_attributes
        RuntimeAttributes? kage_runtime_attributes
        RuntimeAttributes? merge_runtime_attributes
        RuntimeAttributes? concat_runtime_attributes
        RuntimeAttributes? glimpse_phase_runtime_attributes
        RuntimeAttributes? glimpse_sample_runtime_attributes
        Map[String, Int]? chromosome_to_glimpse_command_mem_gb      # for running per-chromosome by choosing large chunk size; this will override glimpse_phase_runtime_attributes
        Int? glimpse_phase_preemptible
        RuntimeAttributes? calculate_metrics_runtime_attributes
    }

    call PreprocessPanelVCF {
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            repeat_mask_bed = repeat_mask_bed,
            segmental_duplications_bed = segmental_duplications_bed,
            simple_repeats_bed = simple_repeats_bed,
            challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            extra_view_args = extra_view_args,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    if (defined(num_short_variants_to_retain)) {
        call ReducePanelVCF {
            input:
                input_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                input_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                output_prefix = output_prefix,
                num_short_variants_to_retain = select_first([num_short_variants_to_retain]),
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }
    }

    String leave_out_output_prefix = output_prefix + ".LO"

    call CreateLeaveOneOutPanelVCF {
        input:
            input_vcf_gz = select_first([ReducePanelVCF.reduced_vcf_gz, PreprocessPanelVCF.preprocessed_panel_vcf_gz]),
            input_vcf_gz_tbi = select_first([ReducePanelVCF.reduced_vcf_gz_tbi, PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi]),
            output_prefix = leave_out_output_prefix,
            leave_out_sample_names = leave_out_sample_names,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    scatter (j in range(length(chromosomes))) {
        String chromosome = chromosomes[j]

        # we repeat preprocessing steps in these per-chromosome tasks
        call KAGEPanelWithPreprocessing.KAGEPanelWithPreprocessing as ChromosomeKAGELeaveOneOutPanel {
            input:
                input_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz,
                input_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz_tbi,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = leave_out_output_prefix,
                chromosomes = [chromosome],
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes,
                medium_runtime_attributes = medium_runtime_attributes,
                large_runtime_attributes = large_runtime_attributes,
                cpu_make_count_model = cpu_make_count_model
        }
    }

    if (do_pangenie) {
        call PanGenieIndex.PanGenieIndex as PanGenieIndex {
            input:
                panel_vcf_gz = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz,
                panel_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_out_panel_vcf_gz_tbi,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                chromosomes = chromosomes,
                output_prefix = leave_out_output_prefix,
                pangenie_docker = pangenie_docker,
                monitoring_script = monitoring_script,
                pangenie_runtime_attributes = pangenie_runtime_attributes
        }
    }

    scatter (j in range(length(leave_out_sample_names))) {
        String leave_out_sample_name = leave_out_sample_names[j]
        File leave_out_cram = leave_out_crams[leave_out_sample_name]

        call KAGECasePerChromosome.KAGECasePerChromosome as KAGECasePerChromosome {
            input:
                input_cram = leave_out_cram,
                sample_name = leave_out_sample_name,
                reference_fasta = case_reference_fasta,
                reference_fasta_fai = case_reference_fasta_fai,
                reference_dict = case_reference_dict,
                chromosomes = [chromosomes],
                panel_index = [ChromosomeKAGELeaveOneOutPanel.index],
                panel_kmer_index_only_variants_with_revcomp = [ChromosomeKAGELeaveOneOutPanel.kmer_index_only_variants_with_revcomp],
                panel_multi_split_vcf_gz = [select_all(ChromosomeKAGELeaveOneOutPanel.preprocessed_panel_multi_split_vcf_gz)],
                panel_multi_split_vcf_gz_tbi = [select_all(ChromosomeKAGELeaveOneOutPanel.preprocessed_panel_multi_split_vcf_gz_tbi)],
                average_coverage = case_average_coverage,
                samtools_docker = samtools_docker,
                kage_docker = kage_docker,
                monitoring_script = monitoring_script,
                kage_runtime_attributes = kage_runtime_attributes
        }

        if (!do_genotype_SVs) {
            scatter (c in range(length(KAGECasePerChromosome.chromosome_kage_vcf_gzs))) {
                call CensorGenotypes {
                    input:
                        input_vcf_gz = KAGECasePerChromosome.chromosome_kage_vcf_gzs[c],
                        input_vcf_gz_tbi = KAGECasePerChromosome.chromosome_kage_vcf_gz_tbis[c],
                        output_prefix = leave_out_sample_name,
                        retained_sv_tsv_gz = select_first([ReducePanelVCF.retained_sv_tsv_gz]),
                        retained_sv_tsv_gz_tbi = select_first([ReducePanelVCF.retained_sv_tsv_gz_tbi]),
                        docker = samtools_docker,
                        monitoring_script = monitoring_script
                }
            }
        }

        if (do_pangenie) {
            call PanGenieGenotype.PanGenieGenotype as PanGenieGenotype {
                input:
                    pangenie_index_chromosome_graphs = select_first([PanGenieIndex.pangenie_index_chromosome_graphs]),
                    pangenie_index_chromosome_kmers = select_first([PanGenieIndex.pangenie_index_chromosome_kmers]),
                    pangenie_index_unique_kmers_map = select_first([PanGenieIndex.pangenie_index_unique_kmers_map]),
                    pangenie_index_path_segments_fasta = select_first([PanGenieIndex.pangenie_index_path_segments_fasta]),
                    index_prefix = leave_out_output_prefix,
                    reference_fasta = case_reference_fasta,
                    reference_fasta_fai = reference_fasta_fai,
                    chromosomes = chromosomes,
                    input_cram = leave_out_cram,
                    sample_name = leave_out_sample_name,
                    docker = docker,
                    pangenie_docker = pangenie_docker,
                    monitoring_script = monitoring_script,
                    pangenie_runtime_attributes = pangenie_runtime_attributes
            }

            # PanGenie evaluation
            call CalculateMetrics as CalculateMetricsPanGenie {
                input:
                    case_vcf_gz = PanGenieGenotype.genotyping_vcf_gz,
                    case_vcf_gz_tbi = PanGenieGenotype.genotyping_vcf_gz_tbi,
                    truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                    truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                    chromosomes = chromosomes,
                    label = "PanGenie",
                    sample_name = leave_out_sample_name,
                    docker = docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = calculate_metrics_runtime_attributes
            }
        }
    }

    call WriteTsv as WriteTsvVcfs {
        input:
            array = if !do_genotype_SVs then select_all(CensorGenotypes.censored_vcf_gz) else KAGECasePerChromosome.chromosome_kage_vcf_gzs,
            docker = docker
    }

    call WriteTsv as WriteTsvTbis {
        input:
            array = if !do_genotype_SVs then select_all(CensorGenotypes.censored_vcf_gz_tbi) else KAGECasePerChromosome.chromosome_kage_vcf_gz_tbis,
            docker = docker
    }

    call WriteLines as WriteLinesSampleNames {
        input:
            array = leave_out_sample_names,
            docker = docker
    }

    call GLIMPSEBatchedCasePerChromosome.GLIMPSEBatchedCasePerChromosome as GLIMPSEBatchedCasePerChromosome {
        input:
            sample_by_chromosome_kage_vcf_gzs_tsv = WriteTsvVcfs.tsv,
            sample_by_chromosome_kage_vcf_gz_tbis_tsv = WriteTsvTbis.tsv,
            sample_names_file = WriteLinesSampleNames.tsv,
            chromosomes = chromosomes,
            genetic_maps = genetic_maps,
            panel_split_vcf_gz = select_all(ChromosomeKAGELeaveOneOutPanel.preprocessed_panel_split_vcf_gz),
            panel_split_vcf_gz_tbi = select_all(ChromosomeKAGELeaveOneOutPanel.preprocessed_panel_split_vcf_gz_tbi),
            kage_merge_batch_size = kage_merge_batch_size,
            glimpse_batch_size = glimpse_batch_size,
            extra_chunk_args = extra_chunk_args,
            extra_phase_args = extra_phase_args,
            output_prefix = "leave-out.batch",
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            kage_docker = kage_docker,
            monitoring_script = monitoring_script,
            merge_runtime_attributes = merge_runtime_attributes,
            concat_runtime_attributes = concat_runtime_attributes,
            glimpse_phase_runtime_attributes = glimpse_phase_runtime_attributes,
            glimpse_sample_runtime_attributes = glimpse_sample_runtime_attributes,
            chromosome_to_glimpse_command_mem_gb = chromosome_to_glimpse_command_mem_gb,
            glimpse_phase_preemptible = glimpse_phase_preemptible
    }

    scatter (j in range(length(leave_out_sample_names))) {
        # KAGE evaluation
        call CalculateMetrics as CalculateMetricsKAGE {
            input:
                case_vcf_gz = GLIMPSEBatchedCasePerChromosome.kage_vcf_gz,
                case_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.kage_vcf_gz_tbi,
                truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                chromosomes = chromosomes,
                label = "KAGE",
                sample_name = leave_out_sample_names[j],
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = calculate_metrics_runtime_attributes
        }

        # KAGE+GLIMPSE evaluation
        call CalculateMetrics as CalculateMetricsKAGEPlusGLIMPSE {
            input:
                case_vcf_gz = GLIMPSEBatchedCasePerChromosome.glimpse_vcf_gz,
                case_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.glimpse_vcf_gz_tbi,
                truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                chromosomes = chromosomes,
                label = "KAGE+GLIMPSE",
                sample_name = leave_out_sample_names[j],
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = calculate_metrics_runtime_attributes
        }

        # KAGE+GLIMPSE+FixVariantCollisions evaluation
        call CalculateMetrics as CalculateMetricsPhasedCollisionless {
            input:
                case_vcf_gz = GLIMPSEBatchedCasePerChromosome.phased_collisionless_vcf_gz,
                case_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.phased_collisionless_vcf_gz_tbi,
                truth_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                truth_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                chromosomes = chromosomes,
                label = "KAGE+GLIMPSE",
                sample_name = leave_out_sample_names[j],
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = calculate_metrics_runtime_attributes
        }
    }

    output {
        # merged
        File kage_vcf_gz = GLIMPSEBatchedCasePerChromosome.kage_vcf_gz
        File kage_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.kage_vcf_gz_tbi
        File glimpse_unphased_vcf_gz = GLIMPSEBatchedCasePerChromosome.glimpse_unphased_vcf_gz
        File glimpse_unphased_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.glimpse_unphased_vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSEBatchedCasePerChromosome.glimpse_vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.glimpse_vcf_gz_tbi
        File phased_collisionless_vcf_gz = GLIMPSEBatchedCasePerChromosome.phased_collisionless_vcf_gz
        File phased_collisionless_vcf_gz_tbi = GLIMPSEBatchedCasePerChromosome.phased_collisionless_vcf_gz_tbi
        # per-sample
        Array[File?] pangenie_vcf_gzs = PanGenieGenotype.genotyping_vcf_gz
        Array[File?] pangenie_vcf_gz_tbis = PanGenieGenotype.genotyping_vcf_gz_tbi

        Array[File] kage_metrics_tsvs = CalculateMetricsKAGE.metrics_tsv
        Array[File] glimpse_metrics_tsvs = CalculateMetricsKAGEPlusGLIMPSE.metrics_tsv
        Array[File] phased_collisionless_metrics_tsvs = CalculateMetricsPhasedCollisionless.metrics_tsv
        Array[File?] pangenie_metrics_tsvs = CalculateMetricsPanGenie.metrics_tsv
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

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
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

task PreprocessPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
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

        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} ~{extra_view_args} -Ou | \
            bcftools norm --no-version -m+ -N -Ou | \
            bcftools plugin fill-tags --no-version -Ou -- -t AF,AC,AN | \
            truvari anno svinfo | \
            bcftools annotate --no-version -a ~{repeat_mask_bed} -c CHROM,FROM,TO -m +RM -Ou | \
            bcftools annotate --no-version -a ~{segmental_duplications_bed} -c CHROM,FROM,TO -m +SD -Ou | \
            bcftools annotate --no-version -a ~{simple_repeats_bed} -c CHROM,FROM,TO -m +SR -Ou | \
            bcftools annotate --no-version -a ~{challenging_medically_relevant_genes_bed} -c CHROM,FROM,TO -m +CMRG -Oz -o ~{output_prefix}.preprocessed.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.vcf.gz

        bcftools view --no-version --min-alleles 3  ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools norm --no-version -m- -N -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.multi.split.vcf.gz

        bcftools norm --no-version -m- -N ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.split.temp.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.split.temp.vcf.gz

        bcftools annotate --no-version -a ~{output_prefix}.preprocessed.multi.split.vcf.gz ~{output_prefix}.preprocessed.split.temp.vcf.gz -m +MULTIALLELIC \
            -Oz -o ~{output_prefix}.preprocessed.split.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.split.vcf.gz
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
    }
}

task ReducePanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String output_prefix
        Int num_short_variants_to_retain

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

        bcftools view ~{input_vcf_gz} | \
            grep -v SVLEN | \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' | \
            shuf | \
            head -n ~{num_short_variants_to_retain} | \
            sort -k1V -k2h | \
            bgzip -c > ~{output_prefix}.retained.short.tsv.gz
        tabix -s1 -b2 -e2 ~{output_prefix}.retained.short.tsv.gz

        bcftools view ~{input_vcf_gz} | \
            grep -E 'SVLEN|#' | \
            bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' | \
            bgzip -c > ~{output_prefix}.retained.sv.tsv.gz
        tabix -s1 -b2 -e2 ~{output_prefix}.retained.sv.tsv.gz

        cat <(bgzip -cd ~{output_prefix}.retained.short.tsv.gz) <(bgzip -cd ~{output_prefix}.retained.sv.tsv.gz) | \
            sort -k1V -k2h | \
            bgzip -c > ~{output_prefix}.retained.tsv.gz
        tabix -s1 -b2 -e2 ~{output_prefix}.retained.tsv.gz

        bcftools norm --no-version -m+any -N ~{input_vcf_gz} \
            -T ~{output_prefix}.retained.tsv.gz \
            -Oz -o ~{output_prefix}.reduced.vcf.gz
        bcftools index -t ~{output_prefix}.reduced.vcf.gz
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
        File retained_short_tsv_gz = "~{output_prefix}.retained.short.tsv.gz"
        File retained_short_tsv_gz_tbi = "~{output_prefix}.retained.short.tsv.gz.tbi"
        File retained_sv_tsv_gz = "~{output_prefix}.retained.sv.tsv.gz"
        File retained_sv_tsv_gz_tbi = "~{output_prefix}.retained.sv.tsv.gz.tbi"
        File retained_tsv_gz = "~{output_prefix}.retained.tsv.gz"
        File retained_tsv_gz_tbi = "~{output_prefix}.retained.tsv.gz.tbi"
        File reduced_vcf_gz = "~{output_prefix}.reduced.vcf.gz"
        File reduced_vcf_gz_tbi = "~{output_prefix}.reduced.vcf.gz.tbi"
    }
}

task CreateLeaveOneOutPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String output_prefix
        Array[String] leave_out_sample_names

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

        bcftools view --no-version ~{input_vcf_gz} -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.vcf.gz

        bcftools norm --no-version -m- -N ~{output_prefix}.preprocessed.LO.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.split.vcf.gz

        # we need to drop multiallelics before LO and trimming, otherwise there may be representation issues in the graph
        bcftools view --no-version --min-alleles 2 --max-alleles 2  ~{input_vcf_gz} -Ou | \
            bcftools view --no-version -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.bi.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.bi.vcf.gz

         bcftools view --no-version --min-alleles 3  ~{input_vcf_gz} -Ou | \
            bcftools norm --no-version -m- -N -Ou | \
            bcftools view --no-version -s ^~{sep=',' leave_out_sample_names} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LO.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LO.multi.split.vcf.gz
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
        File leave_out_panel_vcf_gz = "~{output_prefix}.preprocessed.LO.vcf.gz"
        File leave_out_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.vcf.gz.tbi"
        File leave_out_panel_split_vcf_gz = "~{output_prefix}.preprocessed.LO.split.vcf.gz"
        File leave_out_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.split.vcf.gz.tbi"
        File leave_out_panel_bi_vcf_gz = "~{output_prefix}.preprocessed.LO.bi.vcf.gz"
        File leave_out_panel_bi_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.bi.vcf.gz.tbi"
        File leave_out_panel_multi_split_vcf_gz = "~{output_prefix}.preprocessed.LO.multi.split.vcf.gz"
        File leave_out_panel_multi_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LO.multi.split.vcf.gz.tbi"
    }
}

# set SV GTs to missing
# TODO further subsample retained short variants and set remainder of GTs to missing
task CensorGenotypes {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String output_prefix
        File retained_sv_tsv_gz
        File retained_sv_tsv_gz_tbi

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

        bcftools +setGT -T ~{retained_sv_tsv_gz} --targets-overlap 1 ~{input_vcf_gz} -Oz -o ~{output_prefix}.censored.sv.vcf.gz -- -t a -n .    # TODO set GLs to nan
        bcftools view -T ^~{retained_sv_tsv_gz} --targets-overlap 1 ~{input_vcf_gz} -Oz -o ~{output_prefix}.retained.short.vcf.gz

        bcftools concat ~{output_prefix}.censored.sv.vcf.gz ~{output_prefix}.retained.short.vcf.gz | bcftools sort -Oz -o ~{output_prefix}.censored.vcf.gz
        bcftools index -t ~{output_prefix}.censored.vcf.gz
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
        # File genotyped_short_tsv_gz = "~{output_prefix}.genotyped.short.tsv.gz"
        # File genotyped_short_tsv_gz_tbi = "~{output_prefix}.genotyped.short.tsv.gz.tbi"
        File censored_vcf_gz = "~{output_prefix}.censored.vcf.gz"
        File censored_vcf_gz_tbi = "~{output_prefix}.censored.vcf.gz.tbi"
    }
}

task CalculateMetrics {
    input {
        File case_vcf_gz        # bi+multi split
        File case_vcf_gz_tbi
        File truth_vcf_gz       # bi+multi split
        File truth_vcf_gz_tbi
        Array[String] chromosomes
        String label
        String sample_name

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # split multiallelics in case (may be redundant)
        bcftools view --no-version -r ~{sep="," chromosomes} -s ~{sample_name} ~{case_vcf_gz} | bcftools norm --no-version -m- -N -Oz -o case.split.vcf.gz
        bcftools index -t case.split.vcf.gz

        # mark case variants in panel
        bcftools annotate --no-version -r ~{sep="," chromosomes} -a case.split.vcf.gz -m +CASE ~{truth_vcf_gz} -Oz -o panel.annot.vcf.gz
        bcftools index -t panel.annot.vcf.gz

        # TODO old conda install freezing upon solve, we'll mix conda and pip for now
        # conda install -y seaborn numpy=1.26.4
        conda install -y numpy=1.26.4
        pip install seaborn==0.13.2

        python - --case_vcf_gz case.split.vcf.gz \
                 --truth_vcf_gz panel.annot.vcf.gz \
                 --label ~{label} \
                 --sample_name ~{sample_name} \
                 <<-'EOF'
        import argparse
        import allel
        import numpy as np
        import pandas as pd
        import sklearn.metrics
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns

        matplotlib.use('Agg')

        def load_callset(callset_file, **kwargs):
            return allel.read_vcf(callset_file, **kwargs)

        def get_gt_vsp(callset):
            return allel.GenotypeDaskArray(callset['calldata/GT'])

        def calculate_metrics_and_plot(case_vcf_gz, truth_vcf_gz, label, sample_name, average='macro'):
            samples = [sample_name]

            callset = load_callset(case_vcf_gz, samples=samples, fields='*', alt_number=1)
            gt_vp = get_gt_vsp(callset)[:, 0, :].compute()

            truth_callset = load_callset(truth_vcf_gz, samples=samples, fields='*', alt_number=1)
            is_case_V = truth_callset['variants/CASE']
            truth_gt_vp = get_gt_vsp(truth_callset)[is_case_V, 0, :].compute()

            is_any_missing_case_v = np.any(gt_vp == -1, axis=1)
            is_any_missing_truth_v = np.any(truth_gt_vp == -1, axis=1)
            is_missing_v = is_any_missing_case_v | is_any_missing_truth_v

            is_snp_v = truth_callset['variants/is_snp'][is_case_V]
            is_multiallelic_v = truth_callset['variants/MULTIALLELIC'][is_case_V]
            altfreq_v = truth_callset['variants/AF'][is_case_V]
            is_altfreq_v = [['[0%, 1%)', (0. <= altfreq_v) & (altfreq_v < 0.01)],
                            ['[1%, 5%)', (0.01 <= altfreq_v) & (altfreq_v < 0.05)],
                            ['[5%, 10%)', (0.05 <= altfreq_v) & (altfreq_v < 0.1)],
                            ['[10%, 50%)', (0.1 <= altfreq_v) & (altfreq_v < 0.50)],
                            ['[50%, 100%]', (0.5 <= altfreq_v) & (altfreq_v <= 1.)]]
            altlen_v = truth_callset['variants/altlen'][is_case_V]
            is_altlen_v = [['(-inf,-500]', altlen_v <= -500],
                           ['(-500,-50]', (-500 < altlen_v) & (altlen_v <= -50)],
                           ['(-50,-1]', (-50 < altlen_v) & (altlen_v <= -1)],
                           ['0 (SNP)', is_snp_v],
                           ['[1,50)', (1 <= altlen_v) & (altlen_v < 50)],
                           ['[50,500)', (50 <= altlen_v) & (altlen_v < 500)],
                           ['[500,5000)', (500 <= altlen_v) & (altlen_v < 5000)],
                           ['[5000,inf)', 5000 <= altlen_v]]
            is_sv_v = (altlen_v <= -50) | (altlen_v >= 50)
            num_i = len(is_altfreq_v)
            num_j = len(is_altlen_v)

            metrics_dicts = []

            for context in ['ALL', 'CMRG', 'US', 'RM', 'SD', 'SR']:

                if context == 'ALL':
                    is_context_v = True
                elif context == 'US':
                    is_context_v = ~(truth_callset[f'variants/RM'][is_case_V] |
                                     truth_callset[f'variants/SD'][is_case_V] |
                                     truth_callset[f'variants/SR'][is_case_V])
                else:
                    is_context_v = truth_callset[f'variants/{context}'][is_case_V]

                for allelic in ['biallelic+multiallelic', 'multiallelic', 'biallelic']:
                    num_evals = np.zeros((num_i, num_j))
                    precisions = np.zeros((num_i, num_j))
                    recalls = np.zeros((num_i, num_j))
                    f1s = np.zeros((num_i, num_j))
                    confusion_matrices = np.zeros((num_i, num_j, 3, 3))

                    if allelic == 'biallelic+multiallelic':
                        is_allelic_v = True
                    elif allelic == 'multiallelic':
                        is_allelic_v = is_multiallelic_v
                    else:
                        is_allelic_v = ~is_multiallelic_v

                    truth_missing_count = (is_any_missing_truth_v & is_allelic_v & is_context_v).sum()
                    case_missing_count = (is_any_missing_case_v & is_allelic_v & is_context_v).sum()
                    missing_count = (is_missing_v & is_allelic_v & is_context_v).sum()

                    for i, (filter_name_i, is_v_i) in enumerate(is_altfreq_v):
                        for j, (filter_name_j, is_v_j) in enumerate(is_altlen_v):
                            # is_eval_v = ~is_missing_v & is_v_i & is_v_j & is_allelic_v & is_context_v         # (old behavior: evaluate over non-missing intersection of truth and case)
                            is_eval_v = ~is_any_missing_truth_v & is_v_i & is_v_j & is_allelic_v & is_context_v # evaluate over non-missing truth
                            gt_vp[is_any_missing_case_v, :] = 0                                                 # set missing in case to hom-ref

                            num_eval = np.sum(is_eval_v)
                            enc_gt_n = np.sum(gt_vp[is_eval_v], axis=1)
                            truth_enc_gt_n = np.sum(truth_gt_vp[is_eval_v], axis=1)
                            if num_eval == 0:
                                precision = np.nan
                                recall = np.nan
                                f1 = np.nan
                            else:
                                precision = sklearn.metrics.precision_score(truth_enc_gt_n, enc_gt_n, average=average)
                                recall = sklearn.metrics.recall_score(truth_enc_gt_n, enc_gt_n, average=average)
                                f1 = sklearn.metrics.f1_score(truth_enc_gt_n, enc_gt_n, average=average)
                            confusion_matrix = sklearn.metrics.confusion_matrix(truth_enc_gt_n, enc_gt_n, labels=[0, 1, 2])

                            num_evals[i][j] = num_eval
                            precisions[i][j] = precision
                            recalls[i][j] = recall
                            f1s[i][j] = f1
                            confusion_matrices[i][j] = confusion_matrix

                            metrics_dicts.append({
                                'LABEL': label,
                                'SAMPLE_NAME': sample_name,
                                'CONTEXT': context,
                                'ALLELIC': allelic,
                                'ALTFREQ': filter_name_i,
                                'ALTLEN': filter_name_j,
                                'NUM_EVAL': num_eval,
                                'PRECISION': precision,
                                'RECALL': recall,
                                'F1': f1,
                                'CONFUSION_MATRIX': confusion_matrix
                            })

                    if num_evals.sum() == 0:
                        continue

                    non_sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_truth_enc_gt_n = np.sum(truth_gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_precision = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.precision_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_recall = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.recall_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_f1 = np.nan if non_sv_truth_enc_gt_n.size == 0 else sklearn.metrics.f1_score(non_sv_truth_enc_gt_n, non_sv_enc_gt_n, average=average)
                    non_sv_count = np.sum(~is_sv_v & is_allelic_v & is_context_v)

                    sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_truth_enc_gt_n = np.sum(truth_gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_precision = sklearn.metrics.precision_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_recall = sklearn.metrics.recall_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_f1 = sklearn.metrics.f1_score(sv_truth_enc_gt_n, sv_enc_gt_n, average=average)
                    sv_count = np.sum(is_sv_v & is_allelic_v & is_context_v)

        #             fig, ax = plt.subplots(5, 1, figsize=(12, 20))
                    fig, ax = plt.subplots(4, 1, figsize=(12, 16))

                    ax[0] = sns.heatmap(num_evals, ax=ax[0], linewidths=1, linecolor='k', annot=True,
                                        norm=matplotlib.colors.LogNorm(), cmap='Blues')
                    cbar = ax[0].collections[0].colorbar
                    cbar.ax.set_ylabel('number of alt alleles', rotation=270, labelpad=20)
                    ax[0].set_title(f'{label}\n{sample_name}\ncontext = {context}, {allelic}\n\n' +
                                    f'alt allele count\n' +
                                    f'non-SV = {non_sv_count}, SV = {sv_count}\n' +
                                    f'(missing: truth = {truth_missing_count}, case = {case_missing_count}, union = {missing_count})')

        #             diag_mask = np.tile(np.eye(3), (num_i, num_j)).astype(bool)
        #             ax[1] = sns.heatmap(np.reshape(confusion_matrices, (3 * num_i, 3 * num_j)),
        #                                 ax=ax[1], linecolor='k',
        # #                                 norm=matplotlib.colors.LogNorm(),
        #                                 mask=~diag_mask, cmap='Greens')
        #             ax[1] = sns.heatmap(np.reshape(confusion_matrices, (3 * num_i, 3 * num_j)),
        #                                 ax=ax[1], linecolor='k',
        # #                                 norm=matplotlib.colors.LogNorm(),
        #                                 mask=diag_mask, cmap='Reds', cbar=False)
        #             cbar = ax[1].collections[1].colorbar
        # #             cbar.ax.set_ylabel(f'number of alt alleles', rotation=270, labelpad=40)
        #             ax[1].set_title(f'normalized confusion matrices')

                    ax[1] = sns.heatmap(precisions, ax=ax[1], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[1].collections[0].colorbar
                    cbar.ax.set_ylabel(f'precision ({average})', rotation=270, labelpad=20)
                    ax[1].set_title(f'precision ({average})\n' +
                                    f'non-SV = {non_sv_precision:.4f}, SV = {sv_precision:.4f}')

                    ax[2] = sns.heatmap(recalls, ax=ax[2], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[2].collections[0].colorbar
                    cbar.ax.set_ylabel(f'recall ({average})', rotation=270, labelpad=20)
                    ax[2].set_title(f'recall ({average})\n' +
                                    f'non-SV = {non_sv_recall:.4f}, SV = {sv_recall:.4f}')

                    ax[3] = sns.heatmap(f1s, ax=ax[3], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.5, cmap='Greens')
                    cbar = ax[3].collections[0].colorbar
                    cbar.ax.set_ylabel(f'F1 ({average})', rotation=270, labelpad=20)
                    ax[3].set_title(f'F1 ({average})\n' +
                                    f'non-SV = {non_sv_f1:.4f}, SV = {sv_f1:.4f}')

                    for i, a in enumerate(ax):
                        a.set_ylabel('AF')
                        a.tick_params(bottom=False, left=False)
                        a.set_xticklabels([])
        #                 if i != 1:
                        a.set_yticklabels([filter_name for filter_name, _ in is_altfreq_v], rotation=0)
        #                 else:
        #                     a.set_yticklabels(sum([[None, filter_name, None]
        #                                            for filter_name, _ in is_altfreq_v], []),
        #                                       rotation=0)
                    ax[3].set_xlabel('ALT length - REF length (bp)')
                    ax[3].set_xticklabels([filter_name for filter_name, _ in is_altlen_v], rotation=0)

                    plt.savefig(f'{sample_name}.{context}.{allelic}.{label}.metrics.png')
                    plt.show()

            metrics_df = pd.DataFrame.from_dict(metrics_dicts)
            metrics_df.to_csv(f'{sample_name}.{label}.metrics.tsv', sep='\t', index=False)


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--case_vcf_gz',
                                type=str)

            parser.add_argument('--truth_vcf_gz',
                                type=str)

            parser.add_argument('--label',
                                type=str)

            parser.add_argument('--sample_name',
                                type=str)

            args = parser.parse_args()

            calculate_metrics_and_plot(args.case_vcf_gz,
                                       args.truth_vcf_gz,
                                       args.label,
                                       args.sample_name)

        if __name__ == '__main__':
            main()
        EOF
    >>>

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
        File metrics_tsv = "~{sample_name}.~{label}.metrics.tsv"
        Array[File] metrics_plots = glob("~{sample_name}.*.png")
    }
}

task WriteTsv {
    input {
        Array[Array[String]] array
        String docker
    }

    command <<<
    >>>

    output {
        File tsv = write_tsv(array)
    }

    runtime {
        docker: docker
    }
}

task WriteLines {
    input {
        Array[String] array
        String docker
    }

    command <<<
    >>>

    output {
        File tsv = write_lines(array)
    }

    runtime {
        docker: docker
    }
}
