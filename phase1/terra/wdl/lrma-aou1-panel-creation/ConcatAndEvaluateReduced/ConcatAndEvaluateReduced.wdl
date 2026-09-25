########################################################################################
# ConcatAndEvaluateReduced (the WDL names its workflow PhasedPanelEvaluation)
# Phase 1 stage c (reference panel), launch step 7 (after the per-chromosome step).
#
# Concatenates the per-chromosome outputs of ChromosomePhasedPanelCreationFromHiPhase
# into whole-genome VCFs, including the final panel VCF, and evaluates the panel:
# leave-out genotyping of held-out samples from their short-read CRAMs (KAGE and
# GLIMPSE, optionally PanGenie), and vcfdist and overlap metrics, stratified by
# repeat, segmental-duplication, and medically relevant gene regions.
#
# Runs once per cohort. This copy carries a local patch; see SOURCE.txt.
# Inputs:  the per-chromosome arrays from ChromosomePhasedPanelCreationFromHiPhase,
#          leave_out_sample_names / leave_out_crams / leave_out_chromosomes,
#          stratification BEDs, reference files, Docker images
# Outputs: panel_vcf_gz  the whole-genome panel VCF
#          evaluation_summary_tsvs
# Previous: ChromosomePhasedPanelCreationFromHiPhase. Next: KAGEPanelWithPreprocessing.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/evaluation/phased-panel/VcfdistAndOverlapMetricsEvaluation.wdl"
import "repo/wdl/evaluation/kage/LeaveOutEvaluationReduced.wdl"

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow PhasedPanelEvaluation {    # TODO change name later, easier to share configs for now

    input {
        # common inputs
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String output_prefix
        File? monitoring_script

        File? hiphase_short_vcf_gz      # this can be expensive to evaluate, make it optional
        File? hiphase_short_vcf_gz_tbi
        File hiphase_sv_vcf_gz
        File hiphase_sv_vcf_gz_tbi

        Boolean do_short_read_stages

        # per chromosome inputs
        Array[File] filter_and_concat_vcf_gzs
        Array[File] filter_and_concat_vcf_gz_tbis
        Array[File] before_shapeit4_collisionless_bcfs
        Array[File] before_shapeit4_collisionless_bcf_csis
        Array[File] phased_vcf_gzs
        Array[File] phased_vcf_gz_tbis
        Array[File] collisionless_bcfs
        Array[File] collisionless_bcf_csis
        Array[File] panel_vcf_gzs
        Array[File] panel_vcf_gz_tbis

        # inputs for LeaveOutEvaluation
        File case_reference_fasta
        File case_reference_fasta_fai
        File case_reference_dict
        Array[File] leave_out_genetic_maps
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        Array[String] leave_out_chromosomes
        Array[String] leave_out_sample_names
        Int case_average_coverage
        Int kage_merge_batch_size
        Int glimpse_batch_size
        String? extra_chunk_args
        String? extra_phase_args
        Boolean do_pangenie
        String? extra_view_args
        Map[String, File] leave_out_crams
        String leave_out_docker
        String kage_docker
        String samtools_docker
        String pangenie_docker
        Int? cpu_make_count_model
        RuntimeAttributes? leave_out_runtime_attributes
        RuntimeAttributes? leave_out_medium_runtime_attributes
        RuntimeAttributes? leave_out_large_runtime_attributes
        RuntimeAttributes? pangenie_runtime_attributes
        RuntimeAttributes? kage_runtime_attributes
        RuntimeAttributes? merge_runtime_attributes
        RuntimeAttributes? concat_runtime_attributes
        RuntimeAttributes? glimpse_phase_runtime_attributes
        RuntimeAttributes? glimpse_sample_runtime_attributes
        Map[String, Int]? chromosome_to_glimpse_command_mem_gb      # for running per-chromosome by choosing large chunk size; this will override glimpse_phase_runtime_attributes
        Int? glimpse_phase_preemptible
        RuntimeAttributes? calculate_metrics_runtime_attributes

        # reduced arguments
        Int? num_short_variants_to_retain
        Boolean do_genotype_SVs = true

        # inputs for Ivcfmerge
        String ivcfmerge_docker

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        # inputs for VcfdistAndOverlapMetricsEvaluation
        Array[String] vcfdist_samples
        Array[File] confident_regions_bed_files
        File vcfdist_truth_vcf
        File vcfdist_truth_vcf_idx
        String evaluation_chromosomes_regions_arg
        Array[File] vcfdist_bed_file_array
        Array[String?] vcfdist_extra_args_array       # one for each element in vcfdist_bed_file_array
        Int? vcfdist_mem_gb
        String overlap_metrics_docker

        String summarize_evaluations_docker
    }

    call ConcatVcfs as ConcatVcfsFilterAndConcatVcfs { input:
        vcf_gzs = filter_and_concat_vcf_gzs,
        vcf_gz_tbis = filter_and_concat_vcf_gz_tbis,
        output_prefix = output_prefix + ".filter_and_concat",
        extra_concat_args = "--naive",
        docker = samtools_docker,
        monitoring_script = monitoring_script
    }

#    # cannot naive concat .bcf to .vcf.gz
#    call ConcatVcfs as ConcatVcfsBeforeShapeit4FixVariantCollisions { input:
#        vcf_gzs = before_shapeit4_collisionless_bcfs,
#        vcf_gz_tbis = before_shapeit4_collisionless_bcf_csis,
#        output_prefix = output_prefix + ".before_shapeit4_collisionless",
#        docker = samtools_docker,
#        monitoring_script = monitoring_script
#    }

    call ConcatVcfs as ConcatVcfsShapeit4 { input:
        vcf_gzs = phased_vcf_gzs,
        vcf_gz_tbis = phased_vcf_gz_tbis,
        output_prefix = output_prefix + ".phased",
        extra_concat_args = "--naive",
        docker = samtools_docker,
        monitoring_script = monitoring_script
    }

#    # cannot naive concat .bcf to .vcf.gz
#    call ConcatVcfs as ConcatVcfsFixVariantCollisions { input:
#        vcf_gzs = collisionless_bcfs,
#        vcf_gz_tbis = collisionless_bcf_csis,
#        output_prefix = output_prefix + ".phased.collisionless",
#        docker = samtools_docker,
#        monitoring_script = monitoring_script
#    }

    call ConcatVcfs as ConcatVcfsPanGeniePanelCreation { input:
        vcf_gzs = panel_vcf_gzs,
        vcf_gz_tbis = panel_vcf_gz_tbis,
        output_prefix = output_prefix + ".panel",
        extra_concat_args = "--naive",
        docker = samtools_docker,
        monitoring_script = monitoring_script
    }

    if (do_short_read_stages) {
        call LeaveOutEvaluationReduced.LeaveOutEvaluation { input:
            input_vcf_gz = ConcatVcfsPanGeniePanelCreation.vcf_gz,
            input_vcf_gz_tbi = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi,
            case_reference_fasta = case_reference_fasta,
            case_reference_fasta_fai = case_reference_fasta_fai,
            case_reference_dict = case_reference_dict,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            genetic_maps = leave_out_genetic_maps,
            repeat_mask_bed = repeat_mask_bed,
            segmental_duplications_bed = segmental_duplications_bed,
            simple_repeats_bed = simple_repeats_bed,
            challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
            output_prefix = output_prefix,
            chromosomes = leave_out_chromosomes,
            leave_out_sample_names = leave_out_sample_names,
            case_average_coverage = case_average_coverage,
            kage_merge_batch_size = kage_merge_batch_size,
            glimpse_batch_size = glimpse_batch_size,
            extra_chunk_args = extra_chunk_args,
            extra_phase_args = extra_phase_args,
            do_pangenie = do_pangenie,
            extra_view_args = extra_view_args,
            leave_out_crams = leave_out_crams,
            docker = leave_out_docker,
            samtools_docker = samtools_docker,
            kage_docker = kage_docker,
            pangenie_docker = pangenie_docker,
            monitoring_script = monitoring_script,
            cpu_make_count_model = cpu_make_count_model,
            runtime_attributes = leave_out_runtime_attributes,
            medium_runtime_attributes = leave_out_medium_runtime_attributes,
            large_runtime_attributes = leave_out_large_runtime_attributes,
            pangenie_runtime_attributes = pangenie_runtime_attributes,
            kage_runtime_attributes = kage_runtime_attributes,
            merge_runtime_attributes = merge_runtime_attributes,
            concat_runtime_attributes = concat_runtime_attributes,
            glimpse_phase_runtime_attributes = glimpse_phase_runtime_attributes,
            glimpse_sample_runtime_attributes = glimpse_sample_runtime_attributes,
            chromosome_to_glimpse_command_mem_gb = chromosome_to_glimpse_command_mem_gb,
            glimpse_phase_preemptible = glimpse_phase_preemptible,
            calculate_metrics_runtime_attributes = calculate_metrics_runtime_attributes,
            num_short_variants_to_retain = num_short_variants_to_retain,
            do_genotype_SVs = do_genotype_SVs
        }

#        # naively phased unphased GLIMPSE VCFs
#        call NaivelyPhase as GLIMPSENaivelyPhased { input:
#            vcf_gz = LeaveOutEvaluation.glimpse_unphased_vcf_gz,
#            vcf_gz_tbi = LeaveOutEvaluation.glimpse_unphased_vcf_gz_tbi,
#            docker = kage_docker,
#            monitoring_script = monitoring_script
#        }

        call FixVariantCollisions as GLIMPSEFixVariantCollisions { input:
            phased_bcf = LeaveOutEvaluation.glimpse_vcf_gz,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + ".glimpse.merged.phased.collisionless",
            monitoring_script = monitoring_script
        }

        if (do_pangenie) {
            # naively phase unphased PanGenie VCFs
            scatter (i in range(length(LeaveOutEvaluation.pangenie_vcf_gzs))) {
                call NaivelyPhase as PanGenieNaivelyPhased { input:
                    vcf_gz = select_first([LeaveOutEvaluation.pangenie_vcf_gzs[i]]),
                    vcf_gz_tbi = select_first([LeaveOutEvaluation.pangenie_vcf_gz_tbis[i]]),
                    docker = kage_docker,
                    monitoring_script = monitoring_script
                }
            }

            # merge naively phased PanGenie VCFs
            call Ivcfmerge as PanGenieNaivelyPhasedMergeAcrossSamples { input:
                vcf_gzs = PanGenieNaivelyPhased.naively_phased_vcf_gz,
                vcf_gz_tbis = PanGenieNaivelyPhased.naively_phased_vcf_gz_tbi,
                sample_names = leave_out_sample_names,
                output_prefix = output_prefix + ".pangenie.merged",
                docker = ivcfmerge_docker,
                monitoring_script = monitoring_script
            }
        }
    }

    scatter (i in range(length(vcfdist_bed_file_array))) {
        File vcfdist_bed_file = vcfdist_bed_file_array[i]
        String? vcfdist_extra_args = vcfdist_extra_args_array[i]

        if (defined(hiphase_short_vcf_gz)) {
            # evaluate HiPhase short
            call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateHiPhaseShort { input:
                samples = vcfdist_samples,
                confident_regions_bed_files = confident_regions_bed_files,
                truth_vcf = vcfdist_truth_vcf,
                truth_vcf_idx = vcfdist_truth_vcf_idx,
                eval_vcf = select_first([hiphase_short_vcf_gz]),
                eval_vcf_idx = select_first([hiphase_short_vcf_gz_tbi]),
                region = evaluation_chromosomes_regions_arg,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                vcfdist_bed_file = vcfdist_bed_file,
                vcfdist_extra_args = vcfdist_extra_args,
                vcfdist_mem_gb = vcfdist_mem_gb,
                overlap_phase_tag = "PS",
                overlap_metrics_docker = overlap_metrics_docker
            }
            String hiphase_short_label = "HiPhaseShort"
        }

        # evaluate HiPhase SV
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateHiPhaseSV { input:
            samples = vcfdist_samples,
            confident_regions_bed_files = confident_regions_bed_files,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = hiphase_sv_vcf_gz,
            eval_vcf_idx = hiphase_sv_vcf_gz_tbi,
            region = evaluation_chromosomes_regions_arg,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            vcfdist_bed_file = vcfdist_bed_file,
            vcfdist_extra_args = vcfdist_extra_args,
            vcfdist_mem_gb = vcfdist_mem_gb,
            overlap_phase_tag = "PS",
            overlap_metrics_docker = overlap_metrics_docker
        }

        # evaluate filtered short + SV
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateFiltered { input:
            samples = vcfdist_samples,
            confident_regions_bed_files = confident_regions_bed_files,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = ConcatVcfsFilterAndConcatVcfs.vcf_gz,
            eval_vcf_idx = ConcatVcfsFilterAndConcatVcfs.vcf_gz_tbi,
            region = evaluation_chromosomes_regions_arg,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            vcfdist_bed_file = vcfdist_bed_file,
            vcfdist_extra_args = vcfdist_extra_args,
            vcfdist_mem_gb = vcfdist_mem_gb,
            overlap_phase_tag = "PS",
            overlap_metrics_docker = overlap_metrics_docker
        }

#        # evaluate before Shapeit4 collisionless short + SV
#        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateBeforeShapeit4FixVariantCollisions { input:
#            samples = vcfdist_samples,
#            confident_regions_bed_files = confident_regions_bed_files,
#            truth_vcf = vcfdist_truth_vcf,
#            truth_vcf_idx = vcfdist_truth_vcf_idx,
#            eval_vcf = ConcatVcfsBeforeShapeit4FixVariantCollisions.vcf_gz,
#            eval_vcf_idx = ConcatVcfsBeforeShapeit4FixVariantCollisions.vcf_gz_tbi,
#            region = evaluation_chromosomes_regions_arg,
#            reference_fasta = reference_fasta,
#            reference_fasta_fai = reference_fasta_fai,
#            vcfdist_bed_file = vcfdist_bed_file,
#            vcfdist_extra_args = vcfdist_extra_args,
#            vcfdist_mem_gb = vcfdist_mem_gb,
#            overlap_phase_tag = "PS",
#            overlap_metrics_docker = overlap_metrics_docker
#        }

        # evaluate Shapeit4 short + SV
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateShapeit4 { input:
            samples = vcfdist_samples,
            confident_regions_bed_files = confident_regions_bed_files,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = ConcatVcfsShapeit4.vcf_gz,
            eval_vcf_idx = ConcatVcfsShapeit4.vcf_gz_tbi,
            region = evaluation_chromosomes_regions_arg,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            vcfdist_bed_file = vcfdist_bed_file,
            vcfdist_extra_args = vcfdist_extra_args,
            vcfdist_mem_gb = vcfdist_mem_gb,
            overlap_phase_tag = "NONE",
            overlap_metrics_docker = overlap_metrics_docker
        }

#        # evaluate collisionless short + SV
#        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateFixVariantCollisions { input:
#            samples = vcfdist_samples,
#            confident_regions_bed_files = confident_regions_bed_files,
#            truth_vcf = vcfdist_truth_vcf,
#            truth_vcf_idx = vcfdist_truth_vcf_idx,
#            eval_vcf = ConcatVcfsFixVariantCollisions.vcf_gz,
#            eval_vcf_idx = ConcatVcfsFixVariantCollisions.vcf_gz_tbi,
#            region = evaluation_chromosomes_regions_arg,
#            reference_fasta = reference_fasta,
#            reference_fasta_fai = reference_fasta_fai,
#            vcfdist_bed_file = vcfdist_bed_file,
#            vcfdist_extra_args = vcfdist_extra_args,
#            vcfdist_mem_gb = vcfdist_mem_gb,
#            overlap_phase_tag = "NONE",
#            overlap_metrics_docker = overlap_metrics_docker
#        }

        # evaluate panel short + SV
        call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluatePanel { input:
            samples = vcfdist_samples,
            confident_regions_bed_files = confident_regions_bed_files,
            truth_vcf = vcfdist_truth_vcf,
            truth_vcf_idx = vcfdist_truth_vcf_idx,
            eval_vcf = ConcatVcfsPanGeniePanelCreation.vcf_gz,
            eval_vcf_idx = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi,
            region = evaluation_chromosomes_regions_arg,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            vcfdist_bed_file = vcfdist_bed_file,
            vcfdist_extra_args = vcfdist_extra_args,
            vcfdist_mem_gb = vcfdist_mem_gb,
            overlap_phase_tag = "NONE",
            overlap_metrics_docker = overlap_metrics_docker
        }

        if (do_short_read_stages) {
#            # evaluate naively phased GLIMPSE
#            call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateGLIMPSENaivelyPhased { input:
#                samples = vcfdist_samples,
#                confident_regions_bed_files = confident_regions_bed_files,
#                truth_vcf = vcfdist_truth_vcf,
#                truth_vcf_idx = vcfdist_truth_vcf_idx,
#                eval_vcf = select_first([GLIMPSENaivelyPhased.naively_phased_vcf_gz]),
#                eval_vcf_idx = select_first([GLIMPSENaivelyPhased.naively_phased_vcf_gz_tbi]),
#                region = evaluation_chromosomes_regions_arg,
#                reference_fasta = reference_fasta,
#                reference_fasta_fai = reference_fasta_fai,
#                vcfdist_bed_file = vcfdist_bed_file,
#                vcfdist_extra_args = vcfdist_extra_args,
#                vcfdist_mem_gb = vcfdist_mem_gb,
#                overlap_phase_tag = "NONE",
#                overlap_metrics_docker = overlap_metrics_docker
#            }
#            String glimpse_naively_phased_label = "GLIMPSENaivelyPhased"

#            # evaluate GLIMPSE
#            call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateGLIMPSE { input:
#                samples = vcfdist_samples,
#                confident_regions_bed_files = confident_regions_bed_files,
#                truth_vcf = vcfdist_truth_vcf,
#                truth_vcf_idx = vcfdist_truth_vcf_idx,
#                eval_vcf = select_first([LeaveOutEvaluation.glimpse_vcf_gz]),
#                eval_vcf_idx = select_first([LeaveOutEvaluation.glimpse_vcf_gz_tbi]),
#                region = evaluation_chromosomes_regions_arg,
#                reference_fasta = reference_fasta,
#                reference_fasta_fai = reference_fasta_fai,
#                vcfdist_bed_file = vcfdist_bed_file,
#                vcfdist_extra_args = vcfdist_extra_args,
#                vcfdist_mem_gb = vcfdist_mem_gb,
#                overlap_phase_tag = "NONE",
#                overlap_metrics_docker = overlap_metrics_docker
#            }
#            String glimpse_label = "GLIMPSE"

            # evaluate collisionless GLIMPSE
            call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluateGLIMPSEFixVariantCollisions { input:
                samples = vcfdist_samples,
                confident_regions_bed_files = confident_regions_bed_files,
                truth_vcf = vcfdist_truth_vcf,
                truth_vcf_idx = vcfdist_truth_vcf_idx,
                eval_vcf = select_first([GLIMPSEFixVariantCollisions.collisionless_bcf]),
                eval_vcf_idx = select_first([GLIMPSEFixVariantCollisions.collisionless_bcf_csi]),
                region = evaluation_chromosomes_regions_arg,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                vcfdist_bed_file = vcfdist_bed_file,
                vcfdist_extra_args = vcfdist_extra_args,
                vcfdist_mem_gb = vcfdist_mem_gb,
                overlap_phase_tag = "NONE",
                overlap_metrics_docker = overlap_metrics_docker
            }
            String glimpse_fix_variant_collisions_label = "GLIMPSEFixVariantCollisions"

            if (do_pangenie) {
                # evaluate naively phased PanGenie
                call VcfdistAndOverlapMetricsEvaluation.VcfdistAndOverlapMetricsEvaluation as EvaluatePanGenieNaivelyPhased { input:
                    samples = vcfdist_samples,
                    confident_regions_bed_files = confident_regions_bed_files,
                    truth_vcf = vcfdist_truth_vcf,
                    truth_vcf_idx = vcfdist_truth_vcf_idx,
                    eval_vcf = select_first([PanGenieNaivelyPhasedMergeAcrossSamples.merged_vcf_gz]),
                    eval_vcf_idx = select_first([PanGenieNaivelyPhasedMergeAcrossSamples.merged_vcf_gz_tbi]),
                    region = evaluation_chromosomes_regions_arg,
                    reference_fasta = reference_fasta,
                    reference_fasta_fai = reference_fasta_fai,
                    vcfdist_bed_file = vcfdist_bed_file,
                    vcfdist_extra_args = vcfdist_extra_args,
                    vcfdist_mem_gb = vcfdist_mem_gb,
                    overlap_phase_tag = "NONE",
                    overlap_metrics_docker = overlap_metrics_docker
                }
                String pangenie_naively_phased_label = "PanGenieNaivelyPhased"
            }
        }

        Array[String] labels_per_vcf = select_all([
            hiphase_short_label,
            "HiPhaseSV",
            "ConcatAndFiltered",
#            "BeforeShapeit4FixVariantCollisions",
            "Shapeit4",
#            "FixVariantCollisions",
            "Panel",
#            glimpse_naively_phased_label,
#            glimpse_label,
            glimpse_fix_variant_collisions_label,
            pangenie_naively_phased_label
        ])
        call SummarizeEvaluations { input:
            labels_per_vcf = labels_per_vcf,
            vcfdist_outputs_per_vcf_and_sample = select_all([
                EvaluateHiPhaseShort.precision_recall_summary_per_sample,
                EvaluateHiPhaseSV.precision_recall_summary_per_sample,
                EvaluateFiltered.precision_recall_summary_per_sample,
#                EvaluateBeforeShapeit4FixVariantCollisions.precision_recall_summary_per_sample,
                EvaluateShapeit4.precision_recall_summary_per_sample,
#                EvaluateFixVariantCollisions.precision_recall_summary_per_sample,
                EvaluatePanel.precision_recall_summary_per_sample,
#                EvaluateGLIMPSENaivelyPhased.precision_recall_summary_per_sample,
#                EvaluateGLIMPSE.precision_recall_summary_per_sample,
                EvaluateGLIMPSEFixVariantCollisions.precision_recall_summary_per_sample,
                EvaluatePanGenieNaivelyPhased.precision_recall_summary_per_sample
            ]),
            overlap_metrics_outputs_per_vcf = select_all([
                EvaluateHiPhaseShort.overlap_metrics_outputs,
                EvaluateHiPhaseSV.overlap_metrics_outputs,
                EvaluateFiltered.overlap_metrics_outputs,
#                EvaluateBeforeShapeit4FixVariantCollisions.overlap_metrics_outputs,
                EvaluateShapeit4.overlap_metrics_outputs,
#                EvaluateFixVariantCollisions.overlap_metrics_outputs,
                EvaluatePanel.overlap_metrics_outputs,
#                EvaluateGLIMPSENaivelyPhased.overlap_metrics_outputs,
#                EvaluateGLIMPSE.overlap_metrics_outputs,
                EvaluateGLIMPSEFixVariantCollisions.overlap_metrics_outputs,
                EvaluatePanGenieNaivelyPhased.overlap_metrics_outputs
            ]),
            docker = summarize_evaluations_docker
        }
    }

    output {
        File panel_vcf_gz = ConcatVcfsPanGeniePanelCreation.vcf_gz
        File panel_vcf_gz_tbi = ConcatVcfsPanGeniePanelCreation.vcf_gz_tbi

        Array[File] evaluation_summary_tsvs = SummarizeEvaluations.evaluation_summary_tsv
    }
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix
        String? extra_concat_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    Int disk_size_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools concat ~{sep=' ' vcf_gzs} ~{extra_concat_args} -Oz -o ~{output_prefix}.vcf.gz
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

task FixVariantCollisions {

    input {
        File phased_bcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix

        File? monitoring_script
    }

    command <<<
        set -euxo pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # convert bcf to vcf.gz
        bcftools view ~{phased_bcf} -Oz -o phased.vcf.gz
        bcftools index -t phased.vcf.gz

        java ~{fix_variant_collisions_java} \
            phased.vcf.gz \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools view collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view -Oz -o ~{output_prefix}.vcf.gz
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.vcf.gz
        bcftools view ~{output_prefix}.vcf.gz -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File collisionless_bcf = "~{output_prefix}.bcf"
        File collisionless_bcf_csi = "~{output_prefix}.bcf.csi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}

task SummarizeEvaluations {
    input {
        Array[String] labels_per_vcf
#        Array[Array[VcfdistOutputs]] vcfdist_outputs_per_vcf_and_sample
        Array[Array[VcfdistPrecisionRecallSummary]] vcfdist_outputs_per_vcf_and_sample
        Array[OverlapMetricsOutputs] overlap_metrics_outputs_per_vcf

        String docker
        Int disk_size_gb = 50
        Int boot_disk_size_gb = 15
        Int mem_gb = 8
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        python - --labels_per_vcf_txt ~{write_lines(labels_per_vcf)} \
                 --vcfdist_outputs_per_vcf_and_sample_json ~{write_json(vcfdist_outputs_per_vcf_and_sample)} \
                 --overlap_metrics_outputs_per_vcf_json ~{write_json(overlap_metrics_outputs_per_vcf)} \
                 <<-'EOF'
        import argparse
        import json
        import pandas as pd

        def summarize(labels_per_vcf_txt,
                      vcfdist_outputs_per_vcf_and_sample_json,
                      overlap_metrics_outputs_per_vcf_json):
            with open(labels_per_vcf_txt) as f:
                labels = f.read().splitlines()

            with open(vcfdist_outputs_per_vcf_and_sample_json) as f:
                vcfdist_outputs_per_vcf_and_sample = json.load(f)

            with open(overlap_metrics_outputs_per_vcf_json) as f:
                overlap_metrics_outputs_per_vcf = json.load(f)

            summary_dict = {}
            for i, label in enumerate(labels):
                summary_dict[label] = {}
                summary_dict[label]['NUM_VCFDIST_SAMPLES'] = len(vcfdist_outputs_per_vcf_and_sample[i])
                summary_dict[label].update(summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_vcf_and_sample[i]))
                if overlap_metrics_outputs_per_vcf:
                    summary_dict[label].update(summarize_overlap_metrics_outputs(overlap_metrics_outputs_per_vcf[i]))

            pd.DataFrame.from_dict(summary_dict, orient='index').to_csv('evaluation_summary.tsv', sep='\t', float_format="%.4f")

        def summarize_vcfdist_outputs_over_samples(vcfdist_outputs_per_sample):
            precision_recall_metrics_dict = {}
            for s, vcfdist_outputs in enumerate(vcfdist_outputs_per_sample):
                precision_recall_metrics_dict[s] = {}
                pr_metrics_df = pd.read_csv(vcfdist_outputs['precision_recall_summary_tsv'], sep='\t', index_col=[0, 1])
                for var_type in ['SNP', 'INDEL', 'SV']:
                    var_type_metrics_dict = pr_metrics_df.loc[var_type, 'NONE'][['TRUTH_TP', 'QUERY_TP', 'TRUTH_FN', 'QUERY_FP', 'PREC', 'RECALL', 'F1_SCORE']].add_prefix(f'{var_type}_').to_dict()
                    precision_recall_metrics_dict[s].update(var_type_metrics_dict)
            return pd.DataFrame.from_dict(precision_recall_metrics_dict, orient='index').mean(axis=0)


        def summarize_overlap_metrics_outputs(overlap_metrics_outputs):
            return pd.read_csv(overlap_metrics_outputs['metrics_tsv'], sep='\t').iloc[0].to_dict()

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--labels_per_vcf_txt',
                                type=str)

            parser.add_argument('--vcfdist_outputs_per_vcf_and_sample_json',
                                type=str)

            parser.add_argument('--overlap_metrics_outputs_per_vcf_json',
                                type=str)

            args = parser.parse_args()

            summarize(args.labels_per_vcf_txt,
                      args.vcfdist_outputs_per_vcf_and_sample_json,
                      args.overlap_metrics_outputs_per_vcf_json)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        File evaluation_summary_tsv = "evaluation_summary.tsv"
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        bootDiskSizeGb: boot_disk_size_gb
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

# assumes all VCFs have identical variants
# TODO make this work for identical filenames
task Ivcfmerge {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        Array[String] sample_names
        String output_prefix
        String? region_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        wget https://github.com/iqbal-lab-org/ivcfmerge/archive/refs/tags/v1.0.0.tar.gz
        tar -xvf v1.0.0.tar.gz

        mkdir compressed
        mv ~{sep=' ' vcf_gzs} compressed
        mv ~{sep=' ' vcf_gz_tbis} compressed

        if [ $(ls compressed/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls compressed/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls compressed/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            mkdir decompressed
            ls compressed/*.vcf.gz | xargs -I % sh -c 'bcftools annotate --no-version ~{region_args} -x INFO % --threads 2 -Ov -o decompressed/$(basename % .gz)'
            time python ivcfmerge-1.0.0/ivcfmerge.py <(ls decompressed/*.vcf) ~{output_prefix}.vcf
            bcftools annotate --no-version -S ~{write_lines(sample_names)} -x FORMAT/FT ~{output_prefix}.vcf --threads 2 -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
    }

    output {
        File monitoring_log = "monitoring.log"
        File merged_vcf_gz = "~{output_prefix}.vcf.gz"
        File merged_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 4])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 250]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task NaivelyPhase {
    input {
        File vcf_gz
        File vcf_gz_tbi

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String output_prefix = basename(vcf_gz, ".vcf.gz")

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # naively set all GTs to phased for Vcfdist evaluation
        bcftools +setGT ~{vcf_gz} -Oz -o ~{output_prefix}.naively_phased.vcf.gz -- -t a -n p
        bcftools index -t ~{output_prefix}.naively_phased.vcf.gz
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
        File naively_phased_vcf_gz = "~{output_prefix}.naively_phased.vcf.gz"
        File naively_phased_vcf_gz_tbi = "~{output_prefix}.naively_phased.vcf.gz.tbi"
    }
}
