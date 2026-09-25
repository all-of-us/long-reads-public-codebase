########################################################################################
# FilterIntrasampleDevV11
# Phase 1 stage b (per-sample filtering), step 3 of the callset-integration launch order.
#
# Scores each SV in one sample with an isolation-forest model (GATK
# ExtractVariantAnnotations / TrainVariantAnnotationsModel / ScoreVariantAnnotations).
# Training sites are the sample's SVs that match HPRC assembly-based (dipcall) truth
# calls. The output VCF carries a CALIBRATION_SENSITIVITY value per record; the Phase 1
# strict callset keeps records with CALIBRATION_SENSITIVITY <= 0.7 (applied in
# TruvariIntersample).
#
# The GATK filtering workflow is imported at launch time from a FireCloud method
# (see the import line below), so it is not vendored in repo/.
#
# Runs once per sample.
# Inputs:  truvari_collapsed_vcf_gz and regenotyped_vcf_gz (from the two previous steps),
#          training_resource_vcf_gz and training_resource_bed (HPRC truth), annotations
#          (the INFO/FORMAT fields used as model features), reference_fasta,
#          resolve_python_script (../resolve.py)
# Outputs: scored_vcf_gz  the scored VCF (stored in Terra as v11_scored_kanpig_vcf_gz)
# Previous: GetRegenotypedVcfKanpigV11. Next: FilterAndSplit.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/slee:FilterIntrasampleJointVcfFilteringDev/versions/2/plain-WDL/descriptor" as JointVcfFiltering

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

workflow FilterIntrasampleDevV11 {
    input {
        String sample_name
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_vcf_gz_tbi
        File regenotyped_vcf_gz
        File regenotyped_vcf_gz_tbi
        File training_resource_vcf_gz
        File training_resource_vcf_gz_tbi

        File? monitoring_script

        # PreprocessVCF/PostprocessVCF arguments
        File reference_fasta
        File reference_fasta_fai
        Boolean is_sniffles
        Boolean is_training_resource_split_and_resolved
        File resolve_python_script
        String preprocess_and_postprocess_docker

        # IdentifyTrainingSites arguments
        File training_resource_bed
        String truvari_docker

        # JointVcfFiltering arguments
        Array[String] annotations

        String? model_backend
        File? training_python_script
        File? hyperparameters_json
        File? scoring_python_script

        String? extract_extra_args = "--mode INDEL"
        String? train_extra_args = "--mode INDEL"
        String? score_extra_args = "--mode INDEL"

        String gatk_docker
        File? gatk_override

        RuntimeAttributes? extract_runtime_attributes
        RuntimeAttributes? train_runtime_attributes
        RuntimeAttributes? score_runtime_attributes
    }

    call PreprocessVCF {
        input:
            sample_name = sample_name,
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_vcf_gz_tbi = truvari_collapsed_vcf_gz_tbi,
            regenotyped_vcf_gz = regenotyped_vcf_gz,
            regenotyped_vcf_gz_tbi = regenotyped_vcf_gz_tbi,
            training_resource_vcf_gz = training_resource_vcf_gz,
            training_resource_vcf_gz_tbi = training_resource_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            resolve_python_script = resolve_python_script,
            is_sniffles = is_sniffles,
            is_training_resource_split_and_resolved = is_training_resource_split_and_resolved,
            output_prefix = sample_name,
            docker = preprocess_and_postprocess_docker,
            monitoring_script = monitoring_script
    }

    call IdentifyTrainingSites {
        input:
            preprocessed_vcf_gz = PreprocessVCF.preprocessed_vcf_gz,
            preprocessed_vcf_gz_tbi = PreprocessVCF.preprocessed_vcf_gz_tbi,
            preprocessed_training_resource_vcf_gz = PreprocessVCF.preprocessed_training_resource_vcf_gz,
            preprocessed_training_resource_vcf_gz_tbi = PreprocessVCF.preprocessed_training_resource_vcf_gz_tbi,
            training_resource_bed = training_resource_bed,
            docker = truvari_docker,
            monitoring_script = monitoring_script
    }

    call JointVcfFiltering.JointVcfFiltering {
        input:
            input_vcfs = [PreprocessVCF.preprocessed_vcf_gz],
            input_vcf_idxs = [PreprocessVCF.preprocessed_vcf_gz_tbi],
            sites_only_vcf = PreprocessVCF.preprocessed_vcf_gz,
            sites_only_vcf_idx = PreprocessVCF.preprocessed_vcf_gz_tbi,
            output_prefix = sample_name,
            annotations = annotations,
            resource_vcf = IdentifyTrainingSites.training_vcf_gz,
            resource_vcf_idx = IdentifyTrainingSites.training_vcf_gz_tbi,
            model_backend = model_backend,
            training_python_script = training_python_script,
            hyperparameters_json = hyperparameters_json,
            scoring_python_script = scoring_python_script,
            extract_extra_args = extract_extra_args,
            train_extra_args = train_extra_args,
            score_extra_args = score_extra_args,
            gatk_docker = gatk_docker,
            gatk_override = gatk_override,
            extract_runtime_attributes = extract_runtime_attributes,
            train_runtime_attributes = train_runtime_attributes,
            score_runtime_attributes = score_runtime_attributes,
            monitoring_script = monitoring_script
    }

    call PostprocessVCF {
        input:
            scored_vcf_gz = JointVcfFiltering.scored_vcfs[0],
            scored_vcf_gz_tbi = JointVcfFiltering.scored_vcf_idxs[0],
            output_prefix = sample_name,
            docker = preprocess_and_postprocess_docker,
            monitoring_script = monitoring_script
    }

    output {
        File extracted_annotations_hdf5 = JointVcfFiltering.extracted_annotations_hdf5
        File? extracted_unlabeled_annotations_hdf5 = JointVcfFiltering.extracted_unlabeled_annotations_hdf5
        File extracted_vcf_gz = JointVcfFiltering.extracted_vcf
        File extracted_vcf_gz_tbi = JointVcfFiltering.extracted_vcf_idx

        Array[File] model_files = JointVcfFiltering.model_files

        File? annotations_hdf5 = JointVcfFiltering.annotations_hdf5s[0]
        File? scores_hdf5 = JointVcfFiltering.scores_hdf5s[0]

        File scored_vcf_gz = PostprocessVCF.postprocessed_vcf_gz
        File scored_vcf_gz_tbi = PostprocessVCF.postprocessed_vcf_gz_tbi

        Array[File?] monitoring_logs = flatten(
          [
            [PreprocessVCF.monitoring_log],
            [IdentifyTrainingSites.monitoring_log],
            JointVcfFiltering.monitoring_logs,
            [PostprocessVCF.monitoring_log],
          ])
    }
}

# takes in truvari-collapsed and resolved+regenotyped VCFs and populates the following INFO fields from the corresponding FORMAT fields:
# Sniffles: ['GQ', 'DR', 'DV', 'SUPP_PAV', 'SUPP_SNIFFLES', 'SUPP_PBSV']
# kanpig: ['KS_1', 'KS_2', 'SQ', 'GQ', 'DP', 'AD_NON_ALT', 'AD_ALL', 'SUPP_PAV', 'SUPP_SNIFFLES', 'SUPP_PBSV']
# TODO: make this production-level code
task PreprocessVCF {
    input {
        String sample_name
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_vcf_gz_tbi
        File regenotyped_vcf_gz
        File regenotyped_vcf_gz_tbi
        File training_resource_vcf_gz
        File training_resource_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        File resolve_python_script
        Boolean is_sniffles
        Boolean is_training_resource_split_and_resolved
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String truvari_collapsed_prefix = basename(truvari_collapsed_vcf_gz, ".vcf.gz")
    String training_resource_prefix = basename(training_resource_vcf_gz, ".vcf.gz")

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # TODO clean this up
        # split multiallelics
        bcftools norm -m- ~{truvari_collapsed_vcf_gz} -Oz -o ~{truvari_collapsed_prefix}.split.vcf.gz

        # resolve symbolic alleles
        python ~{resolve_python_script} ~{truvari_collapsed_prefix}.split.vcf.gz ~{reference_fasta} | \
            bcftools view -i "SVTYPE!='BND'" | \
            bcftools sort -Oz -o ~{truvari_collapsed_prefix}.split.resolved.vcf.gz
        bcftools index -t ~{truvari_collapsed_prefix}.split.resolved.vcf.gz
        
        if ~{is_training_resource_split_and_resolved}; then
            cp ~{training_resource_vcf_gz} ~{training_resource_prefix}.preprocessed.vcf.gz
        else
            bcftools norm -m- ~{training_resource_vcf_gz} -Oz -o ~{training_resource_prefix}.split.vcf.gz
            python ~{resolve_python_script} ~{training_resource_prefix}.split.vcf.gz ~{reference_fasta} | \
                bcftools view -i "SVTYPE!='BND'" | \
                bcftools reheader -s <(echo {sample_name}) | \
                bcftools sort -Oz -o ~{training_resource_prefix}.preprocessed.vcf.gz
        fi
        bcftools index -t ~{training_resource_prefix}.preprocessed.vcf.gz

        # create supplemental header for bcftools annotate
        touch format.hdr.txt
        if ~{is_sniffles}; then
            # GQ, DR, DV from Sniffles-regenotyped VCF
            for format_annot in GQ DR DV
            do
                bcftools view -h ~{regenotyped_vcf_gz} | grep ID=$format_annot, | sed -e 's/FORMAT/INFO/g' >> format.hdr.txt
            done
        else
            # KS, SQ, GQ, DP, AD, GT from kanpig-regenotyped VCF
            for format_annot in SQ GQ DP
            do
                bcftools view -h ~{regenotyped_vcf_gz} | grep ID=$format_annot, | sed -e 's/FORMAT/INFO/g' >> format.hdr.txt
            done
            echo '##INFO=<ID=AD_NON_ALT,Number=1,Type=Integer,Description="Coverage for non-alternate alleles">' >> format.hdr.txt
            echo '##INFO=<ID=AD_ALL,Number=1,Type=Integer,Description="Coverage for all alleles">' >> format.hdr.txt
            echo '##INFO=<ID=KS_1,Number=1,Type=Integer,Description="Kanpig score 1">' >> format.hdr.txt
            echo '##INFO=<ID=KS_2,Number=1,Type=Integer,Description="Kanpig score 2">' >> format.hdr.txt
            echo '##INFO=<ID=GT_COUNT,Number=1,Type=Integer,Description="GT converted into an int in {0,1,2}.">' >> format.hdr.txt
        fi
        echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' >> format.hdr.txt
        echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> format.hdr.txt
        echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> format.hdr.txt

        # convert SUPP integer from resolved VCF to SUPP_* binary indicators
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%SUPP]\n' ~{truvari_collapsed_prefix}.split.resolved.vcf.gz > format.supp.tsv
        python <<EOF
import numpy as np
import pandas as pd

print('Converting SUPP to SUPP_* for each caller...')
format_annot_supp_binary_tsv = 'format.supp_binary.tsv'
format_annot_supp_df = pd.read_csv('format.supp.tsv' , sep='\t', header=None,
                                   names=['CHROM', 'POS', 'REF', 'ALT', 'SUPP'],
                                   dtype={'CHROM': str, 'POS': np.int64, 'REF': str, 'ALT': str, 'SUPP': np.uint8})
format_annot_supp_df['SUPP_BITCOUNT'] = np.unpackbits(format_annot_supp_df['SUPP'].values[:, np.newaxis], axis=1).sum(axis=1)
format_annot_supp_df[['SUPP_PAV', 'SUPP_SNIFFLES', 'SUPP_PBSV']] = format_annot_supp_df['SUPP'].map(lambda i: np.array(list(np.binary_repr(i).zfill(3))).astype(np.uint8)).tolist()
np.savetxt('format.supp_binary.tsv', format_annot_supp_df[['SUPP_PAV', 'SUPP_SNIFFLES', 'SUPP_PBSV']], fmt='%d')

EOF

        # replace SUPP column in format.supp.tsv with SUPP_PAV, SUPP_SNIFFLES, SUPP_PBSV columns
        awk -v OFS='\t' 'FNR==NR{a[NR]=$1;b[NR]=$2;c[NR]=$3;next}{$5=a[FNR];$6=b[FNR];$7=c[FNR]}1' format.supp_binary.tsv format.supp.tsv | bgzip -c > format.supp_binary.tsv.gz
        tabix -s1 -b2 -e2 format.supp_binary.tsv.gz
        
        echo 'Annotating SUPP_* for each caller...'
        bcftools annotate -a format.supp_binary.tsv.gz -h format.hdr.txt -c CHROM,POS,REF,ALT,SUPP_PAV,SUPP_SNIFFLES,SUPP_PBSV ~{regenotyped_vcf_gz} -Oz -o format.supp_binary.vcf.gz
        bcftools index -t format.supp_binary.vcf.gz

        # FC> Ensuring that every record has a unique ID. We need to join by
        # CHROM,POS,ID in what follows, since using CHROM,POS,REF,ALT makes
        # bcftools annotate segfault. 
        bcftools view --header-only format.supp_binary.vcf.gz > tmp.vcf
        bcftools view --no-header format.supp_binary.vcf.gz | awk 'BEGIN { i=0; } { printf("%s\t%s\t%d-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,++i,$3,$4,$5,$6,$7,$8,$9,$10); }' >> tmp.vcf
        bgzip tmp.vcf
        tabix -f tmp.vcf.gz
        bcftools view --no-header tmp.vcf.gz | head -n 10 || echo "0"
        
        echo 'Annotating non-SUPP_* annotations...'
        if ~{is_sniffles}; then
            # GQ, DR, DV from Sniffles-regenotyped VCF
            bcftools query -f '%CHROM\t%POS\t%ID\t[%GQ]\t[%DR]\t[%DV]\n' tmp.vcf.gz | bgzip -c > format.tsv.gz
            tabix -s1 -b2 -e2 format.tsv.gz
            bcftools annotate --threads ${N_THREADS} -a format.tsv.gz -h format.hdr.txt -c CHROM,POS,ID,GQ,DR,DV tmp.vcf.gz -Oz -o ~{output_prefix}.preprocessed.vcf.gz
        else
            # KS, SQ, GQ, DP, AD, GT from kanpig-regenotyped VCF; TODO nicer AD
            bcftools query -f '%CHROM\t%POS\t%ID\t[%KS]\t[%SQ]\t[%GQ]\t[%DP]\t[%AD]\t[%GT]\n' tmp.vcf.gz | awk '{ \
                KS_1=0; KS_2=0; \
                p=0; \
                for (i=1; i<=length($4); i++) { \
                    if (substr($4,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { KS_1=$4; KS_2=$4; } \
                else { KS_1=substr($4,1,p-1); KS_2=substr($4,p+1); } \
                if (KS_1==".") KS_1=0; \
                if (KS_2==".") KS_2=0; \
                \
                SQ=$5; \
                if (SQ==".") SQ=0; \
                \
                GQ=$6; \
                if (GQ==".") GQ=0; \
                \
                DP=$7; \
                if (DP==".") DP=0; \
                \
                AD_NON_ALT=0; AD_ALL=0; \
                p=0; \
                for (i=1; i<=length($8); i++) { \
                    if (substr($8,i,1)==",") { p=i; break; } \
                } \
                if (p==0) { AD_NON_ALT=$8; AD_ALL=$8; } \
                else { AD_NON_ALT=substr($8,1,p-1); AD_ALL=substr($8,p+1); } \
                if (AD_NON_ALT==".") AD_NON_ALT=0; \
                if (AD_ALL==".") AD_ALL=0; \
                \
                GT_COUNT=0; \
                if ($9=="0/0" || $9=="0|0" || $9=="./."  || $9==".|." || $9=="./0" || $9==".|0" || $9=="0/." || $9=="0|.") GT_COUNT=0; \
                else if ($9=="0/1" || $9=="0|1" || $9=="1/0" || $9=="1|0" || $9=="./1" || $9==".|1" || $9=="1/." || $9=="1|.") GT_COUNT=1; \
                else if ($9=="1/1" || $9=="1|1") GT_COUNT=2; \
                \
                printf("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",$1,$2,$3,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT); \
            }' | bgzip -c > format.tsv.gz
            tabix -s1 -b2 -e2 format.tsv.gz
            bcftools annotate --threads ${N_THREADS} -a format.tsv.gz -h format.hdr.txt -c CHROM,POS,ID,KS_1,KS_2,SQ,GQ,DP,AD_NON_ALT,AD_ALL,GT_COUNT tmp.vcf.gz -Oz -o ~{output_prefix}.preprocessed.vcf.gz
        fi
        bcftools view --no-header ~{output_prefix}.preprocessed.vcf.gz | head -n 10 || echo "0"
        bcftools index -t ~{output_prefix}.preprocessed.vcf.gz

        # TODO do we still need to hard-filter SVLEN >= 50 for extract, and should we clear existing filters?
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
        File preprocessed_vcf_gz = "~{output_prefix}.preprocessed.vcf.gz"
        File preprocessed_vcf_gz_tbi = "~{output_prefix}.preprocessed.vcf.gz.tbi"
        File preprocessed_training_resource_vcf_gz = "~{training_resource_prefix}.preprocessed.vcf.gz"
        File preprocessed_training_resource_vcf_gz_tbi = "~{training_resource_prefix}.preprocessed.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}

task IdentifyTrainingSites {
    input {
        File preprocessed_vcf_gz
        File preprocessed_vcf_gz_tbi
        File preprocessed_training_resource_vcf_gz
        File preprocessed_training_resource_vcf_gz_tbi
        File training_resource_bed
        String truvari_extra_args = "--sizemin 50 --sizemax 50 --sizefilt 50"

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        truvari bench \
            -b ~{preprocessed_training_resource_vcf_gz} \
            -c ~{preprocessed_vcf_gz} \
            --includebed ~{training_resource_bed} \
            ~{truvari_extra_args} \
            -o truvari/
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
        File training_vcf_gz = "truvari/tp-comp.vcf.gz"
        File training_vcf_gz_tbi = "truvari/tp-comp.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}

# copy SCORE from INFO field to FORMAT field
task PostprocessVCF {
    input {
        File scored_vcf_gz
        File scored_vcf_gz_tbi
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -eou pipefail

        if [ -s ~{monitoring_script} ]; then
          bash ~{monitoring_script} > monitoring.log &
        fi

        # create supplemental header for bcftools annotate
        echo '##FORMAT=<ID=SCORE,Number=1,Type=Float,Description="Score according to the model applied by ScoreVariantAnnotations">' > format.hdr.txt

        # copy SCORE INFO to FORMAT
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SCORE\n' ~{scored_vcf_gz} | bgzip -c > format.score.tsv.gz
        tabix -s1 -b2 -e2 format.score.tsv.gz
        bcftools annotate -a format.score.tsv.gz -h format.hdr.txt -c CHROM,POS,REF,ALT,FORMAT/SCORE ~{scored_vcf_gz} -Oz -o ~{output_prefix}.score.vcf.gz
        bcftools index -t ~{output_prefix}.score.vcf.gz
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
        File postprocessed_vcf_gz = "~{output_prefix}.score.vcf.gz"
        File postprocessed_vcf_gz_tbi = "~{output_prefix}.score.vcf.gz.tbi"
        File? monitoring_log = "monitoring.log"
    }
}