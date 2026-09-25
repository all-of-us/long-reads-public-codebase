########################################################################################
# GraphEvaluation
# Evaluation, not on the Phase 1 launch path.
#
# Compares SV callsets by building haplotype graphs in windows around SVs and aligning
# assembled haplotypes to them (GraphAligner), per chromosome. Windows can be
# restricted by overlap with evaluation BEDs.
#
# Inputs:  vcf_gz / vcf_id (callsets to compare), truth_vcf_id, haps_vs_chm13_csv
#          (haplotype-to-reference BAMs), tandems_bed, reference_fa, window settings
# Outputs: evaluation, analysis_small_overlap, analysis_large_overlap  per chromosome
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


#
#
workflow GraphEvaluation {
    input {
        Array[String] vcf_id
        String? truth_vcf_id
        Array[File] vcf_gz
        Array[File] vcf_tbi
        Int interval_max_length
        Int flank_length
        Int graphaligner_timeout = 120
        File tandems_bed
        File reference_fa
        File haps_vs_chm13_csv
        File? confident_bed
        Float? small_overlap
        Array[File]? evaluation_bed_small_overlap
        Float? large_overlap
        Array[File]? evaluation_bed_large_overlap
        Boolean force_unique_reads
        Int min_sv_length = 10
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        File? force_windows_bed
        String docker
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
        vcf_id: "Distinct name assigned to each VCF to be evaluated"
        truth_vcf_id: "One of the elements of `vcf_id`, assumed to contain the true set of SVs."
        vcf_gz: "The VCF files to be evaluated. Haplotype clusters are built using the first file in this list."
        interval_max_length: "Do not use windows longer than this for evaluation"
        flank_length: "Ensure this amount of non-tandem sequence around each window"
        tandems_bed: "Tandem track"
        haps_vs_chm13_csv: "List of remote haplotype-reference BAMs"
        evaluation_bed_small_overlap: "For every BED file in this list: use only windows with at least `small_overlap` fraction of bases covered by intervals in the BED."
        evaluation_bed_large_overlap: "For every BED file in this list: use only windows with at least `large_overlap` fraction of bases covered by intervals in the BED."
        chromosomes: "Use only these chromosomes. Each chromosome is processed in parallel and produces separate evaluation files."
        force_unique_reads: "Intended for resolving collisions in identically named sequences across samples. If true, then force sequence names to have a suffix _[sample_name] where sample_name is the unique name provided in the BAM CSV for each BAM"
        graphaligner_timeout: "Maximum time in seconds to allow for graph alignment step"
    }

    scatter(chromosome in chromosomes) {
        call EvaluateChromosome {
            input:
                vcf_id = vcf_id,
                truth_vcf_id = truth_vcf_id,
                vcf_gz = vcf_gz,
                vcf_tbi = vcf_tbi,
                interval_max_length = interval_max_length,
                flank_length = flank_length,
                graphaligner_timeout = graphaligner_timeout,
                tandems_bed = tandems_bed,
                reference_fa = reference_fa,
                haps_vs_chm13_csv = haps_vs_chm13_csv,
                confident_bed = confident_bed,
                small_overlap = small_overlap,
                evaluation_bed_small_overlap = evaluation_bed_small_overlap,
                large_overlap = large_overlap,
                evaluation_bed_large_overlap = evaluation_bed_large_overlap,
                min_sv_length = min_sv_length,
                chromosome = chromosome,
                force_unique_reads = force_unique_reads,
                docker = docker,
                force_windows_bed = force_windows_bed,
                n_cpu = n_cpu,
                ram_size_gb = ram_size_gb
        }
    }

    output {
        Array[File] evaluation = EvaluateChromosome.evaluation
        Array[File] analysis_small_overlap = EvaluateChromosome.analysis_small_overlap
        Array[File] analysis_large_overlap = EvaluateChromosome.analysis_large_overlap
        Array[File] monitor_log = EvaluateChromosome.monitor_log
    }
}


task EvaluateChromosome {
    input {
        Array[String] vcf_id
        String? truth_vcf_id
        Array[File] vcf_gz
        Array[File] vcf_tbi
        Int interval_max_length
        Int flank_length
        Int graphaligner_timeout
        File tandems_bed
        File reference_fa
        File haps_vs_chm13_csv
        File? confident_bed
        Float? small_overlap
        Array[File]? evaluation_bed_small_overlap
        Float? large_overlap
        Array[File]? evaluation_bed_large_overlap
        Int min_sv_length
        String chromosome
        Boolean force_unique_reads
        String docker
        File? force_windows_bed
        Int n_cpu
        Int ram_size_gb
    }
    parameter_meta {
    }

    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"

    command <<<
        git --no-pager --git-dir ~{docker_dir}/sv_merge/.git log --decorate=short --pretty=oneline | head -n 1

        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        if ~{defined(confident_bed)}
        then
            for vcf in ~{sep=" " vcf_gz}; do
                # use bcftools to subset the vcf by the confident bed
                bcftools view -T ~{confident_bed} ${vcf} -Ov -o confident.vcf

                # convert to bgzipped vcf and overwrite the input VCF
                bcftools view -Oz -o ${vcf} confident.vcf

                # index the vcf
                bcftools index -t ${vcf}
            done
        fi

        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        # Downloading all the calls in $chromosome$.
        export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
        IDS_FILE=~{write_lines(vcf_id)}
        ADDRESSES_FILE=~{write_lines(vcf_gz)}
        paste -d , ${IDS_FILE} ${ADDRESSES_FILE} > list.txt
        cat list.txt
        VCFS=""; CLUSTER_BY=""; TOOLS=""; LINE="0";
        while read ROW; do
        TOOL_ID=${ROW%,*}
        TOOL_FILE=${ROW#*,}
        bcftools view --output-type z ${TOOL_FILE} ~{chromosome} > tmp.vcf.gz
        tabix -f tmp.vcf.gz
        bcftools norm --multiallelics - --output-type v tmp.vcf.gz > ${TOOL_ID}.vcf
        rm -f tmp.vcf.gz
        LINE=$(( ${LINE} + 1 ))
        if [[ ${LINE} -eq 1 ]]; then
            VCFS=${TOOL_ID}.vcf
            CLUSTER_BY=${TOOL_ID}.vcf
            TOOLS=${TOOL_ID}
        else
            VCFS="${VCFS},${TOOL_ID}.vcf"
            TOOLS="${TOOLS} ${TOOL_ID}"
        fi
        done < list.txt
        if ~{defined(force_windows_bed)} ; then
            awk -v chr="~{chromosome}" '{ if ( $1 == chr ) print $0 }' ~{force_windows_bed} > force.bed
            FORCE_WINDOWS_FLAG="--windows force.bed"
        else
            FORCE_WINDOWS_FLAG=""
        fi

        # Starting resource monitoring
        export MONITOR_MOUNT_POINT=~{work_dir}
        MONITOR_FILE="./~{chromosome}_monitor.log"
        bash ~{docker_dir}/vm_local_monitoring_script.sh &> ${MONITOR_FILE} &
        MONITOR_JOB=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        # Evaluating all VCFs
        EVALUATION_NAME="~{chromosome}_evaluation"
        rm -rf ./${EVALUATION_NAME}
        ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/evaluate \
        --n_threads ${N_THREADS} \
        --output_dir ~{work_dir}/${EVALUATION_NAME} \
        --bam_csv ~{haps_vs_chm13_csv} \
        --vcfs ${VCFS} \
        --cluster_by ${CLUSTER_BY} \
        --tandems ~{tandems_bed} \
        --ref ~{reference_fa} \
        --interval_max_length ~{interval_max_length} \
        --flank_length ~{flank_length} \
        --graphaligner_timeout ~{graphaligner_timeout} \
        ~{if force_unique_reads then "--force_unique_reads" else ""} \
        --min_sv_length ~{min_sv_length} \
        ${FORCE_WINDOWS_FLAG}

        # Ensure write buffers are flushed to disk
        sync

        # Stopping resource monitoring
        kill ${MONITOR_JOB}
        tail -n 100 ${MONITOR_FILE}
        df -h
        ls -laht

        # Analyzing the results of the evaluation
        ANALYSIS_NAME_SMALL="~{chromosome}_analysis_small"
        ANALYSIS_NAME_LARGE="~{chromosome}_analysis_large"
        if ~{defined(truth_vcf_id)} ; then
            TRUTH_ID_FLAG="--truth_id ~{truth_vcf_id}"
        else
            TRUTH_ID_FLAG=""
        fi
        if ~{defined(evaluation_bed_small_overlap)} ; then
            EVALUATION_BEDS_SMALL=~{sep=',' evaluation_bed_small_overlap}
            EVALUATION_BEDS_SMALL=$(echo ${EVALUATION_BEDS_SMALL} | tr ',' ' ')
            rm -rf ./${ANALYSIS_NAME_SMALL}
            ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/analyze_evaluation \
                --input_dir ~{work_dir}/${EVALUATION_NAME} \
                --output_dir ~{work_dir}/${ANALYSIS_NAME_SMALL} \
                --tools ${TOOLS} \
                --beds ${EVALUATION_BEDS_SMALL} \
                --min_bed_coverage ~{small_overlap} \
                ${TRUTH_ID_FLAG} &
        else
            # Make a placeholder so cromwell doesn't whine
            touch ~{work_dir}/${ANALYSIS_NAME_SMALL}.tar.gz
        fi
        if ~{defined(evaluation_bed_large_overlap)} ; then
            EVALUATION_BEDS_LARGE=~{sep=',' evaluation_bed_large_overlap}
            EVALUATION_BEDS_LARGE=$(echo ${EVALUATION_BEDS_LARGE} | tr ',' ' ')
            rm -rf ./${ANALYSIS_NAME_LARGE}
            ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/analyze_evaluation \
                --input_dir ~{work_dir}/${EVALUATION_NAME} \
                --output_dir ~{work_dir}/${ANALYSIS_NAME_LARGE} \
                --tools ${TOOLS} \
                --beds ${EVALUATION_BEDS_LARGE} \
                --min_bed_coverage ~{large_overlap} \
                ${TRUTH_ID_FLAG} &
        else
            # Make a placeholder so cromwell doesn't whine
            touch ~{work_dir}/${ANALYSIS_NAME_LARGE}.tar.gz
        fi
        if ! ~{defined(evaluation_bed_small_overlap)} && ! ~{defined(evaluation_bed_large_overlap)} ; then
            rm -rf ./${ANALYSIS_NAME_LARGE}
            ${TIME_COMMAND} ~{docker_dir}/sv_merge/build/analyze_evaluation \
                --input_dir ~{work_dir}/${EVALUATION_NAME} \
                --output_dir ~{work_dir}/${ANALYSIS_NAME_LARGE} \
                --tools ${TOOLS} \
                ${TRUTH_ID_FLAG} &
        fi

        wait

        # Outputting
        export GZIP=-1
        ${TIME_COMMAND} tar -czf ${EVALUATION_NAME}.tar.gz --exclude='*.fasta' --exclude='*.fa' ./${EVALUATION_NAME} &
        ${TIME_COMMAND} tar -czf ${ANALYSIS_NAME_SMALL}.tar.gz ./${ANALYSIS_NAME_SMALL} &
        ${TIME_COMMAND} tar -czf ${ANALYSIS_NAME_LARGE}.tar.gz ./${ANALYSIS_NAME_LARGE} &
        wait
    >>>

    output {
        File evaluation = work_dir + "/" + chromosome + "_evaluation.tar.gz"
        File analysis_small_overlap = work_dir + "/" + chromosome + "_analysis_small.tar.gz"
        File analysis_large_overlap = work_dir + "/" + chromosome + "_analysis_large.tar.gz"
        File monitor_log = work_dir + "/" + chromosome + "_monitor.log"
    }
    runtime {
        docker: docker
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk 128 HDD"
        preemptible: 0
    }
}