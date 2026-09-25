########################################################################################
# MergeRegenotypedIntersampleVcf
# Phase 1 stage c (reference panel), step 9 of the callset-integration launch order.
#
# Pastes the per-sample re-genotyped VCFs (which all share the same records) side by
# side into one multi-sample cohort VCF. This is the joint SV callset
# (concat_regenotyped.vcf.gz) used to build the reference panel.
#
# Runs once per cohort.
# Inputs:  regenotyped_vcfs_list  a file listing the gs:// paths of the per-sample VCFs
# Outputs: output_vcf_gz  the cohort SV VCF
# Previous: GetRegenotypedVcfKanpigMergedV11 (all samples). Next: PhysicalPhasing.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every original
# sample. The program combines all such re-genotyped single-sample VCFs into
# an inter-sample-merged VCF with updated genotypes.
#
# Remark: every file in $regenotyped_vcfs_list$ is assumed to have exactly the
# same set of calls in the same order. This is usually the case if the files
# are the result of re-genotyping the same inter-sample VCF with different BAMs.
# If this is not the case the program gives a wrong output, and it should be
# replaced with $bcftools merge --merge none$.
#
workflow MergeRegenotypedIntersampleVcf {
    input {
        File regenotyped_vcfs_list
    }
    parameter_meta {
    }
    
    call MergeImpl {
        input:
            regenotyped_vcfs_list = regenotyped_vcfs_list
    }
    
    output {
        File output_vcf_gz = MergeImpl.output_vcf_gz
        File output_tbi = MergeImpl.output_tbi
    }
}


task MergeImpl {
    input {
        File regenotyped_vcfs_list
        Int n_cpus
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        function pasteThread() {
            local THREAD_ID=$1
            
            OUTPUT_FILE="columns_${THREAD_ID}.txt"; touch ${OUTPUT_FILE}
            FIELDS_FILE="fields_${THREAD_ID}.txt"; touch ${FIELDS_FILE}
            TMP_PREFIX="tmp_${THREAD_ID}"
            i="0";
            while read ADDRESS; do
                i=$(( $i + 1 ))
                # Adding the new sample to the set of columns
                ${TIME_COMMAND} gsutil -m cp ${ADDRESS} ${TMP_PREFIX}.vcf.gz
                bcftools view --header-only ${TMP_PREFIX}.vcf.gz > ${TMP_PREFIX}.txt
                N_ROWS=$(wc -l < ${TMP_PREFIX}.txt)
                tail -n 1 ${TMP_PREFIX}.txt | cut -f 10 > sample_${THREAD_ID}.txt
                if [ $i = "1" ]; then
                    mv sample_${THREAD_ID}.txt ${FIELDS_FILE}
                else
                    paste ${FIELDS_FILE} sample_${THREAD_ID}.txt > ${FIELDS_FILE}.prime
                    mv ${FIELDS_FILE}.prime ${FIELDS_FILE}
                fi
                echo "Current fields of thread ${THREAD_ID}:"; cat ${FIELDS_FILE}
                # Adding the new column to the body
                bcftools view --no-header ${TMP_PREFIX}.vcf.gz | cut -f 10 > ${TMP_PREFIX}.txt
                if [ $i = "1" ]; then
                    mv ${TMP_PREFIX}.txt ${OUTPUT_FILE}
                else
                    ${TIME_COMMAND} paste ${OUTPUT_FILE} ${TMP_PREFIX}.txt > ${OUTPUT_FILE}.prime
                    mv ${OUTPUT_FILE}.prime ${OUTPUT_FILE}
                fi
            done < list_${THREAD_ID}
        }
        
        
        # Main program
        
        # Initializing the inter-sample VCF with the first file
        ADDRESS=$(head -n 1 ~{regenotyped_vcfs_list})
        gsutil -m cp ${ADDRESS} first.vcf.gz
        bcftools view --header-only first.vcf.gz > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        tail -n 1 tmp.txt | cut -f 1,2,3,4,5,6,7,8,9 > fields.txt
        bcftools view --no-header first.vcf.gz | cut -f 1,2,3,4,5,6,7,8,9 > calls.txt
        rm -f first.vcf.gz
        
        # Appending all the remaining files
        N_ROWS=$(wc -l < ~{regenotyped_vcfs_list})
        N_ROWS=$(( ${N_ROWS} / ${N_THREADS} ))
        split -d -l ${N_ROWS} ~{regenotyped_vcfs_list} list_
        COLUMNS_FILES=""; FIELDS_FILES=""
        for LIST_FILE in $(find . -maxdepth 1 -name 'list_*' | sort); do
            ID=${LIST_FILE#./list_}
            pasteThread ${ID} &
            COLUMNS_FILES="${COLUMNS_FILES} columns_${ID}.txt"
            FIELDS_FILES="${FIELDS_FILES} fields_${ID}.txt"
        done
        wait
        paste fields.txt ${FIELDS_FILES} > fields_all.txt
        paste calls.txt ${COLUMNS_FILES} > body.txt
        cat header.txt fields_all.txt body.txt > merged.vcf
        rm -f *.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: n_cpus
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}