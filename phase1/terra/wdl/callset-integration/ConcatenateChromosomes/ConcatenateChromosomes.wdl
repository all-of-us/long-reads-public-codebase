########################################################################################
# ConcatenateChromosomes
# Phase 1 stage c (reference panel), step 6 of the callset-integration launch order.
#
# Concatenates the per-chromosome merged VCFs into one cohort VCF, with samples in the
# order given by samples_file.
#
# Runs once per cohort.
# Inputs:  source_dir (TruvariIntersample output), samples_file (one sample ID per line)
# Outputs: vcf_gz  the concatenated cohort SV VCF
# Previous: TruvariIntersample. Next: RemoveSamplesV11.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


# Simply concatenates per-chromosome intersample-merged VCFs into a single VCF.
#
workflow ConcatenateChromosomes {
    input {
        String source_dir
        File samples_file
    }
    parameter_meta {
        source_dir: "The output dir of workflow $TruvariIntersample$, containing sorted per-chromosome VCFs."
        samples_file: "One sample per line. This is the order of the samples in the output VCF."
    }

    call ConcatenateChromosomesImpl {
        input:
            source_dir = source_dir,
            samples_file = samples_file
    }
    
    output {
        File vcf_gz = ConcatenateChromosomesImpl.vcf_gz
        File vcf_gz_tbi = ConcatenateChromosomesImpl.vcf_gz_tbi
    }
}


task ConcatenateChromosomesImpl {
    input {
        String source_dir
        File samples_file
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    Int disk_size_gb = 256  # Arbitrary

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
        
        # Downloading
        while : ; do
            TEST=$(gsutil -m cp ~{source_dir}'/*.vcf.gz' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -maxdepth 1 -name '*.vcf.gz' | sort -V > list.txt
        
        # Ensuring that samples have the same order in all files
        rm -f list_filtered.txt
        while read FILE; do
            ID=$(basename ${FILE} .vcf.gz)
            bcftools view --samples-file ~{samples_file} --output-type z ${FILE} > ${ID}_filtered.vcf.gz
            tabix -f ${ID}_filtered.vcf.gz
            echo ${ID}_filtered.vcf.gz >> list_filtered.txt
            rm -f ${FILE}
        done < list.txt
        
        # Concatenating
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list_filtered.txt --output-type z > concat.vcf.gz
        tabix -f concat.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/concat.vcf.gz"
        File vcf_gz_tbi = work_dir + "/concat.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 16
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}