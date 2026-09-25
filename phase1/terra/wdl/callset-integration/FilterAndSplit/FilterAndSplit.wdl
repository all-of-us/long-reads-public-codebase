########################################################################################
# FilterAndSplit
# Phase 1 stage b (per-sample filtering), step 4 of the callset-integration launch order.
#
# Copies each record's CALIBRATION_SENSITIVITY and caller-support annotations into
# FORMAT fields (so the inter-sample merge keeps them), optionally filters the sample's
# scored VCF with a bcftools expression, and writes one VCF per chromosome to
# destination_dir. The Phase 1 run left filter_string at "none"; the calibration-
# sensitivity filter was applied in TruvariIntersample instead.
#
# Runs once per sample.
# Inputs:  sample_vcf_gz (scored VCF from FilterIntrasampleDevV11), filter_string,
#          chromosomes, destination_dir
# Outputs: none declared; per-chromosome VCFs are written to destination_dir
# Previous: FilterIntrasampleDevV11. Next: TruvariIntersample.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


# Applies a filter to a VCF and sequentially splits it into a given set of 
# chromosomes.
#
workflow FilterAndSplit {
    input {
        String sample_id
        File sample_vcf_gz
        File sample_vcf_gz_tbi
        String filter_string = "none"
        Array[String] chromosomes
        String destination_dir
    }
    parameter_meta {
        sample_vcf_gz: "Assumed to be already sorted"
        filter_string: "String to be used in $bcftools filter --include$. 'none'=no filter."
        destination_dir: "The filtered and split files are stored in this remote directory."
    }

    call FilterAndSplitImpl {
        input:
            sample_id = sample_id,
            sample_vcf_gz = sample_vcf_gz,
            sample_vcf_gz_tbi = sample_vcf_gz_tbi,
            filter_string = filter_string,
            chromosomes = chromosomes,
            destination_dir = destination_dir
    }
    
    output {
    }
}


task FilterAndSplitImpl {
    input {
        String sample_id
        File sample_vcf_gz
        File sample_vcf_gz_tbi
        String filter_string
        Array[String] chromosomes
        String destination_dir
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    Int disk_size_gb = 10*(ceil(size(sample_vcf_gz,"GB")))  # Arbitrary

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Transfering annotations to the GT field, so that they are preserved by
        # the inter-sample merge later.
        echo '##FORMAT=<ID=CALIBRATION_SENSITIVITY,Number=1,Type=Float,Description="Calibration sensitivity according to the model applied by ScoreVariantAnnotations">' > header.txt
        echo '##FORMAT=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> header.txt
        echo '##FORMAT=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by sniffles">' >> header.txt
        echo '##FORMAT=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by pav">' >> header.txt
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CALIBRATION_SENSITIVITY\t%SUPP_PBSV\t%SUPP_SNIFFLES\t%SUPP_PAV\n' ~{sample_vcf_gz} | bgzip -c > annotations.tsv.gz
        tabix -s1 -b2 -e2 annotations.tsv.gz
        bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,FORMAT/CALIBRATION_SENSITIVITY,FORMAT/SUPP_PBSV,FORMAT/SUPP_SNIFFLES,FORMAT/SUPP_PAV ~{sample_vcf_gz} --output-type z > formatted.vcf.gz
        tabix -f formatted.vcf.gz
        
        # Splitting
        CHROMOSOMES=~{sep='-' chromosomes}
        CHROMOSOMES=$(echo ${CHROMOSOMES} | tr '-' ' ')
        FILTER_STRING="~{filter_string}"
        if [ ${FILTER_STRING} = none ]; then
            INCLUDE_STR=""
        else
            INCLUDE_STR="--include ${FILTER_STRING}"
        fi
        for CHROMOSOME in ${CHROMOSOMES}; do
            ${TIME_COMMAND} bcftools filter ${INCLUDE_STR} --regions ${CHROMOSOME} --output-type z formatted.vcf.gz > ~{sample_id}_${CHROMOSOME}_split.vcf.gz
            tabix -f ~{sample_id}_${CHROMOSOME}_split.vcf.gz
        done
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*_split.vcf.gz*' ~{destination_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>

    output {
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}