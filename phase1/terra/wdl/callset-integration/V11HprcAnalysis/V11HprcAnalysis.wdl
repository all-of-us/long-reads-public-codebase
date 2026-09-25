########################################################################################
# V11HprcAnalysis
# Evaluation, not on the Phase 1 launch path (HPRC benchmarking for the manuscript).
#
# Merges the HPRC samples' assembly-based (dipcall) SV calls into a truth VCF, keeping
# events of at least dipcall_min_sv_length bp inside each sample's confident regions,
# and compares the unfiltered, calibration-sensitivity 0.7, and 0.9 callsets against it
# with truvari bench, genome-wide and outside tandem repeats (tr_bed). Results are
# written as tar.gz archives of the truvari output directories.
#
# Runs once. There is no example_inputs.json because no Phase 1 run was found.
# Inputs:  dipcall_vcf_gz / dipcall_confident_bed / dipcall_sample_id (HPRC truth),
#          unfiltered_vcf_gz, cal_sens_07_vcf_gz, cal_sens_09_vcf_gz, tr_bed, reference_fa
# Outputs: none declared; comparison results are written by the Compare task
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
########################################################################################

version 1.0


#
workflow V11HprcAnalysis {
    input {
        File reference_fa
        File reference_fai
        File tr_bed
        
        Array[String] dipcall_sample_id
        Array[File] dipcall_vcf_gz
        Array[File] dipcall_tbi
        Array[File] dipcall_confident_bed
        Int dipcall_min_sv_length = 50
        
        File unfiltered_vcf_gz
        File unfiltered_tbi
        File cal_sens_07_vcf_gz
        File cal_sens_07_tbi
        File cal_sens_09_vcf_gz
        File cal_sens_09_tbi
        
        String truvari_collapse_params = "--sizemin 50 --sizemax 1000000"
        String nontr_fraction = "0.9"
        String truvari_bench_params_global = "--pick multi --sizefilt 50 --sizemax 1000000 --sizemin 50 --pctsize 0.9 --pctseq 0"
        String truvari_bench_params_nontr = "--pick multi --sizefilt 50 --sizemax 1000000 --sizemin 50 --pctsize 0.9 --pctseq 0.9"
    }
    parameter_meta {
    }

    call BcftoolsMergeDipcall {
        input:
            sample_id = dipcall_sample_id,
            input_vcf_gz = dipcall_vcf_gz,
            input_tbi = dipcall_tbi,
            confident_bed = dipcall_confident_bed,
            min_sv_length = dipcall_min_sv_length
    }
    call Compare as CompareUnfiltered {
        input:
            query_vcf_gz = unfiltered_vcf_gz,
            query_tbi = unfiltered_tbi,
            truth_vcf_gz = BcftoolsMergeDipcall.output_vcf_gz,
            truth_tbi = BcftoolsMergeDipcall.output_tbi,
            tr_bed = tr_bed,
            reference_fai = reference_fai,
            truvari_bench_params_global = truvari_bench_params_global,
            truvari_bench_params_nontr = truvari_bench_params_nontr,
            nontr_fraction = nontr_fraction
    }
    call Compare as Compare07 {
        input:
            query_vcf_gz = cal_sens_07_vcf_gz,
            query_tbi = cal_sens_07_tbi,
            truth_vcf_gz = BcftoolsMergeDipcall.output_vcf_gz,
            truth_tbi = BcftoolsMergeDipcall.output_tbi,
            tr_bed = tr_bed,
            reference_fai = reference_fai,
            truvari_bench_params_global = truvari_bench_params_global,
            truvari_bench_params_nontr = truvari_bench_params_nontr,
            nontr_fraction = nontr_fraction
    }
    call Compare as Compare09 {
        input:
            query_vcf_gz = cal_sens_09_vcf_gz,
            query_tbi = cal_sens_09_tbi,
            truth_vcf_gz = BcftoolsMergeDipcall.output_vcf_gz,
            truth_tbi = BcftoolsMergeDipcall.output_tbi,
            tr_bed = tr_bed,
            reference_fai = reference_fai,
            truvari_bench_params_global = truvari_bench_params_global,
            truvari_bench_params_nontr = truvari_bench_params_nontr,
            nontr_fraction = nontr_fraction
    }

    
    output {
    }
}


task BcftoolsMergeDipcall {
    input {
        Array[String] sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        Array[File] confident_bed
        Int min_sv_length
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int n_files = length(input_vcf_gz)
    
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
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        SAMPLE_IDS=~{sep=',' sample_id}
        CONFIDENT_BEDS=~{sep=',' confident_bed}
        rm -f list.txt
        for i in $(seq 1 ~{n_files}); do
            INPUT_FILE=$(echo ${INPUT_FILES} | cut -d , -f ${i})
            SAMPLE_ID=$(echo ${SAMPLE_IDS} | cut -d , -f ${i})
            CONFIDENT_BED=$(echo ${CONFIDENT_BEDS} | cut -d , -f ${i})
            # Enforcing the right sample name
            echo ${SAMPLE_ID} > samples.txt
            bcftools reheader --samples samples.txt ${INPUT_FILE} > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
            # Restricting to the confident BED
            bcftools view --regions-file ${CONFIDENT_BED} --regions-overlap pos --output-type z tmp1.vcf.gz > tmp2.vcf.gz
            tabix -f tmp2.vcf.gz
            rm -f tmp1.vcf.gz*
            # Removing multiallelic records
            bcftools norm --multiallelics - --output-type v tmp2.vcf.gz > tmp3.vcf
            rm -f tmp2.vcf.gz*
            # Keeping only long events
            java -cp ~{docker_dir} Dipcall2VCF tmp3.vcf ~{min_sv_length} ${SAMPLE_ID}.fixed.vcf
            rm -f tmp3.vcf
            bgzip ${SAMPLE_ID}.fixed.vcf
            tabix -f ${SAMPLE_ID}.fixed.vcf.gz
            echo ${SAMPLE_ID}.fixed.vcf.gz >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        # Removing multiallelic records
        bcftools norm --multiallelics - --output-type v tmp1.vcf.gz > tmp2.vcf
        rm -f tmp1.vcf.gz*
        # Keeping only long events
        java -cp ~{docker_dir} Dipcall2VCF tmp2.vcf ~{min_sv_length} merged.vcf
        rm -f tmp2.vcf
        bgzip merged.vcf
        tabix -f merged.vcf.gz
        ls -laht; tree
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}


task Compare {
    input {
        File query_vcf_gz
        File query_tbi
        File truth_vcf_gz
        File truth_tbi
        File tr_bed
        File reference_fai
        String truvari_bench_params_global
        String truvari_bench_params_nontr
        String nontr_fraction
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/aou"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        function format() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_FILE=$2
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --output-type z ${INPUT_VCF_GZ} > tmp.vcf.gz
            tabix tmp.vcf.gz
            bcftools view --header-only tmp.vcf.gz > header.txt
            N_ROWS=$(wc -l < header.txt)
            head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
            echo -e  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> out.vcf
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> out.vcf
            bcftools view --no-header tmp.vcf.gz | awk '{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT\t0/1\n",$1,$2,$3,$4,$5,$6,$7,$8); }' >> out.vcf
            bgzip -c out.vcf > ${OUTPUT_FILE}
            tabix ${OUTPUT_FILE}
            rm -f tmp.vcf* header.txt out.vcf*
        }
        
        
        format ~{query_vcf_gz} query.vcf.gz
        format ~{truth_vcf_gz} truth.vcf.gz
        
        # All calls
        ${TIME_COMMAND} truvari bench ~{truvari_bench_params_global} -b truth.vcf.gz -c query.vcf.gz -o ./truvari_dir_global
        tar -czf output_global.tar.gz ./truvari_dir_global/
        
        # Outside TRs
        bedtools sort -faidx ~{reference_fai} -i ~{tr_bed} > sorted.bed
        bedtools complement -i sorted.bed -g ~{reference_fai} > complement.bed
        bcftools view --header-only query.vcf.gz > query_nontr.vcf
        bedtools intersect -u -f ~{nontr_fraction} -a query.vcf.gz -b complement.bed >> query_nontr.vcf
        bgzip query_nontr.vcf
        tabix -f query_nontr.vcf.gz
        bcftools view --header-only truth.vcf.gz > truth_nontr.vcf
        bedtools intersect -u -f ~{nontr_fraction} -a truth.vcf.gz -b complement.bed >> truth_nontr.vcf
        bgzip truth_nontr.vcf
        tabix -f truth_nontr.vcf.gz
        ${TIME_COMMAND} truvari bench ~{truvari_bench_params_nontr} -b truth_nontr.vcf.gz -c query_nontr.vcf.gz -o ./truvari_dir_nontr
        tar -czf output_nontr.tar.gz ./truvari_dir_nontr/
    >>>
    
    output {
        File truvari_output_global = work_dir + "/output_global.tar.gz"
        File truvari_output_nontr = work_dir + "/output_nontr.tar.gz"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: "32GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}


task TruvariCollapseDipcall {
    input {
        File input_vcf_gz
        File input_tbi
        String truvari_params
        File reference_fa
        File reference_fai
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} truvari collapse ~{truvari_params} --reference ~{reference_fa} --input ~{input_vcf_gz} > truvari_collapsed.vcf | bcftools sort --max-mem $(( ~{ram_size_gb} - 2 ))G --output-type z > truvari.vcf.gz
        tabix -f truvari.vcf.gz
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/truvari.vcf.gz"
        File output_tbi = work_dir + "/truvari.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 2
        memory: ram_size_gb + "GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
