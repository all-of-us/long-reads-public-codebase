########################################################################################
# GetRegenotypedVcfKanpigV11
# Phase 1 stage b (per-sample filtering), step 2 of the callset-integration launch order.
#
# Re-genotypes one sample's collapsed SV VCF from its read alignments with Kanpig,
# using a sex-specific ploidy BED, and records which callers supported each SV
# (SUPP_PBSV, SUPP_SNIFFLES, SUPP_PAV).
#
# Runs once per sample.
# Inputs:  truvari_collapsed_vcf_gz (from TruvariIntrasampleV11), alignments_bam,
#          reference_fa, ploidy_bed_female / ploidy_bed_male, is_male / sex
# Outputs: regenotyped_kanpig  the re-genotyped VCF
# Previous: TruvariIntrasampleV11. Next: FilterIntrasampleDevV11.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


# Regeotypes a single-sample VCF.
#
workflow GetRegenotypedVcfKanpigV11 {
    input {
        Boolean is_male = false
        String sex = "F"
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    call GetRegenotypedVcfImpl {
        input:
            is_male = is_male,
            sex = sex,
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male
    }
    
    output {
        File regenotyped_kanpig = GetRegenotypedVcfImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = GetRegenotypedVcfImpl.regenotyped_kanpig_tbi
    }
}


task GetRegenotypedVcfImpl {
    input {
        Boolean is_male
        String sex
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "kanpig_regenotyped"
    String kanpig_params_singlesample = "--sizemin 20 --sizemax 10000 --chunksize 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 100*ceil(size(truvari_collapsed_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))
    
    Int n_cpu = 8
    Int mem_gb = 16  # 2*n_cpu suggested by Adam
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        df -h
        
        
        # Unpacks truvari's SUPP field into 3 INFO tags.
        #
        function transferSupp() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                if ($6=="0") printf("\t0\t0\t0");
                else if ($6=="1") printf("\t0\t0\t1");
                else if ($6=="2") printf("\t0\t1\t0");
                else if ($6=="3") printf("\t0\t1\t1");
                else if ($6=="4") printf("\t1\t0\t0");
                else if ($6=="5") printf("\t1\t0\t1");
                else if ($6=="6") printf("\t1\t1\t0");
                else if ($6=="7") printf("\t1\t1\t1");
                printf("\n"); \
            }' | bgzip -c > annotations.tsv.gz
            tabix -f -s1 -b2 -e2 annotations.tsv.gz
            rm -f header.txt
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' > header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> header.txt
            bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/SUPP_PAV,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV ${INPUT_VCF_GZ} --output-type z > ${OUTPUT_VCF_GZ}
            tabix -f ${OUTPUT_VCF_GZ}
        }


        # Main program
        touch ~{alignments_bai}
    
        # Formatting the merged VCF
        HAS_SUPP=$(bcftools view --header-only ~{truvari_collapsed_vcf_gz} | grep '##FORMAT=<ID=SUPP,' && echo 1 || echo 0)
        if [ ${HAS_SUPP} -eq 0 ]; then
            mv ~{truvari_collapsed_vcf_gz} merged.vcf.gz
            mv ~{truvari_collapsed_tbi} merged.vcf.gz.tbi
        else
            transferSupp ~{truvari_collapsed_vcf_gz} merged.vcf.gz
        fi

        # Making sure there is just one occurrence of '##fileformat=' in the
        # header (otherwise kanpig complains).
        rm -f cleaned.vcf*
        echo '##fileformat=VCFv4.2' > cleaned.vcf
        bcftools view --header-only merged.vcf.gz | grep -v '##fileformat=' >> cleaned.vcf
        bcftools view --no-header merged.vcf.gz >> cleaned.vcf
        bgzip cleaned.vcf
        tabix -f cleaned.vcf.gz
        rm -f merged.vcf.gz*

        # Kanpig single-sample
        if [ ~{is_male} == "true" -o ~{sex} == "M" ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        export RUST_BACKTRACE="full"
        ~{docker_dir}/kanpig --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} ~{kanpig_params_singlesample} --reference ~{reference_fa} --input cleaned.vcf.gz --bam ~{alignments_bam} --out tmp1.vcf.gz
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_kanpig = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: n_cpu
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
