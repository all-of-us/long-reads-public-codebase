version 1.0


# Workflow for intra-sample merging for AoU.
#
workflow TruvariIntrasample {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call TruvariIntrasampleImpl {
        input:
            sample_id = sample_id,
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File truvari_collapsed = TruvariIntrasampleImpl.truvari_collapsed
    	File truvari_collapsed_idx = TruvariIntrasampleImpl.truvari_collapsed_idx
    	File bcftools_merged = TruvariIntrasampleImpl.bcftools_merged
    	File bcftools_merged_idx = TruvariIntrasampleImpl.bcftools_merged_idx
    }
}


# Other intermediate files created, but likely aren't useful for production are:
# - preprocessed/ directory for each caller's cleaned result
# - ~{sample_id}.removed.vcf.gz variant representations removed during
# collapsing
#
task TruvariIntrasampleImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(pbsv_vcf_gz,"GB")) + ceil(size(sniffles_vcf_gz,"GB")) + ceil(size(pav_vcf_gz,"GB")) + ceil(size(reference_fa,"GB")) ) + 50
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Step 1 - clean up the VCFs
        # - Assigns quality scores to each SV caller's result
        #  - pav 4
        #  - pbsv 3
        #  - sniffles 2
        # - Resolves any symbolic alts (e.g. `<DEL>` with the sequence from the
        #   reference)
        #   - symbolic variants are given quality score of 1
        # - Filters out variants greater than 100kbp (I might want to remove
        #   this)
        # - Fills in blank genotypes with `0`
        # - Filters out BND variants
        # The quality scores are set based on which variant representations we
        # believe are generally more accurate with higher being better.
        mkdir -p preprocessed
        for in_vcf in ~{pav_vcf_gz} ~{pbsv_vcf_gz} ~{sniffles_vcf_gz}
        do
            outname=preprocessed/$(basename $in_vcf)
            python ~{docker_dir}/resolve.py ${in_vcf} ~{reference_fa} \
                | bcftools norm --check-ref s --fasta-ref ~{reference_fa} -N -m-any \
                | bcftools view -i "SVTYPE != 'BND'" -O z -o ${outname}
            tabix $outname
        done

        # Step 2 - merge
        # Pastes the samples together in the order of the preferred genotypes.
        # That is to say, this creates a three sample VCF with sample columns
        # from  pav_sv, pbsv, and sniffles.
        bcftools merge --threads ${N_THREADS} --merge none --force-samples -O z \
            -o ~{sample_id}.bcftools_merged.vcf.gz \
            preprocessed/$(basename ~{pav_vcf_gz}) \
            preprocessed/$(basename ~{pbsv_vcf_gz}) \
            preprocessed/$(basename ~{sniffles_vcf_gz}) 
        tabix ~{sample_id}.bcftools_merged.vcf.gz

        # Step 3 - collapse
        truvari collapse -i ~{sample_id}.bcftools_merged.vcf.gz -c removed.vcf.gz -k maxqual --gt --intra \
            --pctseq 0.90 --pctsize 0.90 --refdist 500 \
            | bcftools sort -O z -o ~{sample_id}.truvari_collapsed.vcf.gz
        tabix ~{sample_id}.truvari_collapsed.vcf.gz
    >>>
    
    output {
    	File truvari_collapsed = "~{work_dir}/~{sample_id}.truvari_collapsed.vcf.gz"
    	File truvari_collapsed_idx = "~{work_dir}/~{sample_id}.truvari_collapsed.vcf.gz.tbi"
    	File bcftools_merged = "~{work_dir}/~{sample_id}.bcftools_merged.vcf.gz"
    	File bcftools_merged_idx = "~{work_dir}/~{sample_id}.bcftools_merged.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
