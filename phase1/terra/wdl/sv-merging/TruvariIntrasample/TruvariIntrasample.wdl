########################################################################################
# TruvariIntrasample
# Helper from the SV-merging comparisons, not on the Phase 1 launch path. The Phase 1
# intra-sample merge is callset-integration/TruvariIntrasampleV11.
#
# Merges one sample's PBSV, Sniffles2, and PAV SV callsets and collapses them with truvari.
#
# Inputs:  pbsv_vcf_gz, sniffles_vcf_gz, pav_vcf_gz, reference_fa, sample_id
# Outputs: truvari_collapsed, bcftools_merged
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

# Truvari intra-merge for AoU SV
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
        
        # Removing multiallelic records from the input
        bcftools norm --multiallelics - --output-type z ~{pbsv_vcf_gz} > pbsv_new.vcf.gz
        tabix pbsv_new.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{sniffles_vcf_gz} > sniffles_new.vcf.gz
        tabix sniffles_new.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{pav_vcf_gz} > pav_new.vcf.gz
        tabix pav_new.vcf.gz
        rm -f tmp.vcf.gz*
        
        # Step 1 - clean up the VCFs
        # - Assigns quality scores to each SV caller's result
        #  - pav 4
        #  - pbsv 3
        #  - sniffles 2
        # - Resolves any symbolic alts (e.g. `<DEL>` with the sequence from the
        #   reference)
        #   - symbolic variants are given quality score of 1
        # - Filters out `<CNV>` from pbsv and `<INS>` from sniffles
        # - Filters out variants greater than 100Mbp
        # - Fills in blank genotypes with `0`
        # - Filters out BND variants
        # - Turn some deletion/insertion pairs into inversions
        # The quality scores are set based on which variant representations we
        # believe are generally more accurate with higher being better.
        mkdir -p preprocessed
        for in_vcf in pav_new.vcf.gz pbsv_new.vcf.gz sniffles_new.vcf.gz
        do
            prename=preprocessed/pre_inv_$(basename $in_vcf)
            python ~{docker_dir}/resolve.py ${in_vcf} ~{reference_fa} \
                | bcftools norm --check-ref s --fasta-ref ~{reference_fa} -N -m-any \
                | bcftools view -i "SVTYPE != 'BND'" -O z -o ${prename}
            tabix $prename
            outname=preprocessed/$(basename $in_vcf)
            ${TIME_COMMAND} python ~{docker_dir}/inversion_guesser.py -i $prename -o $outname
        done

        # Step 2 - merge
        # Pastes the samples together in the order of the preferred genotypes.
        # That is to say, this creates a three sample VCF with sample columns
        # from pbsv, sniffles, pav_sv
        bcftools merge --threads ${N_THREADS} --merge none --force-samples -O z \
            -o tmp.vcf.gz \
            preprocessed/pbsv_new.vcf.gz \
            preprocessed/sniffles_new.vcf.gz \
            preprocessed/pav_new.vcf.gz 
        tabix tmp.vcf.gz
        
        # Removing multiallelic records (we observed that they are created in
        # Step 2 sometimes).
        bcftools norm --multiallelics - --output-type z tmp.vcf.gz > ~{sample_id}.bcftools_merged.vcf.gz
        tabix ~{sample_id}.bcftools_merged.vcf.gz
        rm -f tmp.vcf.gz*

        # Step 3 - collapse
        truvari collapse -i ~{sample_id}.bcftools_merged.vcf.gz -c removed.vcf.gz \
            --sizemin 0 --sizemax 1000000 -k maxqual --gt het --intra --pctseq 0.90 --pctsize 0.90 --refdist 500 \
            | bcftools sort --max-mem 8G -O z -o ~{sample_id}.truvari_collapsed.vcf.gz
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
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
