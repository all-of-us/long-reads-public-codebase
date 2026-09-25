########################################################################################
# GLIMPSEBatchedCasePerChromosomeSingleBatch (the WDL names its workflow GLIMPSEBatchedCasePerChromosome)
# Phase 1 stage d (genotyping and imputation), launch step 10.
#
# For one batch of short-read samples: merges their KAGE genotypes, then phases and
# imputes them against the panel with GLIMPSE (chunk, phase, ligate, sample), per
# chromosome, and resolves colliding alleles in the result.
#
# Runs once per batch of samples.
# Inputs:  chromosome_kage_vcf_gzs_txts (from KAGECasePerChromosomeFlexscattered),
#          sample_names, chromosomes, genetic_maps, panel_split_vcf_gz,
#          fix_variant_collisions_java and its settings, kage_docker
# Outputs: kage_vcf_gz                  the batch's KAGE genotypes
#          glimpse_vcf_gz               GLIMPSE phased and imputed genotypes
#          phased_collisionless_vcf_gz  the same with collisions resolved
# Previous: KAGECasePerChromosomeFlexscattered. Next: HierarchicallyMergeVcfs.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/methods/phasing/HierarchicallyMergeVcfs.wdl" as HierarchicallyMergeVcfs

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

workflow GLIMPSEBatchedCasePerChromosome {
    input {
        Array[File] chromosome_kage_vcf_gzs_txts         # one txt list of per-chromosome outputs per sample
        Array[File] chromosome_kage_vcf_gz_tbis_txts     # one txt list of per-chromosome outputs per sample
        Array[String] sample_names

        # per chromosome
        Array[String]+ chromosomes
        Array[File]+ genetic_maps
        Array[File] panel_split_vcf_gz # for GLIMPSE
        Array[File] panel_split_vcf_gz_tbi

        String extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"
        String extra_phase_args = ""
        String output_prefix

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        String kage_docker
        File? monitoring_script

        RuntimeAttributes merge_runtime_attributes = {"use_ssd": true}
        RuntimeAttributes concat_runtime_attributes = {"use_ssd": true}
        RuntimeAttributes glimpse_phase_runtime_attributes = {}
        RuntimeAttributes glimpse_sample_runtime_attributes = {}
        Map[String, Int]? chromosome_to_glimpse_command_mem_gb      # for running per-chromosome by choosing large chunk size; this will override glimpse_phase_runtime_attributes
        Int? glimpse_phase_preemptible
    }

    call CreateTables {
        input:
            vcfs_txts = chromosome_kage_vcf_gzs_txts,
            tbis_txts = chromosome_kage_vcf_gz_tbis_txts,
            prefix = output_prefix,
            docker = kage_docker
    }

	Int b = 0		# keep filename conventions the same
    scatter (j in range(length(chromosomes))) {
        call GLIMPSEChunk as ChromosomeGLIMPSEChunk {
            input:
                vcf = panel_split_vcf_gz[j],
                tbi = panel_split_vcf_gz_tbi[j],
                region = chromosomes[j],
                prefix = "chunk." + chromosomes[j],
                extra_chunk_args = extra_chunk_args
        }
        
        String chromosome = chromosomes[j]
        Array[String] chromosome_kage_vcf_gzs = transpose(read_tsv(CreateTables.chromosome_kage_vcf_gzs_tsv))[j]
        Array[String] chromosome_kage_vcf_gz_tbis = transpose(read_tsv(CreateTables.chromosome_kage_vcf_gz_tbis_tsv))[j]

		# switched to non-hierarchical merge
        call HierarchicallyMergeVcfs.Ivcfmerge as BatchChromosomeKAGEMergeAcrossSamples {
            input:
                vcf_gzs = chromosome_kage_vcf_gzs,
                vcf_gz_tbis = chromosome_kage_vcf_gz_tbis,
                sample_names = sample_names,
                output_prefix = output_prefix + ".batch-" + b + "." + chromosome + ".kage",
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        Array[String] input_regions = read_lines(ChromosomeGLIMPSEChunk.input_regions)
        Array[String] output_regions = read_lines(ChromosomeGLIMPSEChunk.output_regions)

        if (defined(chromosome_to_glimpse_command_mem_gb)) {
            Int command_mem_gb = select_first([chromosome_to_glimpse_command_mem_gb])[chromosome]
        }
        scatter (k in range(length(input_regions))) {
            call GLIMPSEPhase as BatchChunkedGLIMPSEPhase {
                input:
                    kage_vcf_gz = BatchChromosomeKAGEMergeAcrossSamples.merged_vcf_gz,
                    kage_vcf_gz_tbi = BatchChromosomeKAGEMergeAcrossSamples.merged_vcf_gz_tbi,
                    panel_split_vcf_gz = panel_split_vcf_gz[j],
                    panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi[j],
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    genetic_map = genetic_maps[j],
                    output_prefix = output_prefix + ".batch-" + b + "." + chromosome + ".shard-" + k + ".phased",
                    extra_phase_args = extra_phase_args,
                    docker = kage_docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = glimpse_phase_runtime_attributes,
                    command_mem_gb = command_mem_gb,
                    preemptible = glimpse_phase_preemptible
            }
        }

        call GLIMPSELigate as BatchChromosomeGLIMPSELigate {
            input:
                vcfs = BatchChunkedGLIMPSEPhase.phased_vcf_gz,
                vcf_idxs = BatchChunkedGLIMPSEPhase.phased_vcf_gz_tbi,
                prefix = output_prefix + ".batch-" + b + "." + chromosome + ".ligated",
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        call GLIMPSESample as BatchChromosomeGLIMPSESample {
            input:
                ligated_vcf_gz = BatchChromosomeGLIMPSELigate.ligated_vcf_gz,
                ligated_vcf_gz_tbi = BatchChromosomeGLIMPSELigate.ligated_vcf_gz_tbi,
                output_prefix = output_prefix + ".batch-" + b + "." + chromosome + ".sampled",
                docker = kage_docker,
                monitoring_script = monitoring_script,
                runtime_attributes = glimpse_sample_runtime_attributes
        }

        call FixVariantCollisions as BatchChromosomeGLIMPSEFixVariantCollisions { input:
            vcf_gz = BatchChromosomeGLIMPSESample.glimpse_vcf_gz,
            vcf_gz_tbi = BatchChromosomeGLIMPSESample.glimpse_vcf_gz_tbi,
            fix_variant_collisions_java = fix_variant_collisions_java,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + ".batch-" + b + "." + chromosome + ".kage.glimpse.collisionless",
            monitoring_script = monitoring_script
        }
    }

    # concat across chromosomes
    call ConcatVcfs as KAGEConcatVcfs {
        input:
            vcf_gzs = BatchChromosomeKAGEMergeAcrossSamples.merged_vcf_gz,
            vcf_gz_tbis = BatchChromosomeKAGEMergeAcrossSamples.merged_vcf_gz_tbi,
            output_prefix = output_prefix + ".kage",
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    call ConcatVcfs as GLIMPSEUnphasedConcatVcfs {
        input:
            vcf_gzs = BatchChromosomeGLIMPSELigate.ligated_vcf_gz,
            vcf_gz_tbis = BatchChromosomeGLIMPSELigate.ligated_vcf_gz_tbi,
            output_prefix = output_prefix + ".kage.glimpse.unphased",
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    call ConcatVcfs as GLIMPSEConcatVcfs {
        input:
            vcf_gzs = BatchChromosomeGLIMPSESample.glimpse_vcf_gz,
            vcf_gz_tbis = BatchChromosomeGLIMPSESample.glimpse_vcf_gz_tbi,
            output_prefix = output_prefix + ".kage.glimpse",
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    call ConcatVcfs as PhasedCollisionlessConcatVcfs {
        input:
            vcf_gzs = BatchChromosomeGLIMPSEFixVariantCollisions.collisionless_vcf_gz,
            vcf_gz_tbis = BatchChromosomeGLIMPSEFixVariantCollisions.collisionless_vcf_gz_tbi,
            output_prefix = output_prefix + ".kage.glimpse.collisionless",
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    output {
        File kage_vcf_gz = KAGEConcatVcfs.vcf_gz
        File kage_vcf_gz_tbi = KAGEConcatVcfs.vcf_gz_tbi
        File glimpse_unphased_vcf_gz = GLIMPSEUnphasedConcatVcfs.vcf_gz
        File glimpse_unphased_vcf_gz_tbi = GLIMPSEUnphasedConcatVcfs.vcf_gz_tbi
        File glimpse_vcf_gz = GLIMPSEConcatVcfs.vcf_gz
        File glimpse_vcf_gz_tbi = GLIMPSEConcatVcfs.vcf_gz_tbi
        File phased_collisionless_vcf_gz = PhasedCollisionlessConcatVcfs.vcf_gz
        File phased_collisionless_vcf_gz_tbi = PhasedCollisionlessConcatVcfs.vcf_gz_tbi
    }
}

# input: arrays of per-sample list files (each list file is one per-chromosome path per line)
# output: tsvs of sample x chromosome paths (one row per sample, one column per chromosome)
task CreateTables {
    input {
        Array[File] vcfs_txts
        Array[File] tbis_txts
        String prefix

        String docker

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        touch ~{prefix}.chromosome_kage_vcf_gzs.tsv
        for txt in ~{sep=' ' vcfs_txts}; do
            paste -s $txt >> ~{prefix}.chromosome_kage_vcf_gzs.tsv
        done

        touch ~{prefix}.chromosome_kage_vcf_gzs.tsv
        for txt in ~{sep=' ' tbis_txts}; do
            paste -s $txt >> ~{prefix}.chromosome_kage_vcf_gz_tbis.tsv
        done
    }

    output {
        File chromosome_kage_vcf_gzs_tsv = "~{prefix}.chromosome_kage_vcf_gzs.tsv"
        File chromosome_kage_vcf_gz_tbis_tsv = "~{prefix}.chromosome_kage_vcf_gz_tbis.tsv"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 3]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task CreateBatches {
    input {
        Array[Array[String]] sample_by_chromosome_vcf_gzs
        Array[Array[String]] sample_by_chromosome_vcf_gz_tbis
        Array[String] sample_names
        Int batch_size

        String docker
        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        cat ~{write_tsv(sample_by_chromosome_vcf_gzs)} | split -l ~{batch_size} - chromosome_vcf_gz_batch_
        cat ~{write_tsv(sample_by_chromosome_vcf_gz_tbis)} | split -l ~{batch_size} - chromosome_vcf_gz_tbi_batch_
        cat ~{write_lines(sample_names)} | split -l ~{batch_size} - sample_name_batch_
    }

    output {
        Array[File] chromosome_vcf_gz_batch_files = glob("chromosome_vcf_gz_batch_*")
        Array[File] chromosome_vcf_gz_tbi_batch_files = glob("chromosome_vcf_gz_tbi_batch_*")
        Array[File] sample_name_batch_files = glob("sample_name_batch_*")
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 3]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 10]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task ConcatVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

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

        mkdir inputs
        mv ~{sep=' ' vcf_gzs} inputs
        mv ~{sep=' ' vcf_gz_tbis} inputs

        # TODO FIX LEXICOGRAPHICAL BUG!
        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz) --naive -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
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

task GLIMPSEChunk {
    input {
        File vcf
        File tbi
        String region
        String prefix
        String? extra_chunk_args

        RuntimeAttributes runtime_attributes = {}
    }

    Int disk_size_gb = 2*ceil(size([vcf, tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > input-regions.txt
        cut -f 4 chunks.txt > output-regions.txt
    >>>

    output {
        File input_regions = "input-regions.txt"
        File output_regions = "output-regions.txt"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
        cpu: select_first([runtime_attributes.cpu, 4])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 10])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task GLIMPSEPhase {
    input {
        File kage_vcf_gz
        File kage_vcf_gz_tbi
        File panel_split_vcf_gz       # for GLIMPSE
        File panel_split_vcf_gz_tbi
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String? extra_phase_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
        Int? command_mem_gb = 7
        Int? preemptible = 2
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view --no-version -r ~{input_region},~{output_region} ~{kage_vcf_gz} | \
            sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | \
            bcftools view -Ob -o ~{output_prefix}.kage.nonan.bcf
        bcftools index ~{output_prefix}.kage.nonan.bcf

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        ./GLIMPSE_phase_static \
            -I ~{output_prefix}.kage.nonan.bcf \
            -R ~{panel_split_vcf_gz} \
            --input-region ~{input_region} \
            --output-region ~{output_region} \
            --map ~{genetic_map} \
            --input-GL \
            --thread $(nproc) \
            ~{extra_phase_args} \
            --output ~{output_prefix}.kage.glimpse.raw.vcf.gz

        # take KAGE VCF header and add GLIMPSE INFO and FORMAT lines (GLIMPSE header only contains a single chromosome and breaks bcftools concat --naive)
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^##' > kage.header.txt
        bcftools view --no-version -h ~{output_prefix}.kage.glimpse.raw.vcf.gz | grep -E '^##INFO|^##FORMAT|^##NMAIN|^##FPLOIDY' > glimpse.header.txt
        bcftools view --no-version -h ~{kage_vcf_gz} | grep '^#CHROM' > kage.columns.txt
        cat kage.header.txt glimpse.header.txt kage.columns.txt > header.txt
        bcftools reheader -h header.txt ~{output_prefix}.kage.glimpse.raw.vcf.gz > ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, command_mem_gb]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, preemptible])
        maxRetries: select_first([runtime_attributes.max_retries, 3])       # OOM can resolve due to non-determinism without memory increase, set retry factor close to 1!
    }

    output {
        File monitoring_log = "monitoring.log"
        File phased_vcf_gz = "~{output_prefix}.vcf.gz"          # note that this actually contains unphased */* GTs; these are converted to phased *|* GTs later on in the sample step
        File phased_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}

task GLIMPSELigate {
    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    Int disk_size_gb = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_ligate_static
        chmod +x GLIMPSE_ligate_static

        ./GLIMPSE_ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.vcf.gz --thread $(nproc)
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File ligated_vcf_gz = "~{prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 7]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, disk_size_gb]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 10])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task GLIMPSESample {
    input {
        File ligated_vcf_gz
        File ligated_vcf_gz_tbi
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_sample_static
        chmod +x GLIMPSE_sample_static

        ./GLIMPSE_sample_static --input ~{ligated_vcf_gz} \
            --solve \
            --output ~{output_prefix}.vcf.gz \
            --log ~{output_prefix}.sample.log
        bcftools index -t ~{output_prefix}.vcf.gz
    }

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
        File monitoring_log = "monitoring.log"
        File glimpse_vcf_gz = "~{output_prefix}.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }
}

task FixVariantCollisions {

    input {
        File vcf_gz                         # biallelic
        File vcf_gz_tbi
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix

        File? monitoring_script
    }

    command <<<
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        java ~{fix_variant_collisions_java} \
            ~{vcf_gz} \
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
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File monitoring_log = "monitoring.log"
        File collisionless_vcf_gz = "~{output_prefix}.vcf.gz"
        File collisionless_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: 10
        preemptible:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
