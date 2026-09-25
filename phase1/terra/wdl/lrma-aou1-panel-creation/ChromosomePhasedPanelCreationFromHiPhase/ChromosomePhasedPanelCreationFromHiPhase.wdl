########################################################################################
# ChromosomePhasedPanelCreationFromHiPhase (the WDL names its workflow PhasedPanelEvaluation)
# Phase 1 stage c (reference panel), launch step 7.
#
# Builds the reference panel for one chromosome from the physically phased callsets:
# filters and combines the SNV/indel and SV VCFs, resolves overlapping ("colliding")
# alleles with fix_variant_collisions_java, phases statistically with SHAPEIT4 in
# chunks that are then ligated, resolves collisions again, and creates the PanGenie
# bubble-graph panel VCF.
#
# Runs once per chromosome.
# Inputs:  hiphase_short_vcf_gz and hiphase_sv_vcf_gz (cohort VCFs built from the
#          PhysicalPhasing outputs), chromosome, reference_fasta,
#          genetic_mapping_tsv_for_shapeit4, fix_variant_collisions_java, operation /
#          weight_tag (collision-resolution settings), helper scripts, filter arguments
# Outputs: phased_vcf_gz     the statistically phased chromosome VCF
#          collisionless_bcf the phased VCF with collisions resolved
#          panel_vcf_gz      the PanGenie panel VCF
# Previous: PhysicalPhasing (all samples). Next: ConcatAndEvaluateReduced.
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0

import "repo/wdl/methods/pangenie/PanGeniePanelCreation.wdl"
import "repo/wdl/methods/phasing/Helper.wdl"

workflow PhasedPanelEvaluation {    # TODO change name later, easier to share configs for now

    input {
        # common inputs
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String chromosome
        String output_prefix
        File? monitoring_script

        File hiphase_short_vcf_gz
        File hiphase_short_vcf_gz_tbi
        File hiphase_sv_vcf_gz
        File hiphase_sv_vcf_gz_tbi
        Boolean subset_short_to_sv_windows
        Int window_padding
        String? subset_filter_args
        String? filter_and_concat_short_filter_args
        String? filter_and_concat_sv_filter_args
        String? extra_chunk_args
        File genetic_mapping_tsv_for_shapeit4
        Int shapeit4_num_threads
        Int shapeit4_memory
        String shapeit4_extra_args

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int operation
        String weight_tag
        Int is_weight_format_field

        # inputs for PanGeniePanelCreation
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing
        String panel_creation_docker
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    if (subset_short_to_sv_windows) {
        call SubsetVcfShortInSVWindows as SubsetVcfShortInSVWindowsChromosome { input:
            short_vcf_gz = hiphase_short_vcf_gz,
            short_vcf_tbi = hiphase_short_vcf_gz_tbi,
            sv_vcf_gz = hiphase_sv_vcf_gz,
            sv_vcf_tbi = hiphase_sv_vcf_gz_tbi,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            reference_fasta_dict = reference_fasta_dict,
            prefix = output_prefix + "." + chromosome + ".short.subset.windowed",
            region = chromosome,
            window_padding = window_padding,
            filter_args = subset_filter_args
        }
    }

    call FilterAndConcatVcfs as FilterAndConcatVcfsChromosome { input:
        short_vcf = select_first([SubsetVcfShortInSVWindowsChromosome.subset_short_vcf_gz, hiphase_short_vcf_gz]),
        short_vcf_tbi = select_first([SubsetVcfShortInSVWindowsChromosome.subset_short_vcf_gz_tbi, hiphase_short_vcf_gz_tbi]),
        sv_vcf = hiphase_sv_vcf_gz,
        sv_vcf_tbi = hiphase_sv_vcf_gz_tbi,
        region = chromosome,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        filter_and_concat_short_filter_args = filter_and_concat_short_filter_args,
        filter_and_concat_sv_filter_args = filter_and_concat_sv_filter_args,
        prefix = output_prefix + "." + chromosome + ".filter_and_concat"
    }

    call FixVariantCollisions as BeforeShapeit4FixVariantCollisionsChromosome { input:
        phased_bcf = FilterAndConcatVcfsChromosome.filter_and_concat_vcf,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix + "." + chromosome + ".before_shapeit4_collisionless"
    }

    call CreateShapeit4Chunks { input:
        vcf = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf,
        tbi = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf_csi,
        region = chromosome,
        prefix = output_prefix + "." + chromosome,
        extra_chunk_args = extra_chunk_args
    }

    Array[String] region_list = read_lines(CreateShapeit4Chunks.chunks)

    scatter (j in range(length(region_list))) {
        call Helper.Shapeit4 as Shapeit4 { input:
            vcf_input = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf,
            vcf_index = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf_csi,
            mappingfile = genetic_mapping_dict[chromosome],
            region = region_list[j],
            prefix = output_prefix + "." + chromosome + ".shard-" + j + ".filter_and_concat.phased",
            num_threads = shapeit4_num_threads,
            memory = shapeit4_memory,
            extra_args = shapeit4_extra_args
        }
    }

    call LigateVcfs { input:
        vcfs = Shapeit4.phased_bcf,
        prefix = output_prefix + "." + chromosome + ".phased.ligated"
    }

    call FixVariantCollisions as FixVariantCollisionsChromosome { input:
        phased_bcf = LigateVcfs.ligated_vcf_gz,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix + "." + chromosome + ".phased.collisionless"
    }

    call PanGeniePanelCreation.PanGeniePanelCreation as PanGeniePanelCreationChromosome { input:
        phased_bcf = FixVariantCollisionsChromosome.collisionless_bcf,
        reference_fasta = reference_fasta,
        prepare_vcf_script = prepare_vcf_script,
        add_ids_script = add_ids_script,
        merge_vcfs_script = merge_vcfs_script,
        frac_missing = frac_missing,
        output_prefix = output_prefix + "." + chromosome,
        docker = panel_creation_docker,
        monitoring_script = monitoring_script
    }

    output {
        File filter_and_concat_vcf_gz = FilterAndConcatVcfsChromosome.filter_and_concat_vcf
        File filter_and_concat_vcf_gz_tbi = FilterAndConcatVcfsChromosome.filter_and_concat_vcf_tbi

        File before_shapeit4_collisionless_bcf = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf
        File before_shapeit4_collisionless_bcf_csi = BeforeShapeit4FixVariantCollisionsChromosome.collisionless_bcf_csi

        File phased_vcf_gz = LigateVcfs.ligated_vcf_gz
        File phased_vcf_gz_tbi = LigateVcfs.ligated_vcf_gz_tbi

        File collisionless_bcf = FixVariantCollisionsChromosome.collisionless_bcf
        File collisionless_bcf_csi = FixVariantCollisionsChromosome.collisionless_bcf_csi

        File panel_vcf_gz = PanGeniePanelCreationChromosome.panel_vcf_gz
        File panel_vcf_gz_tbi = PanGeniePanelCreationChromosome.panel_vcf_gz_tbi
    }
}

task SubsetVcfShortInSVWindows {
    input {
        File short_vcf_gz
        File short_vcf_tbi
        File sv_vcf_gz
        File sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File reference_fasta_dict
        String prefix
        String region
        Int window_padding
        String? filter_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([short_vcf_gz, short_vcf_tbi], "GB")) + 2*ceil(size([sv_vcf_gz, sv_vcf_tbi], "GB")) + 1

    command {
        set -euxo pipefail

        bcftools view ~{sv_vcf_gz} \
            -r ~{region} \
            -Oz -o sv.region.vcf.gz
        bcftools index -t sv.region.vcf.gz
        gatk PreprocessIntervals \
            -L sv.region.vcf.gz \
            --reference ~{reference_fasta} \
            --padding ~{window_padding} \
            --bin-length 0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{prefix}.windows.interval_list
        gatk IntervalListToBed \
            -I ~{prefix}.windows.interval_list \
            -O ~{prefix}.windows.bed
        bcftools view ~{short_vcf_gz} \
            -R ~{prefix}.windows.bed \
            ~{filter_args} \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    }

    output {
        File subset_short_vcf_gz = "~{prefix}.vcf.gz"
        File subset_short_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# filter out singletons (i.e., keep MAC >= 2) and concatenate with deduplication
task FilterAndConcatVcfs {

    input {
        File short_vcf         # multiallelic
        File short_vcf_tbi
        File sv_vcf            # biallelic
        File sv_vcf_tbi
        String prefix
        String region
        File reference_fasta
        File reference_fasta_fai
        String? filter_and_concat_short_filter_args = "-i 'MAC>=2'"
        String? filter_and_concat_sv_filter_args = "-i 'MAC>=2'"

        RuntimeAttr? runtime_attr_override
    }

    command {
        set -euxo pipefail

        # filter SV singletons
        bcftools view ~{filter_and_concat_sv_filter_args} ~{sv_vcf} \
            -r ~{region} \
            --write-index -Oz -o ~{prefix}.SV.vcf.gz

        # filter short singletons and split to biallelic
        bcftools view ~{filter_and_concat_short_filter_args} ~{short_vcf} \
            -r ~{region} | \
            bcftools norm -m-any -f ~{reference_fasta} | \
            bcftools sort \
            --write-index -Oz -o ~{prefix}.short.vcf.gz

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat \
            ~{prefix}.SV.vcf.gz \
            ~{prefix}.short.vcf.gz \
            --allow-overlaps --remove-duplicates \
            -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    }

    output {
        File filter_and_concat_vcf = "~{prefix}.vcf.gz"
        File filter_and_concat_vcf_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateShapeit4Chunks {

    input {
        File vcf
        File tbi
        String region
        String prefix
        String? extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 1

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
        cut -f 3 chunks.txt > chunks.regions.txt
    >>>

    output {
        File chunks = "chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task LigateVcfs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        wget https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/ligate_static
        chmod +x ligate_static

        ./ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf_gz = "~{prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FixVariantCollisions {

    input {
        File phased_bcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix
    }

    command <<<
        set -euxo pipefail

        # convert bcf to vcf.gz
        bcftools view ~{phased_bcf} -Oz -o phased.vcf.gz
        bcftools index -t phased.vcf.gz

        java ~{fix_variant_collisions_java} \
            phased.vcf.gz \
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
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.vcf.gz
        bcftools view ~{output_prefix}.vcf.gz -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File collisionless_bcf = "~{output_prefix}.bcf"
        File collisionless_bcf_csi = "~{output_prefix}.bcf.csi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 100 SSD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
