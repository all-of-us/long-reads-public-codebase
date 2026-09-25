version 1.0

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

workflow HierarchicallyMergeVcfs {
    input {
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        Array[String] regions   # bcftools regions, e.g. ["chr1,chr2,chr3", "chr4,chr5,chr6", ...]
        Int batch_size
        String output_prefix
        String extra_merge_args = "--threads $(nproc) --force-single --merge none"       # non-region args
        String extra_concat_args = "--threads $(nproc) --naive"
        Boolean use_ivcfmerge   # requires sample_names
        Array[String]? sample_names

        RuntimeAttributes merge_runtime_attributes = {}
        RuntimeAttributes concat_runtime_attributes = {}

        String docker
        File? monitoring_script
    }

    call CreateBatches {
        input:
            vcf_gzs = vcf_gzs,
            vcf_gz_tbis = vcf_gz_tbis,
            batch_size = batch_size,
            docker = docker
    }

    if (use_ivcfmerge) {
        scatter (i in range(length(CreateBatches.vcf_gz_batch_fofns))) {
            scatter (j in range(length(regions))) {
                call Ivcfmerge as IvcfmergeSingleBatchRegion {
                    input:
                        vcf_gzs = read_lines(CreateBatches.vcf_gz_batch_fofns[i]),
                        vcf_gz_tbis = read_lines(CreateBatches.vcf_gz_tbi_batch_fofns[i]),
                        output_prefix = output_prefix + ".batch-" + i + ".region-" + i,
                        sample_names = select_first([sample_names]),
                        region_args = "-r " + regions[j],
                        docker = docker,
                        monitoring_script = monitoring_script,
                        runtime_attributes = merge_runtime_attributes
                }
            }
        }
    }
    if (!use_ivcfmerge) {
        scatter (i in range(length(CreateBatches.vcf_gz_batch_fofns))) {
            scatter (j in range(length(regions))) {
                call MergeVcfs as MergeVcfsSingleBatchRegion {
                    input:
                        vcf_gzs = read_lines(CreateBatches.vcf_gz_batch_fofns[i]),
                        vcf_gz_tbis = read_lines(CreateBatches.vcf_gz_tbi_batch_fofns[i]),
                        output_prefix = output_prefix + ".batch-" + i + ".region-" + i,
                        extra_args = "-r " + regions[j] + " " + extra_merge_args,
                        docker = docker,
                        monitoring_script = monitoring_script,
                        runtime_attributes = merge_runtime_attributes
                }
            }
        }
    }

    Array[Array[File]] region_by_batch_vcf_gzs = transpose(select_first([IvcfmergeSingleBatchRegion.merged_vcf_gz, MergeVcfsSingleBatchRegion.merged_vcf_gz]))
    Array[Array[File]] region_by_batch_vcf_gz_tbis = transpose(select_first([IvcfmergeSingleBatchRegion.merged_vcf_gz_tbi, MergeVcfsSingleBatchRegion.merged_vcf_gz_tbi]))

    if (use_ivcfmerge) {
        # merge all samples in each region
        scatter (j in range(length(regions))) {
            call Ivcfmerge as IvcfmergeSingleRegion {
                input:
                    vcf_gzs = region_by_batch_vcf_gzs[j],
                    vcf_gz_tbis = region_by_batch_vcf_gz_tbis[j],
                    output_prefix = output_prefix + ".region-" + j,
                    sample_names = select_first([sample_names]),
                    docker = docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = merge_runtime_attributes
            }
        }
    }
    if (!use_ivcfmerge) {
        # merge all samples in each region
        scatter (j in range(length(regions))) {
            call MergeVcfs as MergeVcfsSingleRegion {
                input:
                    vcf_gzs = region_by_batch_vcf_gzs[j],
                    vcf_gz_tbis = region_by_batch_vcf_gz_tbis[j],
                    output_prefix = output_prefix + ".region-" + j,
                    extra_args = extra_merge_args,
                    docker = docker,
                    monitoring_script = monitoring_script,
                    runtime_attributes = merge_runtime_attributes
            }
        }
    }

    # concatenate all regions
    call ConcatVcfs {
        input:
            vcf_gzs = select_first([IvcfmergeSingleRegion.merged_vcf_gz, MergeVcfsSingleRegion.merged_vcf_gz]),
            vcf_gz_tbis = select_first([IvcfmergeSingleRegion.merged_vcf_gz_tbi, MergeVcfsSingleRegion.merged_vcf_gz_tbi]),
            output_prefix = output_prefix,
            extra_args = extra_concat_args,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = concat_runtime_attributes
    }

    output {
        File merged_vcf_gz = ConcatVcfs.concatenated_vcf_gz
        File merged_vcf_gz_tbi = ConcatVcfs.concatenated_vcf_gz_tbi
    }
}

task CreateBatches {
    input {
        Array[String] vcf_gzs
        Array[String] vcf_gz_tbis
        Int batch_size

        String docker
        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euox pipefail

        cat ~{write_lines(vcf_gzs)} | split -l ~{batch_size} - vcf_gz_batch_
        cat ~{write_lines(vcf_gz_tbis)} | split -l ~{batch_size} - vcf_gz_tbi_batch_
    }

    output {
        Array[File] vcf_gz_batch_fofns = glob("vcf_gz_batch_*")
        Array[File] vcf_gz_tbi_batch_fofns = glob("vcf_gz_tbi_batch_*")
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
}

task MergeVcfs {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix
        String? extra_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools merge \
            -l ~{write_lines(vcf_gzs)} \
            ~{extra_args} \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    output {
        File monitoring_log = "monitoring.log"
        File merged_vcf_gz = "~{output_prefix}.vcf.gz"
        File merged_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
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
}

# assumes all VCFs have identical variants
# TODO make this work for identical filenames
task Ivcfmerge {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        Array[String] sample_names
        String output_prefix
        String? region_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        wget https://github.com/iqbal-lab-org/ivcfmerge/archive/refs/tags/v1.0.0.tar.gz
        tar -xvf v1.0.0.tar.gz

        mkdir compressed
        mv ~{sep=' ' vcf_gzs} compressed
        mv ~{sep=' ' vcf_gz_tbis} compressed

        if [ $(ls compressed/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls compressed/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls compressed/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            mkdir decompressed
            ls compressed/*.vcf.gz | xargs -I % sh -c 'bcftools annotate --no-version ~{region_args} -x INFO % --threads 2 -Ov -o decompressed/$(basename % .gz)'
            time python ivcfmerge-1.0.0/ivcfmerge.py <(ls decompressed/*.vcf) ~{output_prefix}.vcf
            bcftools annotate --no-version -S ~{write_lines(sample_names)} -x FORMAT/FT ~{output_prefix}.vcf --threads 2 -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
    }

    output {
        File monitoring_log = "monitoring.log"
        File merged_vcf_gz = "~{output_prefix}.vcf.gz"
        File merged_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 4])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 250]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }
}

task ConcatVcfs {
    input{
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix
        String? extra_args

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {"use_ssd": true}
    }

    command {
        set -euox pipefail

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools concat \
            -f ~{write_lines(vcf_gzs)} \
            ~{extra_args} \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t --threads $(nproc) ~{output_prefix}.vcf.gz
    }

    output {
        File monitoring_log = "monitoring.log"
        File concatenated_vcf_gz = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
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
}
