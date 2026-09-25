version 1.0

task call_cigar_hap {
  input {
    File pav_conf
    File pav_asm
    File trimBed
    File asmGz
    File refGz
    String batch
    String threads
    String sample
    String mem_gb
    String hap
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{trimBed}
    tar zxvf ~{asmGz}
    tar zxvf ~{refGz}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/cigar/batched/insdel_~{hap}_~{batch}.bed.gz temp/~{sample}/cigar/batched/snv.bed_~{hap}_~{batch}.gz
    tar zcvf call_cigar_~{hap}_~{sample}_~{batch}.tgz temp/~{sample}/cigar/batched/insdel_~{hap}_~{batch}.bed.gz temp/~{sample}/cigar/batched/snv.bed_~{hap}_~{batch}.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File snvBed = "call_cigar_~{hap}_~{sample}_~{batch}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
  }
}

task call_cigar_merge_hap {
  input {
    File pav_conf
    File pav_asm
    String sample
    String hap
    Array[File] snvBatch
    String threads
    String mem_gb
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    echo ~{sep=" " snvBatch} | tr " " "\n" | xargs -I '@' tar zxvf @
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/cigar/pre_inv/svindel_insdel_~{hap}.bed.gz temp/~{sample}/cigar/pre_inv/snv_snv_~{hap}.bed.gz
    tar zcvf call_cigar_merge_~{hap}_~{sample}.tgz temp/~{sample}/cigar/pre_inv/svindel_insdel_~{hap}.bed.gz temp/~{sample}/cigar/pre_inv/snv_snv_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File insdelBedMerge = "call_cigar_merge_~{hap}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
  }
}

task call_mappable_bed_hap {
  input {
    File pav_conf
    File pav_asm
    String sample
    File delBed
    File insBed
    File invBed
    File trimBed
    String hap
    String threads
    String mem_gb
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{delBed}
    tar zxvf ~{insBed}
    tar zxvf ~{invBed}
    tar zxvf ~{trimBed}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} results/~{sample}/callable/callable_regions_~{hap}_500.bed.gz
    tar zcvf call_mappable_bed_~{hap}_~{sample}.tgz results/~{sample}/callable/callable_regions_~{hap}_500.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_mappable_bed_~{hap}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
  }
}

task call_integrate_sources_hap {
  input {
    File pav_conf
    File pav_asm
    String sample
    String hap
    File preInvSvBed
    File delBedIn
    File insBedIn
    File invBedIn
    File invBatch
    String threads
    String mem_gb
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{preInvSvBed}
    tar zxvf ~{delBedIn}
    tar zxvf ~{insBedIn}
    tar zxvf ~{invBedIn}
    tar zxvf ~{invBatch}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/bed/integrated/~{hap}/svindel_ins.bed.gz temp/~{sample}/bed/integrated/~{hap}/svindel_del.bed.gz temp/~{sample}/bed/integrated/~{hap}/snv_snv.bed.gz temp/~{sample}/bed/integrated/~{hap}/sv_inv.bed.gz
    tar zcvf call_integrate_sources_~{hap}_~{sample}.tgz temp/~{sample}/bed/integrated/~{hap}/svindel_ins.bed.gz temp/~{sample}/bed/integrated/~{hap}/svindel_del.bed.gz temp/~{sample}/bed/integrated/~{hap}/snv_snv.bed.gz temp/~{sample}/bed/integrated/~{hap}/sv_inv.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File all_vars_bed = "call_integrate_sources_~{hap}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
  }
}

task call_merge_haplotypes_chrom_svindel {
  input {
    String svtype
    File pav_conf
    File pav_asm
    String sample
    File svindel_bed_h1
    File svindel_bed_h2
    File callable_h1
    File callable_h2
    String chrom
    String threads
    String mem_gb
    String zones = "us-central1-c us-central1-b"
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{svindel_bed_h1}
    tar zxvf ~{svindel_bed_h2}
    tar zxvf ~{callable_h1}
    tar zxvf ~{callable_h2}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/bed/bychrom/~{svtype}/~{chrom}.bed.gz
    tar zcvf call_merge_haplotypes_chrom_svindel_~{sample}_~{svtype}_~{chrom}.tgz temp/~{sample}/bed/bychrom/~{svtype}/~{chrom}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_merge_haplotypes_chrom_svindel_~{sample}_~{svtype}_~{chrom}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
      zones: zones
  }
}

task call_merge_haplotypes_chrom_snv {
  input {
    String svtype
    File snvBed_h1
    File pav_conf
    File pav_asm
    File snvBed_h2
    File callable_h1
    String sample
    File callable_h2
    String chrom
    String threads
    String mem_gb
    String zones = "us-central1-c us-central1-b"
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{snvBed_h1}
    tar zxvf ~{snvBed_h2}
    tar zxvf ~{callable_h1}
    tar zxvf ~{callable_h2}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/bed/bychrom/~{svtype}/~{chrom}.bed.gz
    tar zcvf call_merge_haplotypes_chrom_snv_~{sample}_~{svtype}_~{chrom}.tgz temp/~{sample}/bed/bychrom/~{svtype}/~{chrom}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_merge_haplotypes_chrom_snv_~{sample}_~{svtype}_~{chrom}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk " + 100 + " HDD"
      bootDiskSizeGb: 50
      preemptible:    3
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
      zones: zones
  }
}

task call_merge_haplotypes {
  input {
    String svtype
    Array[File] inbed
    File pav_conf
    File pav_asm
    File integrated_h1
    File integrated_h2
    File callable_h1
    File callable_h2
    String sample
    String threads
    String mem_gb
    String zones = "us-central1-c us-central1-b"
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{integrated_h1}
    tar zxvf ~{integrated_h2}
    tar zxvf ~{callable_h1}
    tar zxvf ~{callable_h2}
    echo ~{sep=" " inbed} | tr " " "\n" | xargs -I '@' tar zxvf @
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/bed/merged/~{svtype}.bed.gz
    tar zcvf call_merge_haplotypes_~{svtype}_~{sample}.tgz temp/~{sample}/bed/merged/~{svtype}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_merge_haplotypes_~{svtype}_~{sample}.tgz"
  }
  ############################
  runtime {
      cpu:            threads
      memory:         mem_gb + " GiB"
      disks:          "local-disk 100 HDD"
      bootDiskSizeGb: 50
      preemptible:    0
      maxRetries:     1
      docker:         "us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou"
      zones: zones
  }
}
