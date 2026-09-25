version 1.0

task call_inv_flag_insdel_cluster_indel_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz"
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

task call_inv_flag_insdel_cluster_sv_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String hap
    String threads
    String mem_gb
    String vartype
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{inbed}
    tar zxvf ~{pav_asm}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz temp/~{sample}/inv_caller/flag/insdel_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_flag_insdel_cluster_indel_~{sample}_~{vartype}_~{hap}.tgz"
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

task call_inv_cluster_indel_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_cluster_indel_~{hap}_~{sample}_~{vartype}.tgz temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_cluster_indel_~{hap}_~{sample}_~{vartype}.tgz"
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

task call_inv_cluster_snv_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    File inbed
    String threads
    String mem_gb
    String hap
    String vartype
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{inbed}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
    tar zcvf call_inv_cluster_snv_~{hap}_~{sample}_~{vartype}.tgz temp/~{sample}/inv_caller/flag/cluster_~{vartype}_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_cluster_snv_~{hap}_~{sample}_~{vartype}.tgz"
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

task call_inv_merge_flagged_loci_hap {
  input {
    String sample
    String hap
    File pav_conf
    File pav_asm
    File indelFlag
    File svFlag
    File snvCluster
    File indelCluster
    String threads
    String mem_gb
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{indelFlag}
    tar zxvf ~{svFlag}
    tar zxvf ~{snvCluster}
    tar zxvf ~{indelCluster}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} results/~{sample}/inv_caller/flagged_regions_~{hap}.bed.gz
    tar zcvf call_inv_merge_flagged_loci_~{hap}_~{sample}.tgz results/~{sample}/inv_caller/flagged_regions_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_merge_flagged_loci_~{hap}_~{sample}.tgz"
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

task call_inv_batch_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    String hap
    File trimBed
    File flag
    File asmGz
    File refGz
    String batch
    String threads
    String mem_gb
    String zones = "us-central1-c us-central1-b"
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    tar zxvf ~{trimBed}
    tar zxvf ~{flag}
    tar zxvf ~{asmGz}
    tar zxvf ~{refGz}
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/batch/~{hap}/inv_call_~{batch}.bed.gz
    tar zcvf call_inv_batch_~{hap}_~{batch}_~{sample}.tgz temp/~{sample}/inv_caller/batch/~{hap}/inv_call_~{batch}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_batch_~{hap}_~{batch}_~{sample}.tgz"
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

task call_inv_batch_merge_hap {
  input {
    String sample
    File pav_conf
    File pav_asm
    String hap
    Array[File] invBed
    String threads
    String mem_gb
  }
  command <<<
    source activate lr-pav
    set -eux
    cp ~{pav_conf} ./config.json
    tar zxvf ~{pav_asm}
    echo ~{sep=" " invBed} | tr " " "\n" | xargs -I '@' tar zxvf @
    mv /opt/pav /cromwell_root/
    tree
    snakemake -s pav/Snakefile --cores ~{threads} temp/~{sample}/inv_caller/sv_inv_~{hap}.bed.gz
    tar zcvf call_inv_batch_merge_~{hap}_~{sample}.tgz temp/~{sample}/inv_caller/sv_inv_~{hap}.bed.gz
  >>>
  output {
    Array[File] snakemake_logs = glob(".snakemake/log/*.snakemake.log")
    File bed = "call_inv_batch_merge_~{hap}_~{sample}.tgz "
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
