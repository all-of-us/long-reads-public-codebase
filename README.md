# All of Us Long Read Phase 1 Workflows

This repository contains all of the reproducible WDL workflows used in **Phase 1** of the *All of Us* Long Read (AoU-LR) project. These workflows cover various steps of long-read genomic analysis and are provided for transparency and reuse.  

The intended audience is researchers and bioinformaticians who wish to utilize or reference these workflows in their own analyses.

---

## Workflows

### Terra-Hosted Workflows
*The following workflows are available through Terra. Each links to a Terra workspace; no description is provided here.*

| **Workflow Name** | **Location** | **Description** |
|-----------------|-----------|-------------|
| CalculateCoverage | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/CalculateCoverage/2) |  |
| CleanupVCF | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/CleanupVCF/3) |  |
| ConsensusAsmFromShapeitScaffold | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/ConsensusAsmFromShapeitScaffold/6) |  |
| convert_lower_case | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/convert_lower_case/5) |  |
| kage-lite-dev-LOO-evaluation | [Terra Workflow](https://app.terra.bio/#workflows/slee-aou/kage-lite-dev-LOO-evaluation/3) |  |
| kage-lite-panel-with-preprocessing | [Terra Workflow](https://app.terra.bio/#workflows/slee-aou/kage-lite-panel-with-preprocessing/1) |  |
| MergePhasedVcf | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/MergePhasedVcf/5) |  |
| Shapeit5 | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/Shapeit5/5) |  |
| StatisticalPhasing | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/StatisticalPhasing/5) |  |
| vcfdist | [Terra Workflow](https://app.terra.bio/#workflows/slee-aou/vcfdist/2) |  |
| WhatshapPhasing | [Terra Workflow](https://app.terra.bio/#workflows/hangsu/WhatshapPhasing/24) |  |

---

### Dockstore Workflows (Organized by Functionality)

#### 1. Read and Assembly Processing
| Workflow | Location | Description |
|----------|----------|-------------|
| HiFiBamToFastQ | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/HiFiBamToFastQ:sh_wg_extract_fastq) | Converts PacBio HiFi BAM files to FASTQ format for downstream analysis. |
| MergeFastqs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/MergeFastqs:sh_hificnv) | Merges multiple FASTQ files into a single set per sample. |
| PBAssembleWithHifiasm | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/PBAssembleWithHifiasm:sh_update_hifiasm) | Assembles PacBio HiFi reads using Hifiasm. |
| MapAssemblyContigs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/MapAssemblyContigs:sh_hificnv) | Aligns assembly contigs to a reference to validate assemblies or identify SVs. |
| EvaluateAssemblyHap | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/EvaluateAssembly:kvg_summarize_haplotype_fastas) | Assesses haplotype-resolved assemblies against reference sequences. |

---

#### 2. Small Variant Calling & Summaries
| Workflow | Location | Description |
|----------|----------|-------------|
| 1074.T2T.SmallVariantsBasicMetrics | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SmallVariantsBasicMetrics:sh_variants_ss_metrics) | Computes metrics for small variants on T2T reference. |
| PBCCSWholeGenome | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/PBCCSWholeGenome:main) | Calls small variants genome-wide from PacBio CCS reads. |
| SummarizeDVPSmallVariants | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SummarizeSmallVariants:kvg_summarize_vcf) | Summarizes variants from DeepVariant+Pepper pipeline. |
| SummarizePAVSmallVariants | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SummarizeSmallVariants:kvg_summarize_vcf) | Summarizes small variants discovered in PAV contexts. |

---

#### 3. Structural Variant (SV) Discovery & Integration
| Workflow | Location | Description |
|----------|----------|-------------|
| PAV | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/pav-wdl/pav:sh_more_resources_pete) | Detects presence/absence variants (large insertions/deletions). |
| PAV2SVs | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/callset_integration/PAV2SVs:main) | Converts PAV results to standard SV calls. |
| LRMergeSVVCFs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/LRMergeSVVCFs:kvg_sveval) | Merges multiple SV callsets into one. |
| TruvariCollapse | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/callset_integration/TruvariCollapse:main) | Collapses duplicate/equivalent SVs into consensus calls. |
| TruvariIntersample | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/callset_integration/TruvariIntersample:main) | Compares SVs between samples. |
| TruvariIntrasample | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/TruvariIntrasample:main) | Compares SVs within a single sample. |
| SummarizePAVSVs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SummarizeStructuralVariants:sh_baseoff_kvg_summarize_vcf) | Summarizes PAV structural variant calls. |
| SummarizeSnifflesSVs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SummarizeStructuralVariants:sh_baseoff_kvg_summarize_vcf) | Summarizes SVs called by Sniffles. |
| GraphEvaluation | [Dockstore](https://dockstore.org/workflows/github.com/rlorigro/sv_merge/GraphEvaluation:dev) | Builds and evaluates SV overlap graphs. |

---

#### 4. Joint Calling & Cohort Integration
| Workflow | Location | Description |
|----------|----------|-------------|
| JointCalling | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/JointCalling:main) | Jointly genotypes SVs across a cohort. |
| LRJointCallGVCFs | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/LRJointCallGVCFs:main) | Joint genotyping of GVCFs into a cohort-wide VCF. |
| MergeVCFs | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/MergeVCFs:main) | Merges multiple VCFs into one. |
| MergePhasedVCF | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/MergePhasedVCF:hangsu_phasing) | Merges phased VCFs into one. |
| MergeRegenotypedIntersampleVcf | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/callset_integration/MergeRegenotypedIntersampleVcf:main) | Merges per-sample re-genotyped VCFs into cohort VCF. |
| MergeSVsSNPs | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/callset_integration/MergeSVsSNPs:main) | Combines SVs with SNVs/indels into one file. |
| OverlapGraph | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/OverlapGraph:main) | Builds an overlap graph across callsets. |
| OverlapStats | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/OverlapStats:main) | Computes overlap statistics between callsets. |

---

#### 5. Phasing
| Workflow | Location | Description |
|----------|----------|-------------|
| Shapeit4PhaseWholeGenome | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/Shapeit4PhaseWholeGenome:hangsu_phasing_new) | Phases variants genome-wide with Shapeit4. |
| HybridPhaseWholeGenome | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/HybridPhaseWholeGenome:hangsu_phasing_new) | Whole-genome hybrid phasing with long-reads and statistical methods. |
| HybridPhaseWholeGenomeHiphase | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/HybridPhaseWholeGenomeHiphase:hangsu_phasing_new) | Hybrid whole-genome phasing incorporating HIPHASE. |
| HiphaseJointCall | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/HiphaseJointCall:hangsu_phasing_new) | Joint calling with HIPHASE phasing. |
| HiphaseJointCallTRGT | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/HiphaseJointCallTRGT:hangsu_phasing_new) | HIPHASE joint calling restricted to target regions. |
| Flare | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/Flare:hangsu_phasing_new) | Long-read phasing using FLARE. |
| PhasedMerge | [Dockstore](https://dockstore.org/workflows/github.com/fabio-cunial/sv-merging/PhasedMerge:main) | Merges callsets while preserving phasing. |
| PhasingStatistics | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/PhasingStatistics:hangsu_phasing_new) | Summarizes phasing performance and block sizes. |

---

#### 6. Quality Control & Fingerprinting
| Workflow | Location | Description |
|----------|----------|-------------|
| CollectSingleSampleSVvcfMetrics | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/CollectSingleSampleSVvcfMetrics:sh_gatk_sv_vcfmetrics) | Computes SV metrics per sample. |
| LongReadsContaminationEstimation | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/LongReadsContaminationEstimation:sh_more_atomic_qc) | Estimates contamination in long-read data. |
| BuildTempLocalFpStore | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/BuildTempLocalFpStore:sh_exhaust_au) | Builds temporary fingerprint store for identity checks. |
| VerifyFingerprintCCSSample | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/VerifyFingerprintCCSSample:sh_sample_fp) | Verifies CCS sample identity by fingerprinting. |
| SexCheck | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/long-read-pipelines/SexCheck:sh_more_atomic_qc) | Checks reported vs genetic sex. |
| MainVcfQc | [Dockstore](https://dockstore.org/workflows/github.com/broadinstitute/gatk-sv/17-MainVcfQc:v0.26.2-beta) | Runs quality control checks on final VCFs. |
