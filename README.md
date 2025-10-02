# All of Us Long Read Phase 1 Workflows

This repository contains all of the reproducible WDL workflows used in **Phase 1** of the *All of Us* Long Read (AoU-LR) project. These workflows cover various steps of long-read genomic analysis and are provided for transparency and reuse.  

Please note: this code's organization is in flux.

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

---

## Jupyter Notebooks

This repository also contains Jupyter notebooks for data analysis and visualization, organized by platform:

### Terra Notebooks (`notebooks/terra/`)

These notebooks are designed to run in the Terra cloud platform and focus on data processing, analysis, and quality control:

#### Data Import and Processing
| Notebook | Link | Description |
|----------|------|-------------|
| 20231206_Senegal_2020_Data_Import_And_Processing.ipynb | [GitHub](notebooks/terra/20231206_Senegal_2020_Data_Import_And_Processing.ipynb) | Data import and processing for Senegal 2020 dataset |
| main_init_subset_vds.ipynb | [GitHub](notebooks/terra/main_init_subset_vds.ipynb) | Initialize and subset Variant Dataset (VDS) for analysis |

#### Assembly Analysis
| Notebook | Link | Description |
|----------|------|-------------|
| kvg_examine_assemblies.ipynb | [GitHub](notebooks/terra/kvg_examine_assemblies.ipynb) | Examine and analyze genome assemblies |
| kvg_study_read_length_dists.ipynb | [GitHub](notebooks/terra/kvg_study_read_length_dists.ipynb) | Study read length distributions from sequencing data |

#### Variant Analysis
| Notebook | Link | Description |
|----------|------|-------------|
| kvg_examine_small_variants.ipynb | [GitHub](notebooks/terra/kvg_examine_small_variants.ipynb) | Analyze small variants (SNPs, indels) |
| kvg_examine_structural_variants.ipynb | [GitHub](notebooks/terra/kvg_examine_structural_variants.ipynb) | Examine structural variants (SVs) |
| kvg_sv_callset_inventory.ipynb | [GitHub](notebooks/terra/kvg_sv_callset_inventory.ipynb) | Inventory and catalog structural variant callsets |
| kvg_describe_hail_matrix_tables.ipynb | [GitHub](notebooks/terra/kvg_describe_hail_matrix_tables.ipynb) | Describe Hail matrix tables for genomic data |

#### Population Genetics and Statistics
| Notebook | Link | Description |
|----------|------|-------------|
| kvg_pca.ipynb | [GitHub](notebooks/terra/kvg_pca.ipynb) | Principal Component Analysis for population structure |
| kvg_pca_hgdp_tgp.ipynb | [GitHub](notebooks/terra/kvg_pca_hgdp_tgp.ipynb) | PCA analysis incorporating HGDP and TGP reference populations |
| kvg_recompute_relatedness.ipynb | [GitHub](notebooks/terra/kvg_recompute_relatedness.ipynb) | Recompute relatedness estimates between samples |
| kvg_compute_sfs_grch38.ipynb | [GitHub](notebooks/terra/kvg_compute_sfs_grch38.ipynb) | Compute Site Frequency Spectrum on GRCh38 reference |
| kvg_firth_logistic_regression.ipynb | [GitHub](notebooks/terra/kvg_firth_logistic_regression.ipynb) | Firth logistic regression analysis |

#### Phasing and Panel Analysis
| Notebook | Link | Description |
|----------|------|-------------|
| hangsu_hiphase_results.ipynb | [GitHub](notebooks/terra/hangsu_hiphase_results.ipynb) | Analysis of HIPHASE phasing results |
| slee-aou1-analysis-for-phasing-experiments.ipynb | [GitHub](notebooks/terra/slee-aou1-analysis-for-phasing-experiments.ipynb) | Analysis for phasing experiments |
| slee-aou1-kage-panel-preprocessing-for-phasing-experiments.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-for-phasing-experiments.ipynb) | KAGE panel preprocessing for phasing experiments |
| slee-aou1-kage-panel-preprocessing-for-phasing-experiments-short-plus-sv.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-for-phasing-experiments-short-plus-sv.ipynb) | KAGE panel preprocessing with short reads and SVs |
| slee-aou1-kage-panel-preprocessing-no-TR-5Mbp.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-no-TR-5Mbp.ipynb) | KAGE panel preprocessing excluding TR regions >5Mbp |
| slee-aou1-kage-panel-preprocessing-short-only.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-short-only.ipynb) | KAGE panel preprocessing for short reads only |
| slee-aou1-kage-panel-preprocessing-SV-only-no-TR-5Mbp.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-SV-only-no-TR-5Mbp.ipynb) | KAGE panel preprocessing for SVs only, excluding TR regions >5Mbp |
| slee-aou1-kage-panel-preprocessing-sv-only.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing-sv-only.ipynb) | KAGE panel preprocessing for SVs only |
| slee-aou1-kage-panel-preprocessing.ipynb | [GitHub](notebooks/terra/slee-aou1-kage-panel-preprocessing.ipynb) | General KAGE panel preprocessing |
| slee-aou1-preprocess-for-phasing-experiments.ipynb | [GitHub](notebooks/terra/slee-aou1-preprocess-for-phasing-experiments.ipynb) | Preprocessing for phasing experiments |
| slee-aou1-roh.ipynb | [GitHub](notebooks/terra/slee-aou1-roh.ipynb) | Runs of Homozygosity (ROH) analysis |
| slee-aou1-panel-summary-statistics.ipynb | [GitHub](notebooks/terra/slee-aou1-panel-summary-statistics.ipynb) | Panel summary statistics |
| slee-aou1-panel-summary-statistics-SV.ipynb | [GitHub](notebooks/terra/slee-aou1-panel-summary-statistics-SV.ipynb) | Panel summary statistics for structural variants |

#### Quality Control
| Notebook | Link | Description |
|----------|------|-------------|
| ym_callset_QC_py.ipynb | [GitHub](notebooks/terra/ym_callset_QC_py.ipynb) | Python-based callset quality control |
| ym_callset_QC_R.ipynb | [GitHub](notebooks/terra/ym_callset_QC_R.ipynb) | R-based callset quality control |

#### Manuscript Figures and Tables
| Notebook | Link | Description |
|----------|------|-------------|
| main_figure_01_pca.ipynb | [GitHub](notebooks/terra/main_figure_01_pca.ipynb) | Generate PCA figure for main manuscript |
| main_table_02_sv_summary.ipynb | [GitHub](notebooks/terra/main_table_02_sv_summary.ipynb) | Generate structural variant summary table |
| main_table_02_variant_inventory.ipynb | [GitHub](notebooks/terra/main_table_02_variant_inventory.ipynb) | Generate variant inventory table |

### Researcher Workbench Notebooks (`notebooks/rw/`)

These notebooks are designed to run in the All of Us Researcher Workbench and focus on manuscript figures, tables, and specialized analyses:

#### Manuscript Figures
| Notebook | Link | Description |
|----------|------|-------------|
| main_figure_01_length_distributions.ipynb | [GitHub](notebooks/rw/main_figure_01_length_distributions.ipynb) | Generate read and contig length distribution figures |
| main_figure_01_map.ipynb | [GitHub](notebooks/rw/main_figure_01_map.ipynb) | Generate map figure for manuscript |
| main_figure_01_omop.ipynb | [GitHub](notebooks/rw/main_figure_01_omop.ipynb) | Generate OMOP-related figure |
| main_figure_01_pca.ipynb | [GitHub](notebooks/rw/main_figure_01_pca.ipynb) | Generate PCA figure for main manuscript |
| supp_figure_01_assembly.ipynb | [GitHub](notebooks/rw/supp_figure_01_assembly.ipynb) | Generate supplementary assembly figure |

#### Manuscript Tables
| Notebook | Link | Description |
|----------|------|-------------|
| main_table_01_dataset_summary.ipynb | [GitHub](notebooks/rw/main_table_01_dataset_summary.ipynb) | Generate dataset summary table |
| main_table_02_short_read_svs.ipynb | [GitHub](notebooks/rw/main_table_02_short_read_svs.ipynb) | Generate short read structural variant table |
| main_table_02_variant_inventory.ipynb | [GitHub](notebooks/rw/main_table_02_variant_inventory.ipynb) | Generate variant inventory table |

#### Specialized Analyses
| Notebook | Link | Description |
|----------|------|-------------|
| init_subset_vds.ipynb | [GitHub](notebooks/rw/init_subset_vds.ipynb) | Initialize and subset Variant Dataset |
| JW_CYP2D6.ipynb | [GitHub](notebooks/rw/JW_CYP2D6.ipynb) | CYP2D6 gene analysis |
| JW_repeat_expansion_figures.ipynb | [GitHub](notebooks/rw/JW_repeat_expansion_figures.ipynb) | Repeat expansion analysis and figures |
| kvg_firth_logistic_regression.ipynb | [GitHub](notebooks/rw/kvg_firth_logistic_regression.ipynb) | Firth logistic regression analysis |
| kvg_pmi_skip_participants.ipynb | [GitHub](notebooks/rw/kvg_pmi_skip_participants.ipynb) | PMI participant filtering analysis |