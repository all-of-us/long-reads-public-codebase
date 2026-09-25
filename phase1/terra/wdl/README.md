# Phase 1 WDL workflows

Copies of the Phase 1 WDL workflows, taken from their upstream repositories. Launch order is in [GETTING_STARTED.md](../../../GETTING_STARTED.md#launch-order), and each stage is described in [PIPELINE.md](../../PIPELINE.md).

## Layout

```
<upstream>/<Workflow>/
├── <Workflow>.wdl       # entry point: run or validate this file
├── repo/                # WDLs it imports, at their upstream paths (absent if it imports nothing)
├── SOURCE.txt           # upstream repo, branch, upstream path, Dockstore entry, local patches
└── example_inputs.json  # example inputs from a Phase 1 run (GRCh38)
```

The entry point is the upstream workflow file with two changes: its relative `import` lines are rewritten to point into `repo/`, so it runs from this directory, and a header comment is added at the top. The header says what the workflow does, its pipeline stage and launch step, whether it runs per sample, per chromosome, or per cohort, its main inputs and outputs, and the workflows before and after it. `SOURCE.txt` records the upstream path in `entry:`.

Every entry point passes `womtool validate`. To check them all, run [scripts/validate_wdls.sh](../../../scripts/validate_wdls.sh), which is also run in CI.

## Example inputs

Each `example_inputs.json` holds the inputs of one Phase 1 GRCh38 run of that workflow, in the `Workflow.input` form that `womtool inputs` produces. Fixed values (thresholds, tool arguments, Docker images, reference files) are given as used. Values such as `this.aligned_bam` are Terra references to per-sample data (a column of the table the workflow ran on) and `workspace.<name>` to a workspace-level resource; replace them with your own files. Private bucket names are replaced by labels such as `<workspace-bucket:AoU_DRC_WGS_LongReads_PacBio>`, Terra submission and workflow IDs in paths by `<submission-id>` and `<workflow-id>`, and participant identifiers are removed. The inter-sample SV merge examples follow the callset filtered at calibration sensitivity 0.7.

The examples were exported from the Phase 1 Terra workspaces with `scripts/export_terra_configs.py` (see its help).

## How closely each copy matches what ran

Terra no longer keeps the WDL submitted by 2022–2025 runs, so each copy was checked against the upstream commit on the branch that ran, as of the run date, and against the input names used in its Phase 1 runs. This is a best-effort reconstruction.

| Match | Workflows |
|-------|-----------|
| Same code, or differences only in comments, timing wrappers, or new optional inputs whose defaults match the run | All callset-integration workflows except those below; PBCCSWholeGenome; LRJointCallGVCFs; KAGEPanelWithPreprocessing; KAGECasePerChromosomeFlexscattered; ConcatAndEvaluateReduced |
| Resources only | pav: the main Phase 1 runs used branch `sh_more_resources`. This copy is `sh_more_resources_pete`, which was used for reruns and only raises threads to 64 and memory to 256 GB. |
| Replaced with the version that ran | PAV2SVs (see local patches below) |
| Code differs; which version produced the published result is unconfirmed | ChromosomePhasedPanelCreationFromHiPhase: the GRCh38 panel configuration ran branch `sl_chr6_v1`, whose `bcftools norm` step skips normalization; this copy (`sl_v1.1b`) normalizes against the reference. LRMergeSVVCFs (helper): ran `fc_sveval`, 17 lines differ from this copy (`kvg_sveval`). |
| Inputs match the Phase 1 runs, but the exact commit cannot be verified (2022 runs on branches since rewritten or deleted) | PBAssembleWithHifiasm, SummarizeSmallVariants, HiFiBamToFastQ, MergeFastqs, MapAssemblyContigs, LongReadsContaminationEstimation, BuildTempLocalFpStore |
| Later version than what ran (QC only; the Phase 1 runs used inputs this copy lacks) | SexCheck, VerifyFingerprintCCSSample |
| No example inputs (no Phase 1 run found) | PhysicalPhasing, V11HprcAnalysis, SmallVariantsBasicMetrics |

## Local patches

Four workflows differ from their upstream branch. Each change is recorded in a `local_patch:` line in that workflow's `SOURCE.txt`.

| Workflow | Upstream problem | Fix |
|----------|------------------|-----|
| PAV2SVs | Upstream `main` changed after the Phase 1 run: later commits split multiallelic records, change the Docker image, and drop the SNP output. | Used `PAV2SVs.wdl` from `9a31d42`, the last commit on or before the last Phase 1 run. |
| MergeFastqs, MapAssemblyContigs | On branch `sh_hificnv`, these import `tasks/Utils.wdl`, `tasks/Finalize.wdl`, and `tasks/AlignReads.wdl`. Those paths were removed in a 2023 reorganisation, and the task `MergeFastqGzs` was removed with them. | Added those task files from the last commit with that layout (`ad47b56`). |
| ConcatAndEvaluateReduced | On branch `sl_v1.1b`, `GLIMPSEBatchedCasePerChromosome.wdl` changed its inputs (commit `1159603`) without updating its caller, `LeaveOutEvaluationReduced.wdl`. | Used `GLIMPSEBatchedCasePerChromosome.wdl` from `13edc1a`, the version its caller expects. |

## Index

Descriptions are grouped by upstream repository. Stage letters refer to [PIPELINE.md](../../PIPELINE.md).

### callset-integration ([fabio-cunial/callset_integration](https://github.com/fabio-cunial/callset_integration), branch `main`)

Launch order and runtime requirements: [callset-integration/README.md](callset-integration/README.md).

| Workflow | Stage | Description |
|----------|-------|-------------|
| [TruvariIntrasampleV11](callset-integration/TruvariIntrasampleV11/TruvariIntrasampleV11.wdl) | b | Intra-sample merge of PBSV, Sniffles2, and PAV calls (upstream file `TruvariIntrasampleFabio.wdl`, workflow `TruvariIntrasampleFabio`) |
| [GetRegenotypedVcfKanpigV11](callset-integration/GetRegenotypedVcfKanpigV11/GetRegenotypedVcfKanpigV11.wdl) | b | Intra-sample Kanpig re-genotyping and annotation |
| [FilterIntrasampleDevV11](callset-integration/FilterIntrasampleDevV11/FilterIntrasampleDevV11.wdl) | b | Isolation-forest scoring against HPRC truth |
| [FilterAndSplit](callset-integration/FilterAndSplit/FilterAndSplit.wdl) | b | Calibration-sensitivity filter, then split by chromosome |
| [TruvariIntersample](callset-integration/TruvariIntersample/TruvariIntersample.wdl) | c | Inter-sample merge, per chromosome |
| [ConcatenateChromosomes](callset-integration/ConcatenateChromosomes/ConcatenateChromosomes.wdl) | c | Concatenate per-chromosome merged VCFs |
| [RemoveSamplesV11](callset-integration/RemoveSamplesV11/RemoveSamplesV11.wdl) | c | Prepare the merged VCF for re-genotyping |
| [GetRegenotypedVcfKanpigMergedV11](callset-integration/GetRegenotypedVcfKanpigMergedV11/GetRegenotypedVcfKanpigMergedV11.wdl) | c | Re-genotype the merged VCF per sample with Kanpig |
| [MergeRegenotypedIntersampleVcf](callset-integration/MergeRegenotypedIntersampleVcf/MergeRegenotypedIntersampleVcf.wdl) | c | Merge per-sample re-genotyped VCFs into the cohort VCF |
| [PAV2SVs](callset-integration/PAV2SVs/PAV2SVs.wdl) | a | Convert PAV output to standard SV calls |
| [TruvariCollapse](callset-integration/TruvariCollapse/TruvariCollapse.wdl) | — | Collapse duplicate or equivalent SVs into consensus calls |
| [MergeSVsSNPs](callset-integration/MergeSVsSNPs/MergeSVsSNPs.wdl) | — | Combine SVs with SNVs/indels in one VCF |
| [V11HprcAnalysis](callset-integration/V11HprcAnalysis/V11HprcAnalysis.wdl) | — | HPRC-only benchmarking for the manuscript |

### long-read-pipelines ([broadinstitute/long-read-pipelines](https://github.com/broadinstitute/long-read-pipelines))

| Workflow | Branch | Stage | Description |
|----------|--------|-------|-------------|
| [HiFiBamToFastQ](long-read-pipelines/HiFiBamToFastQ/HiFiBamToFastQ.wdl) | `sh_wg_extract_fastq` | a | Convert PacBio HiFi BAM to FASTQ |
| [MergeFastqs](long-read-pipelines/MergeFastqs/MergeFastqs.wdl) | `sh_hificnv` | a | Merge FASTQ files into one set per sample (patched; see above) |
| [PBAssembleWithHifiasm](long-read-pipelines/PBAssembleWithHifiasm/PBAssembleWithHifiasm.wdl) | `sh_update_hifiasm` | a | Assemble HiFi reads with hifiasm |
| [MapAssemblyContigs](long-read-pipelines/MapAssemblyContigs/MapAssemblyContigs.wdl) | `sh_hificnv` | a | Align assembly contigs to the reference (patched; see above) |
| [PBCCSWholeGenome](long-read-pipelines/PBCCSWholeGenome/PBCCSWholeGenome.wdl) | `main` | a | Small variants (DeepVariant) and SVs (PBSV, Sniffles2) from aligned HiFi BAMs |
| [SummarizeSmallVariants](long-read-pipelines/SummarizeSmallVariants/SummarizeSmallVariants.wdl) | `kvg_fix_t2t_joint_calls` | a | Summarize small-variant calls. Earlier docs call it SummarizeDVPSmallVariants or SummarizePAVSmallVariants. |
| [SmallVariantsBasicMetrics](long-read-pipelines/SmallVariantsBasicMetrics/SmallVariantsBasicMetrics.wdl) | `sh_variants_ss_metrics` | a (QC) | Small-variant metrics (run as 1074.T2T.SmallVariantsBasicMetrics) |
| [CollectSingleSampleSVvcfMetrics](long-read-pipelines/CollectSingleSampleSVvcfMetrics/CollectSingleSampleSVvcfMetrics.wdl) | `sh_gatk_sv_vcfmetrics` | a (QC) | Per-sample SV metrics |
| [LongReadsContaminationEstimation](long-read-pipelines/LongReadsContaminationEstimation/LongReadsContaminationEstimation.wdl) | `sh_more_atomic_qc` | a (QC) | Contamination estimate |
| [BuildTempLocalFpStore](long-read-pipelines/BuildTempLocalFpStore/BuildTempLocalFpStore.wdl) | `sh_exhaust_au` | a (QC) | Temporary fingerprint store for identity checks |
| [VerifyFingerprintCCSSample](long-read-pipelines/VerifyFingerprintCCSSample/VerifyFingerprintCCSSample.wdl) | `sh_sample_fp` | a (QC) | Verify sample identity by fingerprinting (upstream `VerifyFingerprint`) |
| [SexCheck](long-read-pipelines/SexCheck/SexCheck.wdl) | `sh_more_atomic_qc` | a (QC) | Reported versus genetic sex (upstream `SexCheckNaive`) |
| [LRJointCallGVCFs](long-read-pipelines/LRJointCallGVCFs/LRJointCallGVCFs.wdl) | `main` | c | Joint-genotype gVCFs into a cohort VCF with GLnexus |
| [LRMergeSVVCFs](long-read-pipelines/LRMergeSVVCFs/LRMergeSVVCFs.wdl) | `kvg_sveval` | — | Merge several SV callsets into one |
| [MergePhasedVCF](long-read-pipelines/MergePhasedVCF/MergePhasedVCF.wdl) | `hangsu_phasing` | — | Merge phased VCFs |

### pav-wdl ([broadinstitute/pav-wdl](https://github.com/broadinstitute/pav-wdl))

| Workflow | Branch | Stage | Description |
|----------|--------|-------|-------------|
| [pav](pav-wdl/pav/pav.wdl) | `sh_more_resources_pete` | a | Assembly-based SV and indel calling with PAV |

### lrma-aou1-panel-creation ([broadinstitute/lrma-aou1-panel-creation](https://github.com/broadinstitute/lrma-aou1-panel-creation))

| Workflow | Branch | Stage | Description |
|----------|--------|-------|-------------|
| [PhysicalPhasing](lrma-aou1-panel-creation/PhysicalPhasing/PhysicalPhasing.wdl) | `hs_physicalphasing` | c | Physically phase SNVs/indels and SVs in one sample with HiPhase (the WDL names its workflow `PhasicalPhasing`) |
| [ChromosomePhasedPanelCreationFromHiPhase](lrma-aou1-panel-creation/ChromosomePhasedPanelCreationFromHiPhase/ChromosomePhasedPanelCreationFromHiPhase.wdl) | `sl_v1.1b` | c | Per chromosome: statistical phasing with SHAPEIT4, removal of colliding variants, and creation of the bubble-graph panel |
| [ConcatAndEvaluateReduced](lrma-aou1-panel-creation/ConcatAndEvaluateReduced/ConcatAndEvaluateReduced.wdl) | `sl_v1.1b` | c | Concatenate per-chromosome panels; leave-out and vcfdist evaluations (patched; see above). Earlier docs call it ConcatAndEvaluate. |
| [KAGEPanelWithPreprocessing](lrma-aou1-panel-creation/KAGEPanelWithPreprocessing/KAGEPanelWithPreprocessing.wdl) | `sl_v1.1b` | d | Per chromosome, build the KAGE k-mer index and count model from the panel |
| [KAGECasePerChromosomeFlexscattered](lrma-aou1-panel-creation/KAGECasePerChromosomeFlexscattered/KAGECasePerChromosomeFlexscattered.wdl) | `sl_v1.1b` | d | Genotype one sample against the panel with KAGE |
| [GLIMPSEBatchedCasePerChromosomeSingleBatch](lrma-aou1-panel-creation/GLIMPSEBatchedCasePerChromosomeSingleBatch/GLIMPSEBatchedCasePerChromosomeSingleBatch.wdl) | `sl_v1.1b` | d | Phase and impute a batch of genotyped samples with GLIMPSE |
| [HierarchicallyMergeVcfs](lrma-aou1-panel-creation/HierarchicallyMergeVcfs/HierarchicallyMergeVcfs.wdl) | `sl_v1.1b` | d | Hierarchically merge cohort VCFs with bcftools or ivcfmerge |

### sv-merging ([fabio-cunial/sv-merging](https://github.com/fabio-cunial/sv-merging), branch `main`)

| Workflow | Stage | Description |
|----------|-------|-------------|
| [TruvariIntrasample](sv-merging/TruvariIntrasample/TruvariIntrasample.wdl) | — | Compare SVs within a single sample |
| [JointCalling](sv-merging/JointCalling/JointCalling.wdl) | — | Jointly genotype SVs across a cohort |
| [MergeVCFs](sv-merging/MergeVCFs/MergeVCFs.wdl) | — | Merge VCFs |
| [OverlapGraph](sv-merging/OverlapGraph/OverlapGraph.wdl) | — | Build an overlap graph across callsets |
| [OverlapStats](sv-merging/OverlapStats/OverlapStats.wdl) | — | Overlap statistics between callsets |

### gatk-sv and sv-merge

| Workflow | Upstream | Branch | Stage | Description |
|----------|----------|--------|-------|-------------|
| [MainVcfQc](gatk-sv/MainVcfQc/MainVcfQc.wdl) | [broadinstitute/gatk-sv](https://github.com/broadinstitute/gatk-sv) | `v0.26.2-beta` | c (QC) | Quality control of cohort VCFs |
| [GraphEvaluation](sv-merge/GraphEvaluation/GraphEvaluation.wdl) | [rlorigro/sv_merge](https://github.com/rlorigro/sv_merge) | `dev` | — | Build and evaluate SV overlap graphs |

A stage of "—" marks a helper or evaluation workflow that is not on the main launch path.

## Not vendored

The workflows below were used in Phase 1 but could not be copied here. Their Dockstore entries are unpublished and their Git branches have been deleted.

| Workflow | Dockstore name | Removed branch |
|----------|----------------|----------------|
| EvaluateAssemblyHap | EvaluateAssembly | `kvg_summarize_haplotype_fastas` |
| SummarizePAVSVs, SummarizeSnifflesSVs | SummarizeStructuralVariants | `sh_baseoff_kvg_summarize_vcf` |

Docker images, the FireCloud method that `FilterIntrasampleDevV11` imports, and inputs supplied only at launch time are not vendored.
