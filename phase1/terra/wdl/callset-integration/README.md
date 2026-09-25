# Inter-sample SV VCF v1.1 workflows

Copied from [fabio-cunial/callset_integration](https://github.com/fabio-cunial/callset_integration) (branch `main`), also on [Dockstore](https://dockstore.org/organizations/callset_integration). Steps 1–4 are [stage b](../../../PIPELINE.md#b-per-sample-filtering) and steps 5–9 are [stage c](../../../PIPELINE.md#c-reference-panel-creation).

## Launch order (AoU + HPRC)

| Step | Workflow | Purpose | Key inputs | Key outputs |
|------|----------|---------|------------|-------------|
| 1 | [TruvariIntrasampleV11](TruvariIntrasampleV11/TruvariIntrasampleV11.wdl) | Intra-sample merge (PBSV + Sniffles2 + PAV) | Per-caller VCFs, reference | `truvari_collapsed` |
| 2 | [GetRegenotypedVcfKanpigV11](GetRegenotypedVcfKanpigV11/GetRegenotypedVcfKanpigV11.wdl) | Intra-sample Kanpig re-genotyping and annotation | Collapsed VCF, BAM | `regenotyped_kanpig` |
| 3 | [FilterIntrasampleDevV11](FilterIntrasampleDevV11/FilterIntrasampleDevV11.wdl) | Isolation-forest scoring and filtering | Collapsed and re-genotyped VCFs, HPRC training resource, [resolve.py](resolve.py) | `v11_scored_kanpig_vcf_gz` |
| 4 | [FilterAndSplit](FilterAndSplit/FilterAndSplit.wdl) | CAL_SENS filter, then split by chromosome | Scored VCF | Per-chromosome files under `destination_dir` |
| 5 | [TruvariIntersample](TruvariIntersample/TruvariIntersample.wdl) | Inter-sample merge, in parallel by chromosome | Split VCF directory, `density_counter_py` | Per-chromosome merged VCFs |
| 6 | [ConcatenateChromosomes](ConcatenateChromosomes/ConcatenateChromosomes.wdl) | Concatenate chromosomes | Per-chromosome merged directory, samples file | `concat.vcf.gz` |
| 7 | [RemoveSamplesV11](RemoveSamplesV11/RemoveSamplesV11.wdl) | Prepare for inter-sample re-genotyping | Concatenated VCF | `concat_for_regenotyping.vcf.gz` |
| 8 | [GetRegenotypedVcfKanpigMergedV11](GetRegenotypedVcfKanpigMergedV11/GetRegenotypedVcfKanpigMergedV11.wdl) | Re-genotype the merged VCF per sample | Single-column VCF, BAM | Per-sample re-genotyped VCF |
| 9 | [MergeRegenotypedIntersampleVcf](MergeRegenotypedIntersampleVcf/MergeRegenotypedIntersampleVcf.wdl) | Merge re-genotyped single-sample VCFs | List of re-genotyped VCF paths | `concat_regenotyped.vcf.gz` |

HPRC-only benchmarking and the manuscript analysis start from the intra-sample `truvari_collapsed` VCFs. They run FilterAndSplit (unfiltered split), TruvariIntersample, and ConcatenateChromosomes, then [V11HprcAnalysis](V11HprcAnalysis/V11HprcAnalysis.wdl).

## Runtime requirements not in this repository

- Pre-built Docker images (`truvari_intrasample`, `kanpig_experiments`, `hgsvc2`). Their build definitions are in the upstream repository.
- `density_counter_py` for TruvariIntersample (a chunking helper script)
- HPRC dipcall training VCFs and BEDs for FilterIntrasampleDevV11
- The GATK isolation-forest filtering method that FilterIntrasampleDevV11 imports from FireCloud at launch time (method `slee:FilterIntrasampleJointVcfFilteringDev`; the URL is in the WDL)

[resolve.py](resolve.py) is the symbolic-allele resolver. It is baked into the Docker images as `/truvari_intrasample/resolve.py` and passed to FilterIntrasampleDevV11 as `resolve_python_script`.
