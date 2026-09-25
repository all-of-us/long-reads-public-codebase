# All of Us Long Read Project

Public code and documentation for long-read sequencing in the [All of Us Research Program](https://allofus.nih.gov/) (AoU-LR).

This repository holds the WDL workflows, Jupyter notebooks, and analysis scripts behind population-scale structural variant (SV) discovery, reference-panel construction, short-read genotyping and imputation, and downstream interpretation.

**To install and run the code, start with [GETTING_STARTED.md](GETTING_STARTED.md).** It covers access, installation, and each stage in order.

## Where the code runs

AoU-LR data are processed in two secure cloud environments. Individual-level data never leave them, and none are included in this repository.

| Environment | Who runs it | What runs there | Code in this repository |
|-------------|-------------|-----------------|-------------------------|
| [Terra](https://app.terra.bio) | AoU [Data and Research Center (DRC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11157478/) | Primary processing: alignment, assembly, variant calling, cohort merging, reference panel, imputation, callset QC | WDL workflows, QC notebooks |
| [Researcher Workbench](https://workbench.researchallofus.org) (Controlled Tier) | AoU-LR Working Group | Phenotype analyses, figures, and tables on data copied from Terra | Notebooks, console scripts |

## Repository layout

Each phase directory is split by the environment its code runs in:

```
phase1/
├── README.md          # Phase 1 overview, pipeline figure, upstream repositories
├── PIPELINE.md        # Stages a–e: inputs, outputs, and workflows in order
├── SOFTWARE.md        # Tool versions cited in the Methods, and where each runs
├── figures/           # pipeline overview figure
├── terra/
│   ├── wdl/           # one directory per workflow (see terra/wdl/README.md)
│   └── notebooks/     # callset QC notebooks
└── rw/
    ├── notebooks/     # phenotype, figure, and table notebooks
    └── scripts/       # console scripts (the one conda install)
phase2/                # not yet released
scripts/               # maintenance tools (WDL checks, config export)
```

WDL validity is checked on every push by [.github/workflows/validate-wdls.yml](.github/workflows/validate-wdls.yml). To run the same check locally, use `scripts/validate_wdls.sh`.

## Phases

- **Phase 1** ([phase1/](phase1/README.md)): PacBio HiFi sequencing (~8×) of 1,027 participants who self-identified as Black or African American.
- **Phase 2** ([phase2/](phase2/README.md)): PacBio HiFi sequencing of 12,261 participants, plus Oxford Nanopore sequencing, across six computed continental ancestry groups. Code is not yet released.

## Citing

Citation metadata for this code, including authors and ORCIDs, are in [CITATION.cff](CITATION.cff). GitHub shows a "Cite this repository" button generated from it.

## License

See [LICENSE](LICENSE).
