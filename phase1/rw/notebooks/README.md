# Phase 1 phenotype and figure notebooks (Researcher Workbench)

These notebooks run in a JupyterLab environment in the *All of Us* Researcher Workbench Controlled Tier ([getting started](../../../GETTING_STARTED.md#researcher-workbench)). They use Controlled Tier EHR and cohort tables, plus the Terra outputs the DRC [releases to the Workbench](../../../GETTING_STARTED.md#3-data-release-to-the-researcher-workbench). Each notebook states the machine size to use, and Hail notebooks need an environment with Hail (Spark) support. Callset QC notebooks are in [terra/notebooks/](../../terra/notebooks/README.md).

The copies in this repository have outputs removed, and data paths are written relative to the `AOU_LR_DATA` environment variable, which you set before running. See [GETTING_STARTED.md](../../../GETTING_STARTED.md#4-run-the-workbench-notebooks).

| Notebook | Needs | What it does |
|----------|-------|--------------|
| [00_init_subset_vds.ipynb](00_init_subset_vds.ipynb) | Controlled Tier short-read VDS | Subsets the v7 VDS to the 1,027 long-read participants |
| [01_dataset_summary_table.ipynb](01_dataset_summary_table.ipynb) | Cohort Builder, QC metrics | Illumina, PacBio, and ONT summary table |
| [02_sampling_map.ipynb](02_sampling_map.ipynb) | ZIP3 tables | Sampling-location map |
| [03_length_distributions.ipynb](03_length_distributions.ipynb) | Terra hifiasm reports | Figure 1 assembly-length panel |
| [04_assembly_figure.ipynb](04_assembly_figure.ipynb) | Terra assemblies and an example BAM | Assembly metrics and an example alignment |
| [05_pca_figure.ipynb](05_pca_figure.ipynb) | Precomputed ancestry loadings | Figure 1 ancestry PCA |
| [06_pmi_skip_participants.ipynb](06_pmi_skip_participants.ipynb) | Cohort Builder, ancestry PCA | Participants recorded as PMI: Skip on the ancestry PCA |
| [07_variant_inventory_table.ipynb](07_variant_inventory_table.ipynb) | Short-read and long-read callsets | Variant-class counts |
| [08_short_read_svs_table.ipynb](08_short_read_svs_table.ipynb) | GATK-SV VCF | Short-read SV counts for the overlapping samples |
| [09_omop_phenotypes.ipynb](09_omop_phenotypes.ipynb) | EHR / OMOP | Condition prevalence for the 1,027 participants |
| [10_firth_logistic_regression.ipynb](10_firth_logistic_regression.ipynb) | EHR, survey, Locityper genotypes | Phenotype table for Firth regression / PheWAS |
| [11_sv_disease_associations.ipynb](11_sv_disease_associations.ipynb) | EHR and SV genotypes | Fisher tests of GWAS-LD SVs against OMOP conditions |
| [12_cyp2d6.ipynb](12_cyp2d6.ipynb) | Haplotype assemblies | CYP2D6 star-allele assignment. Several cells are `sys.argv` scripts pasted into the notebook, so this copy is not a step-by-step run. The procedure is in the Methods *CYP2D6 variation* section. |
| [13_repeat_expansion_figures.ipynb](13_repeat_expansion_figures.ipynb) | Assemblies and BAMs | FMR1 / HTT haplotig extraction. Clones `EichlerLab/AoU_WDL` at runtime. |
