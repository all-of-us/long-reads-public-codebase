# Phase 1 callset QC notebooks (Terra)

These notebooks run in the Terra workspace that holds the workflow outputs. They read Terra products: joint-call VCFs, Hail MatrixTables, and assembly tables. Hail notebooks need a Hail environment when you create the cloud environment. Phenotype and figure notebooks are in [rw/notebooks/](../../rw/notebooks/README.md).

The copies in this repository have outputs removed, and data paths are written relative to the `AOU_LR_DATA` environment variable, which you set before running. See [GETTING_STARTED.md](../../../GETTING_STARTED.md#4-run-the-workbench-notebooks).

| Notebook | Needs | What it does |
|----------|-------|--------------|
| [00_init_subset_vds.ipynb](00_init_subset_vds.ipynb) | Terra joint-call products, short-read VDS | Subsets the short-read v7 VDS to the Phase 1 long-read cohort and runs Hail `sample_qc` |
| [01_describe_hail_matrix_tables.ipynb](01_describe_hail_matrix_tables.ipynb) | Terra GLnexus MatrixTables | Sample and variant counts, QC, LD pruning, PCA, relatedness |
| [02_read_length_distributions.ipynb](02_read_length_distributions.ipynb) | Terra read-length histograms | PacBio length and coverage curves |
| [03_examine_assemblies.ipynb](03_examine_assemblies.ipynb) | Terra assembly table | Contig count, N50, duplication ratio, and length versus coverage |
| [04_examine_small_variants.ipynb](04_examine_small_variants.ipynb) | Per-caller small-variant QC tables | Ti/Tv and SNP counts for DeepVariant, PAV, and Clair3 |
| [05_examine_structural_variants.ipynb](05_examine_structural_variants.ipynb) | Per-caller SV stats | PAV, PBSV, and Sniffles counts by type and length |
| [06_callset_qc_python.ipynb](06_callset_qc_python.ipynb) | Joint SV VCFs | Writes SVTYPE and SVLEN tables for the R notebook |
| [07_callset_qc_r.ipynb](07_callset_qc_r.ipynb) | Tables from `06_callset_qc_python.ipynb` | SVs per genome, cohort counts, length densities |
| [08_sv_callset_inventory.ipynb](08_sv_callset_inventory.ipynb) | Callset locations | Starting inventory of long-read SV callsets |
| [09_hiphase_results.ipynb](09_hiphase_results.ipynb) | HiPhase VCFs from Terra | Spot-check of phased versus unphased genotypes |
| [10_compute_sfs_grch38.ipynb](10_compute_sfs_grch38.ipynb) | Terra GLnexus MatrixTable | Site-frequency spectra for autosomal biallelic SNVs and indels |
| [11_pca_hgdp_tgp.ipynb](11_pca_hgdp_tgp.ipynb) | HGDP+1kGP MatrixTable | Global reference PCA used as a projection panel |
| [12_pca.ipynb](12_pca.ipynb) | Terra joint-call MatrixTable | PCA of Phase 1 samples projected onto gnomAD and the short-read AoU PCA |
| [13_pca_figure.ipynb](13_pca_figure.ipynb) | Precomputed ancestry loadings | Figure 1 ancestry PCA |
| [14_sv_summary_table.ipynb](14_sv_summary_table.ipynb) | Joint SV VCFs | Short-read GATK-SV versus long-read joint SV counts |
| [15_variant_inventory_table.ipynb](15_variant_inventory_table.ipynb) | Joint SV VCFs and Hail tables | SNV, indel, and SVTYPE counts |

`00_init_subset_vds`, `13_pca_figure`, and `15_variant_inventory_table` have Researcher Workbench counterparts with the same analysis name (`00_init_subset_vds`, `05_pca_figure`, and `07_variant_inventory_table` in `rw/notebooks/`). Those versions differ and read Controlled Tier inputs.
