# Structural variant interpretation in the All of Us Program using long-read sequencing

These analyses are **stage e** of the Phase 1 pipeline. They run after Terra stages a–d. PheWAS (section 5) queries EHR tables and must run in a Researcher Workbench terminal. The other sections run wherever their inputs are available: the imputed SV genotypes plus public MAGE and 1000 Genomes files. See [GETTING_STARTED.md](../../../../GETTING_STARTED.md) and [PIPELINE.md stage e](../../../PIPELINE.md#e-downstream-analyses--figures).

Each script starts with a comment block describing its inputs, outputs, and usage.

## Install

These scripts are the only code in this repository that needs an install. From the repository root:

```bash
conda env create -f phase1/rw/scripts/AoU_SV_analysis/environment.yml
conda activate aou-sv-analysis
export SCRIPTS=$PWD/phase1/rw/scripts/AoU_SV_analysis
```

Then `cd` into a working directory that holds your input files. Scripts without command-line flags read and write the fixed file names given below in the current directory. The commands below call the scripts through `$SCRIPTS`. Inputs are Controlled Tier or public callsets, not files in git. Some steps need an input prepared outside these scripts; each such input is described where it is used.

The environment provides Python (`pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `statsmodels`, `upsetplot`, `adjustText`, `google-cloud-bigquery`, `pandas-gbq`) plus `bedtools`, `bcftools`, `vcftools`, and `plink2`.

Two tools are not in that environment:

| Tool | Role | How to run it |
|------|------|----------------|
| [CADD-SV](https://github.com/kircherlab/CADD-SV) | PHRED-scaled SV impact scores | Run the published CADD-SV Snakemake workflow upstream. This repository only plots the resulting BED. |
| [CAVIAR](https://github.com/fhormoz/caviar) | Fine-mapping of SV and SNV signals | Build the upstream tool and put its `CAVIAR` binary on your `PATH`; both `run_caviar.sh` scripts call it by that name. The manuscript uses CAVIAR default settings. |

## 1. CADD-SV plot

1. Produce `CADD_output.bed` with the upstream CADD-SV Snakemake workflow.
2. Reshape scores so that the fourth column is `ID|...|score`:

```bash
grep -v CADD-SV_PHRED-score CADD_output.bed | awk '{OFS="\t"}{print "chr"$1, $2, $3, $5"|"$4"|"$6}' > SV_CADD.score.bed
```

3. Plot:

```bash
python $SCRIPTS/plot_cadd_sv.py \
  --cadd SV_CADD.score.bed \
  --summary sv_summary.tsv \
  --known_ids known_sv_ids.txt \
  --output Cadd_sv_sample_summary.png
```

`--summary` is tab-separated with columns `Variant_ID`, `Sample_Count`, `Score`, `Sensitivity`, `Sample_IDs`, `SVTYPE`; the plot uses only `Variant_ID` and `Sample_Count`. `--known_ids` lists SVs seen in earlier callsets, one per line (Supplementary Methods: *Variant annotation and comparison to external datasets*).

## 2. Population-scale summaries

Samples are 2,504 unrelated 1000 Genomes participants. Genotypes are the KAGE + GLIMPSE imputed SV callset, `SV_genotype.vcf.gz`. `1KG_population` lists each sample's superpopulation (columns `Sample_name`, `Superpopulation_code`; from the [1000 Genomes data portal](https://www.internationalgenome.org/data/)).

1. Genotype table. Writes `SV_genotype.stats.gz`, which steps 2, 3 and 6 read:

```bash
bash $SCRIPTS/population-scale_analysis/stats.sh
```

2. Per-SV sample counts, written to `SV_summary.txt`:

```bash
python $SCRIPTS/population-scale_analysis/sv_sample_counts.py
```

3. SV counts per participant and a violin plot across the five superpopulations:

```bash
python $SCRIPTS/population-scale_analysis/sv_count_per_sample.py
```

4. Coding-sequence overlap, then the UpSet plot. `CDS_annotation.bed` has four columns: chromosome, start, end, gene.

```bash
bedtools intersect -a SV_genotype.vcf.gz -b CDS_annotation.bed -wao \
  | awk '{OFS="\t"}{if($NF!=0) print $1, $2, $3, $(NF-3), $(NF-2), $(NF-1)}' \
  | sort | uniq > SV_CDS_ID.txt
python $SCRIPTS/population-scale_analysis/Gene_upset.py
```

5. Regulatory-element overlap, then the UpSet plot. `regulatory_element_annotation.bed` has four columns: chromosome, start, end, element ID. The script reads two columns, the SV ID and the element ID:

```bash
bedtools intersect -a SV_genotype.vcf.gz -b regulatory_element_annotation.bed -wao \
  | awk '{OFS="\t"}{if($NF!=0) print $3, $(NF-1)}' \
  | sort | uniq > SV_regulatory_element_ID.txt
python $SCRIPTS/population-scale_analysis/regulatory_element_upset.py
```

6. Allele-frequency comparison, then Fst. `AF_comparison.py` writes `SV_AF_AFR_vs_NonAFR.txt`. `AFR_samples` and `Otherpop_samples` list the AFR and non-AFR sample names, one per line. `Fst_distribution.py` also reads `SV_IDs`, a CSV with a column `SV_ID` listing the SVs in the same order as `SV_genotype.vcf.gz`, and writes `fst_histogram.png`:

```bash
python $SCRIPTS/population-scale_analysis/AF_comparison.py
vcftools --gzvcf SV_genotype.vcf.gz --weir-fst-pop AFR_samples --weir-fst-pop Otherpop_samples --out AFR_vs_nonAFR
python $SCRIPTS/population-scale_analysis/Fst_distribution.py
```

## 3. eQTL mapping

Cohort: 731 1000 Genomes samples with matched RNA-seq (26 populations). Expression `inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz` and covariates `eQTL_covariates.tab.gz` are from [MAGE](https://github.com/mccoy-lab/MAGE). SNV genotypes: [1000 Genomes high-coverage phased panel](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/). Samples are matched by name across files, and only samples present in all inputs are used.

1. Inputs for the association test.
   - SV–gene pairs within 1 Mb. `gene_annotation.bed` has the gene ID in its last column:

     ```bash
     bedtools window -w 1000000 -a SV_genotype.vcf.gz -b gene_annotation.bed \
       | awk 'BEGIN{print "SV,gene"}{print $3","$NF}' > SV_gene_1Mb
     ```

   - SV genotypes as a CSV with SV IDs in the first column and one GT string per sample:

     ```bash
     (printf 'VariantID'; bcftools query -l SV_genotype.vcf.gz | awk '{printf ",%s", $1}'; echo) > SV_genotypes.csv
     bcftools query -f '%ID[,%GT]\n' SV_genotype.vcf.gz >> SV_genotypes.csv
     ```

   - Expression as a CSV with gene IDs in the first column and one column per sample (`TMM_expression.csv`), converted from the MAGE BED.

2. Association test. The covariates file is used as distributed (tab-separated, covariate names in the first column, one column per sample):

```bash
python $SCRIPTS/eQTL_mapping/run_eQTL.py \
  --geno SV_genotypes.csv \
  --expr TMM_expression.csv \
  --covar eQTL_covariates.tab.gz \
  --pairs SV_gene_1Mb \
  --out-pairs eQTL_result.csv \
  --out-bh eQTL_result.bh.csv
```

3. Summary figure. `gene_names.tsv` maps gene IDs to names (tab-separated, no header); `omim_genes.txt` lists OMIM gene names, one per line:

```bash
python $SCRIPTS/eQTL_mapping/eQTL_summary.py \
  --eqtl eQTL_result.csv \
  --eqtl_bh eQTL_result.bh.csv \
  --gene-name gene_names.tsv \
  --gene-list omim_genes.txt \
  --out-png eQTL_summary.png
```

4. Fine-mapping. `caviar_inputs.py` takes per-variant eQTL z-scores for SVs and SNPs, which `run_eQTL.py` does not produce. Each z-score file is tab-separated with no header: variant ID, gene ID, z-score, p-value, position. The genotype files hold numeric dosages (0, 1, 2), whitespace-separated, with a header row of `VariantID` and sample names. The SNP files are gzip-compressed. `caviar_inputs.py` writes `variant.zscore`, `variant.ld`, and `variant.list` into one folder per gene under `caviar/`, which `run_caviar.sh` reads:

```bash
python $SCRIPTS/eQTL_mapping/caviar_inputs.py \
  --sv-z eQTL_SV_zscore.txt \
  --sv-gt SV_dosages.txt \
  --snp-z eQTL_SNP_zscore.txt.gz \
  --snp-gt SNP_dosages.txt.gz \
  --out-dir caviar/
bash $SCRIPTS/eQTL_mapping/run_caviar.sh
```

5. Case studies. `Genotypes_expression.py` reads fixed file names: `eQTL_result.bh.csv` from step 2, `genotypes.csv` (the SV genotype CSV from step 1), `TMM_expression.csv`, and `gene_name` (the gene ID-to-name table from step 3). `Manhattan_plot.py` is set up for *BID* (settings near the end of the script) and reads `eQTL_SV_zscore.BID.txt`, `eQTL_SNP_zscore.BID.txt.gz`, `SV_genotypes.BID.txt`, and `SNP_genotypes.BID.txt.gz` (the step 4 formats, restricted to the locus), plus `genes.bed.gz` and `BID_exons.bed`:

```bash
python $SCRIPTS/eQTL_mapping/Genotypes_expression.py --gene GENE_ID --SV SV_ID
python $SCRIPTS/eQTL_mapping/Manhattan_plot.py
```

## 4. SVs in LD with GWAS-significant SNPs

GWAS Catalog v1.0 genome-wide significant variants: [GWAS Catalog downloads](https://www.ebi.ac.uk/gwas/docs/file-downloads). Convert them to a BED file, `GWAS_snps.bed.gz`, with columns `chromosome`, `SNP_start`, `SNP_end`, `variantID|MappedGene|Trait|P-value|RiskAllele`. SV genotypes are the imputed callset on 2,504 unrelated 1000 Genomes samples used in section 2. SNP genotypes are the high-coverage panel linked in section 3.

1. Pairs within 100 kb (SV ID, GWAS info):

```bash
bedtools window -a SV_genotype.vcf.gz -b GWAS_snps.bed.gz -w 100000 \
  | awk -F"\t" '{print $3"\t"$NF}' | sort | uniq > SVs_near_GWAS_SNPs_100kb.txt
```

2. LD. The genotype files are tab-separated with genotypes coded `0`, `1`, or `2`. The first column is `VariantID` in the SV file and `Variant_ID` in the SNP file; SNP IDs must match the first field of the GWAS info. Pairs with r² ≥ 0.5 are written with r² appended:

```bash
python $SCRIPTS/LD_GWAS_catalog/LD_computation.py \
  --sv-gt SV_genotype.tsv \
  --snp-gt SNP_genotype.tsv.gz \
  --pairs SVs_near_GWAS_SNPs_100kb.txt \
  --out-pairs-ld sv_snp_ld.txt
```

3. Disease and disorder traits. `disease_terms.txt` lists disease or disorder terms, one per line:

```bash
cut -f2 sv_snp_ld.txt | cut -d'|' -f3 | sort -u > gwas_traits.txt
python $SCRIPTS/LD_GWAS_catalog/match_traits.py \
  --ref disease_terms.txt \
  --gwas gwas_traits.txt \
  --out matched_traits.txt
```

4. Summary figure. The script file name is `assocaition_summary.py`. `--all` is the output of step 2. `--genic` has the same format, restricted to pairs where the SV lies in the GWAS variant's mapped gene:

```bash
python $SCRIPTS/LD_GWAS_catalog/assocaition_summary.py \
  --all sv_snp_ld.txt \
  --genic sv_snp_ld.genic.txt \
  --disease matched_traits.txt \
  --out-fig ld_summary.png
```

5. LD heatmap. `Input_preparation.sh` is set up for the SV `chr7-75493495-allele831405-64` (edit the variables at its top for another locus). It reads `GWAS_snps.bed.gz`, `SNP_genotype.tsv.gz`, `SV_genotype.tsv` (the step 2 genotype files), and `SV_pos` (tab-separated `CHROM`, `POS`, `SV_ID`, with a header line). The heatmap script also reads `Genes.bed` (chromosome, start, end, gene, strand) and `associated_snps.txt` (SNP IDs to outline, one per line):

```bash
bash $SCRIPTS/LD_GWAS_catalog/Input_preparation.sh
python $SCRIPTS/LD_GWAS_catalog/LD_heatmap_compuation.py
```

## 5. Phenome-wide association

These steps use the Researcher Workbench. `extract_AoU_condition.py` and `sample_info.py` query Workbench BigQuery through the `WORKSPACE_CDR` environment variable. `sample_IDs.txt` lists the study participants' `person_id` values, one per line.

1. Extract conditions. `extract_AoU_condition.py` reads `condition_occurrence` for the participants in `sample_IDs.txt`, maps concept IDs through `athena_ohdsi_SNOMED/CONCEPT.csv` ([OHDSI Athena](https://athena.ohdsi.org/vocabulary/list)), keeps standard SNOMED disorder concepts, and writes `AoU_samples_EHR_disorder.csv`:

```bash
python $SCRIPTS/pheWAS/extract_AoU_condition.py
```

2. Participant demographics (sex at birth, date of birth, race, ethnicity), written to `AoU_samples_info.csv`:

```bash
python $SCRIPTS/pheWAS/sample_info.py
```

3. PLINK2 fileset of the SV genotypes, read by step 6:

```bash
plink2 --vcf input.vcf.gz --make-pgen --out AoU_SV
```

4. Phenotype file. `disease_concepts.txt` lists the condition concept IDs to test, one per line; each becomes a phenotype column (2 = case, 1 = control):

```bash
python $SCRIPTS/pheWAS/phenotype.py \
  --disease-file disease_concepts.txt \
  --sample-file sample_IDs.txt \
  --ehr-file AoU_samples_EHR_disorder.csv \
  --out-pheno pheno.txt \
  --fid-value 0 \
  --case-value 2 \
  --control-value 1
```

5. Covariates. `ancestry_preds.tsv` is the All of Us ancestry prediction table (tab-separated, with columns `research_id` and `pca_features`):

```bash
python $SCRIPTS/pheWAS/covar.py \
  --sample-file sample_IDs.txt \
  --pc-file ancestry_preds.tsv \
  --info-file AoU_samples_info.csv \
  --out-covar covar.txt \
  --reference-date 2020-01-01 \
  --n-pcs 10
```

6. Association test. `association_testing.sh` runs PLINK2 Firth regression of every SV against each condition listed in `disease.txt`, with covariates `age`, `sex`, and `PC1`–`PC10`, and writes `<condition>.<condition>.glm.firth.gz`. For a sex-specific run, drop `sex` from `--covar-name` inside that script.

```bash
cp disease_concepts.txt disease.txt
bash $SCRIPTS/pheWAS/association_testing.sh
```

7. Fine-mapping of SV–condition associations. This step needs three inputs that the scripts above do not produce:
   - `variant_disease_pairs.tsv.gz`: tab-separated; SV ID, condition, then the IDs of the SNPs to fine-map with it.
   - `AoU_SV_SNPs.pgen/.pvar/.psam`: a PLINK2 fileset holding both the SVs and those SNPs.
   - SNP association results for each pair, in `snp_disease_plink/<SV>_<condition>.<condition>.glm.firth.gz`. These come from the same Firth model as step 6, run on the pair's SNPs, for example:

     ```bash
     plink2 --pfile AoU_SV_SNPs --extract SNPS_FOR_PAIR.txt --pheno pheno.txt --pheno-name CONDITION --covar covar.txt --covar-name age sex PC1-PC10 --covar-variance-standardize --glm firth hide-covar --out snp_disease_plink/SV_ID_CONDITION
     ```

   `compute_ld_r.sh` writes LD for each pair to `LD_r/`. `caviar_input.py` reads `SV_pheWAS_summary.txt` (tab-separated, with columns `SV`, `Disease`, and `Z_disc`, the SV's z-score), the SNP results, and the LD files. It writes one folder per pair under `caviar_input/`, which `run_caviar.sh` reads:

```bash
bash $SCRIPTS/pheWAS/compute_ld_r.sh
python $SCRIPTS/pheWAS/caviar_input.py \
  --sv-pheWAS-file SV_pheWAS_summary.txt \
  --snp-dir snp_disease_plink \
  --ld-dir LD_r \
  --top-n 1000 \
  --out-dir caviar_input
bash $SCRIPTS/pheWAS/run_caviar.sh
```
