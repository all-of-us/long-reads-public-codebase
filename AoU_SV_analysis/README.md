# Structural variant interpretation in the All of Us Program using long-read sequencing
## Project Overview
This project investigates the functional and clinical relevance of structural variants (SVs) identified from long-read sequencing data in the All of Us (AoU) Research Program. To evaluate their broader impact, these variants were imputed into both the 1000 Genomes Project and short-read AoU cohorts. Integration with RNA-seq and electronic health record (EHR) data enabled downstream analyses, such as eQTL mapping, linkage disequilibrium (LD) with GWAS variants, and phenome-wide association studies (PheWAS), to comprehensively characterize the regulatory and disease-related effects of variants.

### Score the effect of SVs ###
#### CADD-SV annotation
```bash
conda activate run.caddsv
snakemake  --use-conda --configfile config.yml -j 4 -n
grep -v CADD-SV_PHRED-score CADD_output.bed |awk '{OFS="\t"}{print "chr"$1, $2, $3, $5"|"$4"|"$6}' > SV_CADD.score.bed
```
#### Distribution of CADD-SV score for shared SVs in the AoU strict cohort
`plot_cadd_sv.py` 
- `--cadd`: CADD results (`SV_CADD.score.bed`)
- `--summary`: Tab-separated file with columns (`Variant_ID`, `Sample_Count`, `Score`, `Sensitivity`, `Sample_IDs`, `SVTYPE`) <br>
- `--known_ids`: IDs of SVs detected in previous callsets (Comparison strategy can refer to *Section: Variant annotation and comparison to external datasets* in the Supplementary file.) `

### Population-scale SV analysis
- **Samples:**
2,540 unrelated 1000 Genomes Project (1KGP) samples
 
- **Genotypes:**
SVs genotyped and imputed into 1KG samples using KAGE and GLIMPSE to obtain high-quality callset.

#### SV summary
`stats.sh`: SV information from the VCF file.
`sv_sample_counts.py.py`: Number of samples per SV and population summary.
'sv_count_per_sample.py':  SV counts per participant and violin plot of SV distributions across five continental groups.

#### Gene intersection
```bash
bedtools intersect -a SV_genotype.vcf.gz  -b CDS_annotation.bed  -wao |awk '{OFS="\t"}{if($NF!=0) print $1, $2, $3, $(NF-3), $(NF-2), $(NF-1)}'  |sort |uniq > SV_CDS_ID.txt
```
`Gene_upset.py`: UpSet plot of protein-coding genes with coding regions overlapping SVs across five continental groups.

#### regulatory element intersection
```bash
bedtools intersect -a SV_genotype.vcf.gz  -b regulatory_element_annotation.bed  -wao |awk '{OFS="\t"}{if($NF!=0) print $1, $2, $3, $(NF-3), $(NF-2), $(NF-1)}'  |sort |uniq > SV_regulatory_element_ID.txt
```
`regulatory_element_upset.py`: UpSet plot of regulatory elements overlapping SVs across five continental groups.

#### African vs. Non-African comparisons
`AF comparions.py`: Comparison of SV allele frequencies between African and non-African samples. <br>
F<sub>ST</sub> computation: Genetic differentiation between African and non-African groups.
```bash
vcftools --gzvcf SV_genotype.vcf.g --weir-fst-pop AFR_samples --weir-fst-pop Otherpop_samples --out AFR_vs_nonAFR
```
`Fst_distribution.py`: genome-wide F<sub>ST</sub> distribution.


### eQTL mapping ###
- **Cohort:**  
Matched DNA-seq and RNA-seq data from 731 1KG individuals, representing 26 globally distributed populations across five continents.

- **Genotypes:**  
AoU SVs genotyped and imputed using KAGE and Glimpse.
SNPs genotype data are available from: [https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)

- **Expression:**  
The expression file `inverse_normal_TMM.filtered.TSS.MAGE.v1.0.bed.gz` contains TMM-normalized expression values transformed to a normal distribution.

- **Covariates:**  
The corresponding covariate file is `eQTL_covariates.tab.gz`.

- **Data source:**  
Both expression and covariate files are available from the MAGE project([https://github.com/mccoy-lab/MAGE](https://github.com/mccoy-lab/MAGE))  

#### eQTL 
##### Generate SV–gene pairs within 1 Mb
```bash
bedtools window -w 1000000 -a SV_genotype.vcf.gz -b gene_annotation.bed |awk 'BEGIN{print "SV,gene"}{print $3","$NF}' > SV_gene_1Mb
```
##### Run eQTL analysis
`run_eqtl.py` 
- `--geno`: SV genotypes for each sample (e.g., `1|1`, `0|1`), with columns: `VariantID`, `Sample1`, `Sample2`, ...
- `--expr`: Gene expression matrix, with columns: `gene`, `Sample1`, `Sample2`, ...
- `--covar`: Sample covariate file, with columns: `id`, `Sample1`, `Sample2`, ...
- `--pairs` List of SV–gene pairs (e.g., variant_gene_1Mb)
- `---out-pairs`: Raw association results for each SV–gene pair
- `--out-bh`: Benjamini–Hochberg FDR–corrected results
##### Visualize eQTL results
`eQTL_summary.py` 
- `--eqtl`: Raw eQTL results (`eQTL_result.csv`)
- `--eqtl_bh`: FDR-corrected results (`eQTL_result.bh.csv`)
- `--gene-name`: Tab-delimited file mapping gene IDs to gene names.
- `--gene-list`: Medically relevant gene names.   
- `--out-png`: Figure summarizing eQTL findings
  
#### Fine-mapping 
##### Generate CAVIAR input files
`caviar_inputs.py`
- `--sv-z`: SV-eQTL summary statistics with columns: `VariantID`, `GeneID`, `Zscore`, `pscore`, `pos`.
- `--sv-gt`: SV genotype file with columns: `VariantID`, `Sample1`, `Sample2`, ...
- `--snp-z`: SNP-eQTL summary statistics with columns: `VariantID`, `GeneID`, `Zscore`, `pscore`, `pos`.
- `--snp-gt`: SNP genotype file with columns: `VariantID`, `Sample1`, `Sample2`, ...
- `variant.zscore`: Variant IDs and corresponding zscores.  
- `variant.ld`: LD r² matrix.  
- `variant.list`: Ordered list of variant IDs used in the analysis.
            
##### Casual variant identification 
`run_caviar.sh`: Process all variant–gene pairs using the generated input files.

##### Case study
`Genotypes_expression.py`: Relationship between genotypes and gene expression for a specific pair
- `--gene`: Gene ID in the SV–gene pair.
- `--SV`: SV ID in the SV–gene pair.

`Manhattan_plot.py`: Generates Manhattan plots for *BID*. The genotype and eQTL significance files follow the same format as `caviar_inputs.py`, but include data for the selected gene only to improve computational efficiency.


### SVs in LD with GWAS-significant SNPs  ###
- **SNVs & trait from GWAS catalog**
Extract genome-wide significant variants and their associated traits from GWAS Catalog v1.0 ([https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/docs/file-downloads)](https://www.ebi.ac.uk/gwas/docs/file-downloads) and convert them into a BED file with the following columns: `chromosome`,`SNP_start`,`SNP_end`,`variantID|MappedGene|Trait|P-value|RiskAllele`

- **SV genotypes:**  
  AoU SVs imputed using 2,504 unrelated 1KG samples.

- **SNP genotypes:**  
  Extracted from the same 1000 Genomes samples:  
  [https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)

#### Identification of SVs near GWAS-significant SNPs
```bash
bedtools window -a SV_1KG_imputed.vcf.gz -b GWAS_snps.bed.gz -w 100000| awk -F"\t" '{print $3"\t"$NF}' | sort | uniq > SVs_near_GWAS_SNPs_100kb.txt
```

##### LD calculation
`LD_computation.py`: Computes linkage disequilibrium (r²) between SV–SNP pairs.
- `--sv-gt`: SV genotypes for each sample (`0` = homozygous reference, `1` = heterozygous, `2` = homozygous alternate), with columns: `VariantID`, `Sample1`, `Sample2`, ...
- `--sbp-gt`: SNP genotypes for each sample (`0` = homozygous reference, `1` = heterozygous, `2` = homozygous alternate), with columns: `VariantID`, `Sample1`, `Sample2`, ...
- `--pairs`: SV-SNP pairs(e.g., `SVs_near_GWAS_SNPs_100kb.txt`).
- `--out-pairs-ld`: SV–SNP pairs with calculated LD (r²) values. 

##### Disease/disorder-related trait detection
`match_traits.py`
- `--ref`: Reference list of disease and disorder terms from EMBL-EBI’s EFO ontology and SNOMED condition-domain vocabulary.
- `--gwas`: GWAS trait file containing traits identified from SV–SNP pairs in linkage disequilibrium.
- `--out`: Output file with traits matched to disease/disorder terms.
 
##### Summary of SV-SNP assocaitions
`assocaition_summary.py`
- `--all`: All identified SV–SNP pairs.
- `--genic`: Subset of pairs where SVs are located in trait-associated genes.
- `--disease`: Traits matched to disease/disorder terms
- `--out_fig`: Histogram showing the number of SVs in LD (r² ≥ 0.5) with GWAS SNPs.

#####  LD heatmap visualization
`Input_preparation.sh`: Prepare input genotype and variant information files for LD visualization.
`LD_heatmap_compuation.py`: Heatmap showing LD patterns (r²) between an SV and surrounding SNPs from the GWAS catalog.

### Note
EHR-related codes (e.g., PheWAS) are available in the All of Us Researcher Workbench. Related Jupyter notebooks (`.ipynb`) can be found in the [long-reads-public-codebase](https://github.com/all-of-us/long-reads-public-codebase/tree/main/notebooks/rw/) repository.
