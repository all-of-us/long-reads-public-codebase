#!/usr/bin/env bash
# Extract SNP and SV genotypes around one SV for the LD heatmap example
# (LD_heatmap_compuation.py). README.md section 4, step 5.
#
# Set up for SV chr7-75493495-allele831405-64 and a +/- 100 kb window; edit the
# variables below for another locus.
# Inputs:
#   GWAS_snps.bed.gz     GWAS Catalog SNPs: chrom, start, end,
#                        "variantID|MappedGene|Trait|P-value|RiskAllele"
#   SNP_genotype.tsv.gz  SNP genotypes: header line, then one row per SNP
#   SV_pos               SV positions: CHROM, POS, SV_ID (with a header line)
#   SV_genotype.tsv      SV genotypes: header line, then one row per SV
# Outputs:
#   snp_ids.chr7_75493495_100kb.txt  IDs of the GWAS SNPs in the window
#   Selected_SNP_genotype.txt        their genotypes
#   Selected_SNP.bed                 their positions, labelled "variantID|rsID"
#   Selected_SV_pos, Selected_SV_genotype.txt  position and genotypes of the SV
# Requires bedtools.
#
# Usage: bash Input_preparation.sh

CHR="chr7"
CENTER=75493495
WINDOW=100000
# Input files
SNPS_BED="GWAS_snps.bed.gz"
SNP_GT="SNP_genotype.tsv.gz"
SV_POS="SV_pos" #CHROM   POS     SV_ID
SV_GT="SV_genotype.tsv"

# Outputs
SNP_ID_LIST="snp_ids.chr7_75493495_100kb.txt"
FILTERED_SNP_GT="Selected_SNP_genotype.txt"
FILTERED_SNP_BED="Selected_SNP.bed"
SELECTED_SV_POS="Selected_SV_pos"
FILTERED_SV_GT="Selected_SV_genotype.txt"

# SV of interest
SV_OF_INTEREST="chr7-75493495-allele831405-64"

# 1. IDs of the GWAS SNPs inside the window.
zcat "$SNPS_BED" \
  | awk -v chr="$CHR" -v center="$CENTER" -v w="$WINDOW" '
      $1 == chr && $2 > center - w && $3 <= center + w { print $4 }
    ' \
  | awk -F"|" '{print $1}' \
  | sort -u \
  > "$SNP_ID_LIST"

# 2. Genotypes of those SNPs (header line kept).
zcat "$SNP_GT" | head -n 1 > "$FILTERED_SNP_GT"
zcat "$SNP_GT"  | grep -w -Ff "$SNP_ID_LIST" -   >> "$FILTERED_SNP_GT"


# 3. BED of those SNPs, labelled "variantID|rsID" (allele suffix removed).
zcat "$SNPS_BED" \
  | grep -w -Ff "$SNP_ID_LIST" - \
  | sort -u \
  | bedtools sort -i - \
  | awk -F"\t" 'BEGIN{OFS="\t"}{
        split($4, a, "|");
        # a[1] = VariantID (e.g., 7:pos:ref:alt)
        # a[5] = rsID (based on your original command)
        print $1, $2, $3, a[1] "|" a[5]
    }' \
  | awk -F"-" '{print $1}' \
  | sort -u \
  | bedtools sort -i - \
  > "$FILTERED_SNP_BED"


# 4. Position and genotypes of the SV of interest.
head -n 1 "$SV_POS" > "$SELECTED_SV_POS"
grep -w "$SV_OF_INTEREST" "$SV_POS" >> "$SELECTED_SV_POS"


head -n 1 "$SV_GT" > "$FILTERED_SV_GT"
awk '{print $3}' "$SELECTED_SV_POS" \
  | grep -w -Ff - "$SV_GT" \
  >> "$FILTERED_SV_GT"
