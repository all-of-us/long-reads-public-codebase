#!/usr/bin/env bash
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

zcat "$SNPS_BED" \
  | awk -v chr="$CHR" -v center="$CENTER" -v w="$WINDOW" '
      $1 == chr && $2 > center - w && $3 <= center + w { print $4 }
    ' \
  | awk -F"|" '{print $1}' \
  | sort -u \
  > "$SNP_ID_LIST"

zcat "$SNP_GT" | head -n 1 > "$FILTERED_SNP_GT"
zcat "$SNP_GT"  | grep -w -Ff "$SNP_ID_LIST" -   >> "$FILTERED_SNP_GT"


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


head -n 1 "$SV_POS" > "$SELECTED_SV_POS"
grep -w "$SV_OF_INTEREST" "$SV_POS" >> "$SELECTED_SV_POS"


head -n 1 "$SV_GT" > "$FILTERED_SV_GT"
awk '{print $3}' "$SELECTED_SV_POS" \
  | grep -w -Ff - "$SV_GT" \
  >> "$FILTERED_SV_GT"
