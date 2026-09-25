#!/bin/bash

# Convert the imputed SV VCF into the tab-separated genotype table read by the
# population-scale scripts (README.md section 2, step 1).
#
# Input:  SV_genotype.vcf.gz  KAGE + GLIMPSE imputed SV callset on 1000 Genomes samples
# Output: SV_genotype.stats.gz  header CHROM POS ID REF ALT <one column per sample>,
#         then one row per SV with each sample's GT string (for example 0|1)
# Requires bcftools.
#
# Usage: bash stats.sh

sample_names=$(bcftools query -l SV_genotype.vcf.gz| tr '\n' '\t')
echo -e "CHROM\tPOS\tID\tREF\tALT\t${sample_names}" > SV_genotype.stats
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' SV_genotype.vcf.gz >> SV_genotype.stats
gzip SV_genotype.stats
