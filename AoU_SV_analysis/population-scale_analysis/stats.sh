#!/bin/bash

sample_names=$(bcftools query -l SV_genotype.vcf.gz| tr '\n' '\t')
echo -e "CHROM\tPOS\tID\tREF\tALT\t${sample_names}" > SV_genotype.stats
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%GT\t]\n' SV_genotype.vcf.gz >> SV_genotype.stats
gzip SV_genotype.stats
