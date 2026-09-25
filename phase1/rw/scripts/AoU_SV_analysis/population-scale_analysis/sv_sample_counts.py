"""For each SV, list the 1000 Genomes samples that carry it and their superpopulations.

README.md section 2, step 2. Reads and writes fixed file names in the current directory.

Inputs:
  SV_genotype.stats.gz  genotype table from stats.sh
  1KG_population        tab-separated; columns Sample_name, Superpopulation_code
                        (from https://www.internationalgenome.org/data/)
Output:
  SV_summary.txt        tab-separated; Variant_ID, Sample_Count, Sample_IDs and
                        Populations (both comma-separated). SVs carried by no sample
                        are omitted. Read by Gene_upset.py and regulatory_element_upset.py.

Usage:
  python sv_sample_counts.py
"""
import gzip
import pandas as pd

input_file = 'SV_genotype.stats.gz'
output_file = 'SV_summary.txt'
population_file = '1KG_population' # Sample_name\tSuperpopulation_code,from https://www.internationalgenome.org/data/

population_df = pd.read_csv(population_file, sep='\t')
population_dict = dict(zip(population_df['Sample_name'], population_df['Superpopulation_code']))

with gzip.open(input_file, 'rt') as f:
    header = f.readline().strip().split('\t')
    sample_names = header[5:]
    variant_counts = []
    line_number = 0
    for line in f:
        fields = line.strip().split('\t')
        variant_id = fields[2]
        genotypes = fields[5:]
        # A sample carries the SV if its GT string contains a 1 (het or hom-alt).
        matching_samples = [sample_names[i] for i, gt in enumerate(genotypes) if '1' in gt]
        count = len(matching_samples)
        matching_populations = {population_dict.get(sample, 'Unknown') for sample in matching_samples}
        matching_populations_str = ",".join(matching_populations)
        if count > 0:
            variant_counts.append((variant_id, count, ",".join(matching_samples), matching_populations_str))

df = pd.DataFrame(variant_counts, columns=['Variant_ID',  'Sample_Count', 'Sample_IDs', 'Populations'])
df.to_csv(output_file, sep='\t', index=False)
