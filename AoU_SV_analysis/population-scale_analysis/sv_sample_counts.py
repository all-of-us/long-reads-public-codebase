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
        matching_samples = [sample_names[i] for i, gt in enumerate(genotypes) if '1' in gt]
        count = len(matching_samples)
        matching_populations = {population_dict.get(sample, 'Unknown') for sample in matching_samples}
        matching_populations_str = ",".join(matching_populations)
        if count > 0:
            variant_counts.append((variant_id, count, ",".join(matching_samples), matching_populations_str))

df = pd.DataFrame(variant_counts, columns=['Variant_ID',  'Sample_Count', 'Sample_IDs', 'Populations'])
df.to_csv(output_file, sep='\t', index=False)
