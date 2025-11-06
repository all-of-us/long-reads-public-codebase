import gzip
import pandas as pd

input_file = 'SV_genotype.stats.gz'
output_file = 'SV_AF_AFR_vs_NonAFR.txt'
population_file = '1KG_population'
population_df = pd.read_csv(population_file, sep='\t')
population_dict = dict(zip(population_df['Sample_name'], population_df['Superpopulation_code']))

target_populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

with gzip.open(input_file, 'rt') as f:
    header = f.readline().strip().split('\t')
    sample_names = header[5:]
    variant_allele_frequencies = []
    for line in f:
        fields = line.strip().split('\t')
        variant_id = fields[2]
        genotypes = fields[5:]
        population_counts = {pop: {'1_count': 0, 'total_count': 0} for pop in target_populations}

        for sample_name, genotype in zip(sample_names, genotypes):
            population = population_dict.get(sample_name, 'Unknown')
            if population in target_populations:
                num_ones = genotype.count('1')
                population_counts[population]['1_count'] += num_ones
                population_counts[population]['total_count'] += 2
        afr_counts = population_counts['AFR']
        non_afr_counts = {
            '1_count': sum(population_counts[pop]['1_count'] for pop in target_populations if pop != 'AFR'),
            'total_count': sum(population_counts[pop]['total_count'] for pop in target_populations if pop != 'AFR')
        }
        print ( afr_counts['1_count'], afr_counts['total_count'], non_afr_counts['1_count'] , non_afr_counts['total_count'], variant_id)
        afr_af = afr_counts['1_count'] / afr_counts['total_count'] if afr_counts['total_count'] > 0 else 'NA'
        non_afr_af = non_afr_counts['1_count'] / non_afr_counts['total_count'] if non_afr_counts['total_count'] > 0 else 'NA'
        row = [variant_id, afr_af, non_afr_af]
        variant_allele_frequencies.append(row)
      
df = pd.DataFrame(variant_allele_frequencies, columns=['SV_ID', 'AFR_AF', 'Non_AFR_AF'])
df.to_csv(output_file, sep='\t', index=False)
