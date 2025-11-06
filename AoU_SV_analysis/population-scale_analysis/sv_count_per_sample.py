import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import pandas as pd

# Input and output file paths
input_file = 'SV_genotype.stats.gz'
output_file = 'SV_per_sample.txt'
population_file = '1KG_population' # Sample_name\tSuperpopulation_code,from https://www.internationalgenome.org/data/

population_df = pd.read_csv(population_file, delim_whitespace=True)
population_dict = dict(zip(population_df['Sample_name'], population_df['Superpopulation_code']))

with gzip.open(input_file, 'rt') as f:
    header = f.readline().strip().split('\t')
    sample_names = header[5:]

    variant_counts = {sample: {'Hom': 0, 'Het': 0} for sample in sample_names}

    for line in f:
        fields = line.strip().split('\t')
        for i, gt in enumerate(fields[5:], start=0):
            if '1' in gt:
                if gt in ['0/1', '1/0','1|0', '0|1']:
                    variant_counts[sample_names[i]]['Het'] += 1
                elif gt in ['1/1','1|1']:
                    variant_counts[sample_names[i]]['Hom'] += 1

with open(output_file, 'w') as out:
    out.write("Sample\tPopulation\tHom\tHet\n")  # Write header
    for sample, counts in variant_counts.items():
        population = population_dict.get(sample, 'Unknown') 
        out.write(f"{sample}\t{population}\t{counts['Hom']}\t{counts['Het']}\n")


df = pd.read_csv(output_file, delim_whitespace=True)

population_order = ['AFR', 'AMR',  'EUR', 'EAS', 'SAS']

df['Population'] = pd.Categorical(df['Population'], categories=population_order, ordered=True)
df = df.sort_values(by='Population')

df['Total_SVs'] = df['Hom'] + df['Het']

plt.figure(figsize=(5.5, 4))
sns.violinplot(x='Population', y='Total_SVs', data=df, order=population_order)

plt.ylabel('Number of SVs', fontsize=12)

plt.xlabel('Population', fontsize=12)

plt.tight_layout()
plt.savefig("1KG_spl_SV_violin.png", dpi=600)
plt.savefig("1KG_spl_SV_violin.pdf", dpi=600, format='pdf')
