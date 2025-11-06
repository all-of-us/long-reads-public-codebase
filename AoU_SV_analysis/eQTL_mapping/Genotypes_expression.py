import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gene', required=True)
parser.add_argument('--SV', required=True)
args = parser.parse_args()

gene = args.gene
SV = args.SV
significant_eGenes_df = pd.read_csv('eQTL_result.bh.csv')
genotypes_df = pd.read_csv('genotypes.csv', index_col=0)

expression_df = pd.read_csv('TMM_expression.csv', index_col=0)
gene_name_df = pd.read_csv('gene_name', sep='\t', header=None, names=['gene_id', 'gene_name'])
gene_name_mapping = pd.Series(gene_name_df.gene_name.values, index=gene_name_df.gene_id).to_dict()

def convert_gt(gt):
    if gt in ['0/0', '0|0']:
        return 0
    elif gt in ['0/1', '1/0', '1|0', '0|1']:
        return 1
    elif gt in ['1/1', '1|1']:
        return 2
    else:
        return 3

genotypes_df = genotypes_df.applymap(convert_gt)
variant_gt = genotypes_df.loc[str(SV)]
gene_expression = expression_df.loc[gene]
plot_data = {
    'GT': variant_gt,
    'Expression': gene_expression
}
plot_df = pd.DataFrame(plot_data).dropna()
color_map = {0: 'green', 1: 'orange', 2: 'blue'}
positions = [0, 0.7, 1.4]
fig, ax = plt.subplots(figsize=(6.5, 5))
for gt, group_df in plot_df.groupby('GT'):
    pos = positions[gt]
    ax.boxplot(group_df['Expression'], widths=0.4, positions=[gt], patch_artist=False, showfliers=False)
    jitter = np.random.normal(0, 0.04, size=len(group_df))
    ax.scatter(gt + jitter, group_df['Expression'], color=color_map[gt], alpha=1, s=3)

bh_pval = significant_eGenes_df.loc[significant_eGenes_df['gene'] == gene, 'bh_corrected_pval'].values[0]
gene_name = gene_name_mapping.get(gene, gene)
ax.set_title(f'{gene_name}', fontsize=12)
ax.set_xticks([0, 1, 2])
ax.tick_params(axis='y', labelsize=12)
ax.set_xticklabels(['AA', 'AB', 'BB'],fontsize=12)
ax.set_xlabel('Genotype',fontsize=12)
ax.set_ylabel('Expression',fontsize=12)
plt.savefig(SV + "_" + gene + '.png', dpi=600)
