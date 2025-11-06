import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fst_file = 'AFR_vs_nonAFR.weir.fst'
fst_data = pd.read_csv(fst_file, sep='\t') 

sv_id_file = 'SV_IDs' # SV IDs in the same order as in the imputation VCF (SV_genotype.vcf.gz)
sv_ids = pd.read_csv(sv_id_file)

fst_data['SV_ID'] = sv_ids['SV_ID']
fst_data = fst_data[fst_data['WEIR_AND_COCKERHAM_FST'].notna()]

af_file = 'SV_AF_AFR_vs_NonAFR.txt'
af_data = pd.read_csv(af_file, sep='\t')

merged_data = af_data.merge(fst_data[['SV_ID', 'WEIR_AND_COCKERHAM_FST']], on='SV_ID')
bin_edges = np.arange(0, 0.8 + 0.05, 0.05)

afr_greater = []
afr_lesser = []
for i in range(len(bin_edges) - 1):
    bin_start = bin_edges[i]
    bin_end = bin_edges[i + 1]
    in_bin = (merged_data['WEIR_AND_COCKERHAM_FST'] > bin_start) & (merged_data['WEIR_AND_COCKERHAM_FST'] <= bin_end)
    afr_greater_count = ((in_bin) & (merged_data['AFR_AF'] > merged_data['Non_AFR_AF'])).sum()
    afr_lesser_count = ((in_bin) & (merged_data['AFR_AF'] <= merged_data['Non_AFR_AF'])).sum()
    afr_greater.append(afr_greater_count)
    afr_lesser.append(afr_lesser_count)

bar_width = 0.018
group_spacing = 0.02
x_positions = np.arange(len(bin_edges) - 1) * (bar_width * 2 + group_spacing)
afr_greater_pos = x_positions - bar_width / 2
afr_lesser_pos = x_positions + bar_width / 2

fst_gt_0_25 = merged_data[merged_data['WEIR_AND_COCKERHAM_FST'] > 0.15]
bins_inset = np.arange(0.15, 0.8 + 0.01, 0.05)
afr_greater_inset = []
afr_lesser_inset = []
for i in range(len(bins_inset) - 1):
    bin_start = bins_inset[i]
    bin_end = bins_inset[i + 1]
    in_bin = (fst_gt_0_25['WEIR_AND_COCKERHAM_FST'] > bin_start) & (fst_gt_0_25['WEIR_AND_COCKERHAM_FST'] <= bin_end)
    afr_greater_inset.append(((in_bin) & (fst_gt_0_25['AFR_AF'] > fst_gt_0_25['Non_AFR_AF'])).sum())
    afr_lesser_inset.append(((in_bin) & (fst_gt_0_25['AFR_AF'] <= fst_gt_0_25['Non_AFR_AF'])).sum())

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(afr_greater_pos, afr_greater, width=bar_width, edgecolor='black', color='skyblue', label='AFR AF > Non-AFR AF')
ax.bar(afr_lesser_pos, afr_lesser, width=bar_width, edgecolor='black', color='coral', label='AFR AF ≤ Non-AFR AF')

ax.set_xticks(x_positions)
ax.set_xticklabels([f'({bin_edges[i]:.2f}, {bin_edges[i+1]:.2f}]' for i in range(len(bin_edges) - 1)], rotation=45, fontsize=10)
ax.set_xlabel('Fixation index', fontsize=12)
ax.set_ylabel('Number of SVs', fontsize=12)

inset_ax = fig.add_axes([0.45, 0.53, 0.5, 0.4])
inset_x_positions = np.arange(len(bins_inset) - 1) * (bar_width * 2 + group_spacing)
inset_afr_greater_pos = inset_x_positions - bar_width / 2
inset_afr_lesser_pos = inset_x_positions + bar_width / 2
inset_ax.bar(inset_afr_greater_pos, afr_greater_inset, width=bar_width, edgecolor='black', color='skyblue', label='AFR AF > Non-AFR AF')
inset_ax.bar(inset_afr_lesser_pos, afr_lesser_inset, width=bar_width, edgecolor='black', color='coral', label='AFR AF ≤ Non-AFR AF')
inset_ax.set_xticks(inset_x_positions)
inset_ax.set_xticklabels([f'({bins_inset[i]:.2f}, {bins_inset[i+1]:.2f}]' for i in range(len(bins_inset) - 1)], rotation=45, fontsize=8)
inset_ax.legend(loc='upper right', fontsize=8, frameon=False)

plt.tight_layout()
plt.savefig("fst_histogram.png", dpi=600)
