#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators

# Populations to process
populations = ["EAS", "SAS", "AFR", "AMR", "EUR"]

# Input files
variants_df = pd.read_csv("SV_summary.txt", sep="\t")
cds_df = pd.read_csv(
    "SV_CDS_ID.txt",
    sep="\t",
    header=None,
    names=["chr", "pos", "Variant_ID", "start", "end", "gene"]
)

genes_data = {}

for pop in populations:
    # Select variants present in that population
    selected_variants = variants_df[
        (variants_df["Sample_Count"] >= 1) &
        (variants_df["Populations"].str.contains(pop))
    ]
    selected_variant_ids = selected_variants["Variant_ID"].tolist()

    selected_genes = cds_df[cds_df["Variant_ID"].isin(selected_variant_ids)]

    genes_data[pop] = set(selected_genes["gene"])

all_genes = sorted(set().union(*genes_data.values()))

data = pd.DataFrame(
    {
        pop: [gene in genes_data[pop] for gene in all_genes]
        for pop in populations
    },
    index=all_genes
)

upset_input = from_indicators(populations, data[populations])

upset = UpSet(upset_input, subset_size='count', sort_by='cardinality')

fig = plt.figure(figsize=(11, 4.7), dpi=600)
upset.plot(fig=fig)

plt.ylabel("Number of genes")
plt.tight_layout()
plt.savefig("cds_genes_upset.png", dpi=600)
plt.savefig("cds_genes_upset.pdf", format="pdf", dpi=600)
plt.close(fig)
