#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genic", required=True, help="Genic SV–GWAS pairs.")
    parser.add_argument("--all", required=True, help="All SV–GWAS pairs.")
    parser.add_argument("--disease", required=True, help="disease/disorder-related traits.")
    parser.add_argument("--out-fig", default="Association_barplot.png")
    return parser.parse_args()

def main():
    args = parse_args()
    genic_file = args.genic
    all_file = args.all
    disease_file = args.disease
    out_fig = args.out

    with open(disease_file) as f:
        disease_traits = {line.strip() for line in f if line.strip()}

    genic_df = pd.read_csv(genic_file, sep="\t", header=None, names=["SV_ID", "GWAS_Info", "LD"])
    genic_df[["GWAS_variant", "GENE", "Trait", "Pval", "SNP_ID"]] = genic_df["GWAS_Info"].str.split("|", expand=True)
    genic_keys = set(zip(genic_df["SV_ID"], genic_df["Trait"]))

    all_df = pd.read_csv(all_file, sep="\t", header=None, names=["SV_ID", "GWAS_Info", "LD"])
    all_df[["GWAS_variant", "GENE", "Trait", "Pval", "SNP_ID"]] = all_df["GWAS_Info"].str.split("|", expand=True)
    all_df = all_df.drop_duplicates(subset=["SV_ID", "GWAS_variant", "Trait"])

    records = []
    for _, row in all_df.iterrows():
        try:
            ld = float(row["LD"])
        except (ValueError, TypeError):
            continue
        key = (row["SV_ID"], row["Trait"])
        category = "Genic" if key in genic_keys else "Intergenic"
        trait_type = "Disease" if row["Trait"] in disease_traits else "Trait"
        records.append((ld, category, trait_type))

    df = pd.DataFrame(records, columns=["LD", "Category", "TraitType"])
    bins = np.arange(0.5, 1.05, 0.05)
    df["LD_bin"] = pd.cut(df["LD"], bins=bins, include_lowest=True)
    counts = df.groupby(["Category", "LD_bin", "TraitType"]).size().unstack(fill_value=0)

    category_order = ["Intergenic", "Genic"]
    ld_bin_labels = counts.index.levels[1]

    color_map = {"Intergenic": "#1f77b4", "Genic": "#ff7f0e"}
    fig, ax = plt.subplots(figsize=(7.5, 4.9))
    bar_width = 0.35
    x = np.arange(len(ld_bin_labels))

    for i, category in enumerate(category_order):
        cat_data = counts.loc[category].reindex(ld_bin_labels, fill_value=0)
        disease_vals = cat_data.get("Disease", pd.Series(0, index=ld_bin_labels)).values
        trait_vals = cat_data.get("Trait", pd.Series(0, index=ld_bin_labels)).values
        xpos = x + (i - 0.5) * bar_width
        base_color = color_map[category]

        ax.bar(xpos, trait_vals, width=bar_width, color=base_color,
               edgecolor="black", linewidth=0.35)
        ax.bar(xpos, disease_vals, width=bar_width, bottom=trait_vals,
               color=base_color, hatch="///", edgecolor="black", linewidth=0.35)

    bin_labels = [f"{round(b.left,2)}–{round(b.right,2)}" for b in ld_bin_labels]
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, rotation=45)
    ax.set_xlabel("Linkage disequilibrium (r²)", fontsize=12)
    ax.set_ylabel("Number of SV–GWAS associations", fontsize=12)

    legend_elements = [
        Patch(facecolor=color_map["Genic"], edgecolor="black", label="SVs in trait-associated genes"),
        Patch(facecolor=color_map["Intergenic"], edgecolor="black", label="Other SVs"),
        Patch(facecolor="white", hatch="///", edgecolor="black", label="Disease/disorder-related", linewidth=0.5),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=12)

    plt.tight_layout()
    plt.savefig(out_fig, dpi=600)


if __name__ == "__main__":
    main()
