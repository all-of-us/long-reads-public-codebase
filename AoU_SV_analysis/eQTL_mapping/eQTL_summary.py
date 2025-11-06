#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text


def plot_eqtl(args):
    eqtl_bh = pd.read_csv(args.bh_file)
    eqtl = pd.read_csv(args.eqtl_file)

    eqtl_bh["SV"] = eqtl_bh["SV"].astype(str)
    eqtl["SV"] = eqtl["SV"].astype(str)
    eqtl["gene_id"] = eqtl["gene"].str.split(".").str[0]
    eqtl.drop(columns=["gene"], inplace=True)

    filtered_bh = eqtl_bh[eqtl_bh["bh_corrected_pval"] <= 0.05].copy()
    filtered_bh["gene_id"] = filtered_bh["gene"].str.split(".").str[0]
    filtered_bh = filtered_bh[["gene_id", "SV", "bh_corrected_pval"]]

    merged = pd.merge(filtered_bh, eqtl, on=["gene_id", "SV"], how="left")
    merged["-log_bh_corrected_pval"] = -np.log10(merged["bh_corrected_pval"])

    df_gene_map = pd.read_csv(args.gene_name_file, sep="\t",
                              header=None, names=["gene_id", "gene_name"])
    df_gene_map["gene_id"] = df_gene_map["gene_id"].str.split(".").str[0]
    merged = pd.merge(merged, df_gene_map, on="gene_id", how="left")

    omim_genes = set()
    with open(args.gene_list_file, "r") as f:
        for line in f:
            gene = line.strip().split()[0] if line.strip() else None
            if gene:
                omim_genes.add(gene)

    plt.figure(figsize=(5.5, 4.07))
    neg_data = merged[merged["coef"] < 0]
    pos_data = merged[merged["coef"] >= 0]

    plt.scatter(neg_data["coef"], neg_data["-log_bh_corrected_pval"],
                color="blue", alpha=0.7, s=11)
    plt.scatter(pos_data["coef"], pos_data["-log_bh_corrected_pval"],
                color="red", alpha=0.7, s=11)

    plt.xlabel("Effect Size (β)", fontsize=11)
    plt.ylabel("-log10(q-val)", fontsize=11)
    plt.xlim(-1.9, 1.9)

    auto_omim_genes = set(
        merged.loc[
            (merged["gene_name"].isin(omim_genes)) &
            (merged["-log_bh_corrected_pval"] >= 60),
            "gene_name"
        ].dropna()
    )

    texts, points = [], []
    for _, row in merged.iterrows():
        if row["gene_name"] in auto_omim_genes:
            x, y = row["coef"], row["-log_bh_corrected_pval"]
            txt = plt.text(x, y, row["gene_name"], fontsize=8,
                           ha="center", va="bottom")
            texts.append(txt)
            points.append((x, y))

    adjust_text(
        texts,
        arrowprops=None,
        force_points=0.3,
        force_text=0.3,
        only_move={'text': 'y', 'points': 'y'},
        expand_points=(0.6, 0.6),
        expand_text=(0.6, 0.6)
    )

    threshold = 7
    for txt, (x_pt, y_pt) in zip(texts, points):
        x_txt, y_txt = txt.get_position()
        if np.hypot(x_txt - x_pt, y_txt - y_pt) > threshold:
            plt.plot([x_pt, x_txt], [y_pt, y_txt],
                     color='gray', linewidth=1)

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=600)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eqtl_bh", dest="eqtl_bh_results", default="eQTL_result.bh.csv")
    parser.add_argument("--eqtl", dest="eqtl_results",default="eQTL_result.csv")
    parser.add_argument("--gene-name", dest="gene_name_file",default="gene_name")
    parser.add_argument("--gene-list", dest="gene_list_file",default="OMIM_genes")
    parser.add_argument("--out-png", dest="out_png",default="eQTL_summary.png")
    args = parser.parse_args()

    plot_eqtl(args)


if __name__ == "__main__":
    main()
