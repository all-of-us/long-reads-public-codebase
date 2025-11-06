#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main():
    parser = argparse.ArgumentParser(description="CADD-SV score vs number of samples.")
    parser.add_argument("--cadd", required=True, help="Path to CADD BED file")
    parser.add_argument("--summary", required=True, help="Path to SV summary file")
    parser.add_argument("--known_ids", required=True, help="Path to SV IDs file")
    parser.add_argument("--output", default="Cadd_sv_sample_summary.png", help="Output PNG filename")
    args = parser.parse_args()

    df_cadd_raw = pd.read_csv(args.cadd, sep="\t", header=None, names=["chr", "start", "end", "info"])
    def parse_cadd_info(x):
        parts = x.split("|")
        sv_id = parts[0].strip() if len(parts) > 0 else None
        sv_score = float(parts[2].strip()) if len(parts) > 2 else None
        return pd.Series([sv_id, sv_score])
    df_cadd_raw[["VariantID", "CADD_Score"]] = df_cadd_raw["info"].apply(parse_cadd_info)
    df_cadd = df_cadd_raw[["VariantID", "CADD_Score"]].dropna()

    df_summary = pd.read_csv(args.summary, sep="\t")
    df_summary["VariantID"] = df_summary["Variant_ID"].astype(str).str.strip()
    df_summary["SampleCount"] = df_summary["Sample_Count"].astype(int)
    df_summary = df_summary[["VariantID", "SampleCount"]]

    with open(args.known_ids, "r") as f:
        known_ids = set(line.strip() for line in f if line.strip())
    df_merged = pd.merge(df_cadd, df_summary, on="VariantID", how="left")
    df_merged.dropna(subset=["SampleCount"], inplace=True)
    df_merged["Novelty"] = df_merged["VariantID"].apply(lambda x: "Yes" if x in known_ids else "No")

    plt.figure(figsize=(8.2, 4.9))
    colors = {"Yes": "orange", "No": "blue"}
    novelty_order = ["Yes", "No"]

    for novelty in novelty_order:
        subset = df_merged[df_merged["Novelty"] == novelty]
        plt.scatter(
            subset["SampleCount"],
            subset["CADD_Score"],
            color=colors[novelty],
            marker="o",
            alpha=0.7,
            s=10,
        )

    plt.xlabel("Number of Samples", fontsize=13)
    plt.ylabel("CADD-SV Score", fontsize=13)
    plt.axhline(y=20, color="gray", linestyle="--", linewidth=0.8)

    x_min, x_max = df_merged["SampleCount"].min(), df_merged["SampleCount"].max()
    plt.xlim(left=x_min - 5, right=x_max + 5)
    plt.tick_params(axis='both', labelsize=12)
    novelty_handles = [
        plt.Line2D([0], [0], marker="o", color="w", label="Known",
                   markerfacecolor="orange", markersize=11),
        plt.Line2D([0], [0], marker="o", color="w", label="AoU detected",
                   markerfacecolor="blue", markersize=11),
    ]
    plt.legend(handles=novelty_handles, title="SV", title_fontsize=14, fontsize=14,
               loc="upper right", frameon=False)

    plt.tight_layout()
    plt.savefig(args.output, dpi=600)

if __name__ == "__main__":
    main()
