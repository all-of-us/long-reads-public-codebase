#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches


def calculate_ld(genotype_sv, genotype_snp):
    try:
        r, _ = pearsonr(genotype_sv, genotype_snp)
        return r ** 2 if not np.isnan(r) else np.nan
    except Exception:
        return np.nan
      
def load_zscore_file(path, compressed=False):
    names = ["VariantID", "GeneID", "Zscore", "pscore", "pos"]
    df = pd.read_csv(path, sep="\t", header=None, compression="gzip" if compressed else None, names=names)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df["pscore"] = pd.to_numeric(df["pscore"], errors="coerce")
    return df.dropna(subset=["pos", "pscore"])

def load_genotypes(path, compressed=False):
    df = pd.read_csv(path, delim_whitespace=True, compression="gzip" if compressed else None)
    if df.columns[0] != "VariantID":
        df = df.rename(columns={df.columns[0]: "VariantID"})
    df.iloc[:, 1:] = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    return df

def genotypes_to_dict(df): return df.set_index("VariantID").T.to_dict("list")
def load_gene_bed(path): return pd.read_csv(path, sep="\t", header=None, names=["chr", "start", "end", "GeneID"])

def get_gene_region(gene_name, bed_df):
    row = bed_df[bed_df["GeneID"] == gene_name]
    return (int(row.iloc[0]["start"]), int(row.iloc[0]["end"])) if not row.empty else (None, None)

def add_enhancer(ax, start, end, y0=-0.2, height=0.73, alpha=0.3):
    if start and end and end > start:
        ax.add_patch(mpatches.Rectangle((start, y0), end - start, height, facecolor="green", alpha=alpha, edgecolor="none", zorder=2))

def add_exons_with_connectors(ax, bed_path, chrom_filter, x_min, x_max, y0=-0.2, height=0.73):
    try: exons = pd.read_csv(bed_path, sep=r"\s+", header=None, names=["chr", "start", "end"], comment="#")
    except Exception as e: print(f"Could not read exon BED file '{bed_path}': {e}"); return
    exons = exons[(exons["chr"] == chrom_filter) & (exons["end"] >= x_min) & (exons["start"] <= x_max)]
    if exons.empty: return
    exons.sort_values("start", inplace=True)
    center_y, coords = y0 + height * 0.5, []
    for _, r in exons.iterrows():
        ex_start, ex_end = max(int(r["start"]), x_min), min(int(r["end"]), x_max)
        if ex_end > ex_start:
            ax.add_patch(mpatches.Rectangle((ex_start, y0), ex_end - ex_start, height, facecolor="gray", edgecolor="black", linewidth=0.3, zorder=8))
            coords.append((ex_start, ex_end))
    for i in range(len(coords) - 1):
        prev_end, next_start = coords[i][1], coords[i + 1][0]
        if next_start > prev_end: ax.plot([prev_end, next_start], [center_y, center_y], linestyle="-", linewidth=0.3, color="black", zorder=9)


def create_manhattan_plot(gene_name, chrom_label, sv_df, snp_df, sv_gt_dict, snp_gt_dict, bed_df, exon_bed, enh_chr, enh_start, enh_end, x_center, half_window, out_png, out_pdf):
    gene_start, gene_end = get_gene_region(gene_name, bed_df)
    if not gene_start or not gene_end: print(f"Gene {gene_name} not found."); return
    x_min, x_max = x_center - half_window, x_center + half_window
    sv_f = sv_df[(sv_df["pos"] >= x_min) & (sv_df["pos"] <= x_max)]
    snp_f = snp_df[(snp_df["pos"] >= x_min) & (snp_df["pos"] <= x_max)]

    snp_to_ld = {}
    for _, snp_row in snp_f.iterrows():
        snp_id = snp_row["VariantID"]
        if snp_id not in snp_gt_dict: continue
        g_snp = snp_gt_dict[snp_id]
        for sv_id in sv_f["VariantID"]:
            if sv_id not in sv_gt_dict: continue
            g_sv = sv_gt_dict[sv_id]
            if len(g_sv) == len(g_snp):
                ld = calculate_ld(g_sv, g_snp)
                if 0 <= ld <= 1: snp_to_ld[snp_id] = max(ld, snp_to_ld.get(snp_id, 0.0))

    cmap = plt.cm.get_cmap("RdYlBu_r")
    fig = plt.figure(figsize=(6.5, 3.6), dpi=600)
    gs = gridspec.GridSpec(nrows=1, ncols=2, width_ratios=[20, 0.5], wspace=0.03)
    ax = fig.add_subplot(gs[0])

    for _, r in snp_f.iterrows(): ax.scatter(r["pos"], -np.log10(r["pscore"]), c=[snp_to_ld.get(r["VariantID"], 0)], cmap=cmap, vmin=0, vmax=1, marker=".", s=25, zorder=5)
    for _, r in sv_f.iterrows():  ax.scatter(r["pos"], -np.log10(r["pscore"]), c=[1], cmap=cmap, vmin=0, vmax=1, marker="*", s=60, zorder=6)

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlim(x_min, x_max)
    ax.set_xticks([x_min + (x_max - x_min) * f for f in (0.25, 0.5, 0.75)])
    ax.set_xlabel("Coordinate"); ax.set_ylabel("-log10(p-value)")
    ax.text(0.0, -0.06, chrom_label, fontsize=10.5, transform=ax.transAxes, ha="right")

    if enh_chr == chrom_label and (enh_start <= x_max and enh_end >= x_min):
        add_enhancer(ax, max(enh_start, x_min), min(enh_end, x_max))
    add_exons_with_connectors(ax, exon_bed, chrom_label, x_min, x_max)

    cax = fig.add_subplot(gs[1])
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1)), cax=cax)
    cbar.ax.tick_params(labelsize=10); cbar.set_label("r² to lead variant", fontsize=11)
    ax.scatter([], [], color="black", marker="*", s=30, label="SV"); ax.scatter([], [], color="black", marker=".", s=30, label="SNP")
    ax.legend(frameon=False, loc="upper right")
    plt.tight_layout(); plt.savefig(out_png, dpi=600); plt.savefig(out_pdf, dpi=600, format="pdf"); plt.close(fig)


sv_zscore_file = "eQTL_SV_zscore.BID.txt"
snp_zscore_file = "eQTL_SNP_zscore.BID.txt.gz"
sv_genotype_file = "SV_genotypes.BID.txt"
snp_genotype_file = "SNP_genotypes.BID.txt.gz"
gene_bed_file = "genes.bed.gz"
exon_bed_file = "BID_exons.bed"
gene_name = "ENSG00000015475.19"
chrom_label = "chr22"
enh_chr, enh_start, enh_end = "chr22", 17800733, 17801999
x_center, half_window = 17801130, 130000
out_png, out_pdf = "BID_manhattan_plot.png", "BID_manhattan_plot.pdf"

sv_df = load_zscore_file(sv_zscore_file)
snp_df = load_zscore_file(snp_zscore_file, compressed=True)
sv_gt = genotypes_to_dict(load_genotypes(sv_genotype_file))
snp_gt = genotypes_to_dict(load_genotypes(snp_genotype_file, compressed=True))
gene_bed = load_gene_bed(gene_bed_file)

create_manhattan_plot(gene_name, chrom_label, sv_df, snp_df, sv_gt, snp_gt, gene_bed, exon_bed_file, enh_chr, enh_start, enh_end, x_center, half_window, out_png, out_pdf)
