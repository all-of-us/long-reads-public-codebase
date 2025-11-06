#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.transforms as transforms
from scipy.stats import pearsonr
import matplotlib

# Input files
snp_gt_file = "Selected_SNP_GT.txt"      
snp_pos_file = "Selected_SNP.bed"      
sv_gt_file  = "Selected_SV_GT.txt"    
sv_pos_file = "Selected_SV_pos" 
gene_file   = "Genes.bed"      # genes within the region, with columns: chr,start,end,gene_ID,strand

highlight_snp_file = "associated_snps.txt"  # six SNP IDs that show high LD with SVs.

chrom_label = "chr7"
x_min = 75410000
x_max = 75585000

# Output prefix (for LD matrix and figures)
out_prefix = "LD_heatmap_example"

# SV of interest
variant_of_interest = "chr7-75493495-allele831405-64"
variant_of_interest_label = f"{chrom_label}:75493495-INS-63"


def calculate_ld(genotype1, genotype2):
    valid_idx = (~np.isnan(genotype1)) & (~np.isnan(genotype2))
    if valid_idx.sum() > 1:
        r, _ = pearsonr(genotype1[valid_idx], genotype2[valid_idx])
        return r**2
    return 0.0

def main():
    snp_gt_df = pd.read_csv(snp_gt_file, sep="\t")
    snp_pos_df = pd.read_csv(
        snp_pos_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "info"],
    )

    snp_pos_df["GenotypeID"] = snp_pos_df["info"].str.split("|").str[0]
    snp_pos_df["RSID"]       = snp_pos_df["info"].str.split("|").str[1]
    snp_pos_df["start"]      = pd.to_numeric(snp_pos_df["start"], errors="coerce")

    snp_pos_mapping = snp_pos_df.set_index("GenotypeID")["start"].to_dict()
    id_to_rsid_map  = snp_pos_df.set_index("GenotypeID")["RSID"].to_dict()

    snp_gt_df["Position"] = snp_gt_df["Variant_ID"].map(snp_pos_mapping)
    snp_gt_df["Position"] = pd.to_numeric(snp_gt_df["Position"], errors="coerce")
    snp_gt_df.dropna(subset=["Position"], inplace=True)
    snp_gt_df.reset_index(drop=True, inplace=True)

    sv_gt_df = pd.read_csv(sv_gt_file, sep="\t")

    sv_pos_df = pd.read_csv(
        sv_pos_file,
        sep="\t",
        header=None,
        names=["chr", "pos", "Variant_ID"],
    )
    sv_pos_df["pos"] = pd.to_numeric(sv_pos_df["pos"], errors="coerce")
    sv_pos_mapping   = sv_pos_df.set_index("Variant_ID")["pos"].to_dict()

    sv_gt_df["Position"] = sv_gt_df["VariantID"].map(sv_pos_mapping)
    sv_gt_df["Position"] = pd.to_numeric(sv_gt_df["Position"], errors="coerce")
    sv_gt_df.dropna(subset=["Position"], inplace=True)
    sv_gt_df.reset_index(drop=True, inplace=True)

    common_samples = list(
        set(snp_gt_df.columns[1:-1]).intersection(sv_gt_df.columns[1:-1])
    )

    snp_gt_df = snp_gt_df[["Variant_ID", "Position"] + common_samples]
    sv_gt_df  = sv_gt_df[["VariantID",  "Position"] + common_samples]

    snp_gt_df.rename(columns={"Variant_ID": "ID"}, inplace=True)
    sv_gt_df.rename(columns={"VariantID":  "ID"}, inplace=True)

    all_variants_df = pd.concat([snp_gt_df, sv_gt_df], ignore_index=True)

    all_variants_df["Position"] = pd.to_numeric(
        all_variants_df["Position"], errors="coerce"
    )
    all_variants_df.dropna(subset=["Position"], inplace=True)
    all_variants_df.sort_values(by="Position", inplace=True)
    all_variants_df.reset_index(drop=True, inplace=True)

    all_variants_df.iloc[:, 2:] = (
        all_variants_df.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").fillna(0)
    )

    num_variants = len(all_variants_df)
    genotypes    = all_variants_df.iloc[:, 2:].values.astype(float)

    ld_matrix = np.zeros((num_variants, num_variants), dtype=float)

    for i in range(num_variants):
        for j in range(i, num_variants):
            ld_val = calculate_ld(genotypes[i], genotypes[j])
            ld_matrix[i, j] = ld_val
            ld_matrix[j, i] = ld_val

    # Save LD matrix
    ld_matrix_df = pd.DataFrame(
        ld_matrix,
        index=all_variants_df["ID"],
        columns=all_variants_df["ID"],
    )
    ld_matrix_out = f"{out_prefix}.LD_matrix.tsv"
    ld_matrix_df.to_csv(ld_matrix_out, sep="\t")

    gene_df = pd.read_csv(
        gene_file,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "gene", "strand"],
    )

    with open(highlight_snp_file) as f:
        selected_genotypes = [
            line.strip() for line in f
            if line.strip() and not line.startswith("#")
        ]

    fig = plt.figure(figsize=(10, 8))

    ld_ax = fig.add_axes([-0.1, 0.228, 1.12, 1.12])
    gene_ax = fig.add_axes([0.1, 0.05, 0.75, 0.18])
    cbar_ax = fig.add_axes([0.88, 0.6, 0.02, 0.25])

    tri_mask  = np.triu(np.ones_like(ld_matrix, dtype=bool), k=0)
    ld_masked = np.ma.masked_where(tri_mask, ld_matrix)

    cmap = matplotlib.cm.get_cmap("Reds").copy()
    cmap.set_bad(color="white")

    im = ld_ax.imshow(
        ld_masked,
        cmap=cmap,
        origin="lower",
        interpolation="nearest",
        extent=(0, num_variants, 0, num_variants),
        vmin=0,
        vmax=1,
    )

    cb = plt.colorbar(im, cax=cbar_ax)
    cb.set_label(r"$r^2$", fontsize=10)

    ld_ax.set_xticks([])
    ld_ax.set_yticks([])
    for spine in ld_ax.spines.values():
        spine.set_visible(False)

    # Rotate heatmap by -45 degrees
    rot_transform = transforms.Affine2D().rotate_deg(-45)
    im.set_transform(rot_transform + ld_ax.transData)

    ld_ax.set_xlim(0, 1.4 * num_variants)
    ld_ax.set_ylim(0, 1.4 * num_variants)

    for geno_id in selected_genotypes:
        idx_list = all_variants_df.index[all_variants_df["ID"] == geno_id].tolist()
        if not idx_list:
            print(f"Warning: {geno_id} not found in all_variants_df.")
            continue

        i = idx_list[0]
        display_label = id_to_rsid_map.get(geno_id, geno_id)

        col_rect = patches.Rectangle(
            (i, 0),
            1,
            num_variants,
            linewidth=1,
            edgecolor="black",
            facecolor="none",
        )
        col_rect.set_transform(rot_transform + ld_ax.transData)
        ld_ax.add_patch(col_rect)

        row_rect = patches.Rectangle(
            (0, i),
            num_variants,
            1,
            linewidth=1,
            edgecolor="black",
            facecolor="none",
        )
        row_rect.set_transform(rot_transform + ld_ax.transData)
        ld_ax.add_patch(row_rect)

        ld_ax.text(
            -3.2,
            i + 0.5,
            display_label,
            transform=(rot_transform + ld_ax.transData),
            color="black",
            fontsize=10,
            ha="center",
            va="center",
            rotation=-45,
        )

        ld_ax.text(
            i + 0.5,
            45.5,
            display_label,
            transform=(rot_transform + ld_ax.transData),
            color="black",
            fontsize=10,
            ha="center",
            va="center",
            rotation=45,
        )

    idx_list = all_variants_df.index[all_variants_df["ID"] == variant_of_interest].tolist()
    if idx_list:
        i = idx_list[0]
        col_rect = patches.Rectangle(
            (i, 0),
            1,
            num_variants,
            linewidth=1,
            edgecolor="#586cb2",
            facecolor="none",
        )
        col_rect.set_transform(rot_transform + ld_ax.transData)
        ld_ax.add_patch(col_rect)

        row_rect = patches.Rectangle(
            (0, i),
            num_variants,
            1,
            linewidth=1,
            edgecolor="#586cb2",
            facecolor="none",
        )
        row_rect.set_transform(rot_transform + ld_ax.transData)
        ld_ax.add_patch(row_rect)

        ld_ax.text(
            i + 0.5,
            i + 26.2,
            variant_of_interest_label,
            transform=(rot_transform + ld_ax.transData),
            color="#586cb2",
            fontsize=11,
            ha="center",
            va="center",
            rotation=45,
        )
        ld_ax.text(
            -6.2,
            i + 0.5,
            variant_of_interest_label,
            transform=(rot_transform + ld_ax.transData),
            color="#586cb2",
            fontsize=11,
            ha="center",
            va="center",
            rotation=-45,
        )

    gene_ax.set_xlim(x_min, x_max)
    gene_ax.set_ylim(0, 1)

    for spine in gene_ax.spines.values():
        spine.set_visible(False)

    chrom_line_y = 0.7
    gene_ax.hlines(
        y=chrom_line_y,
        xmin=x_min,
        xmax=x_max,
        colors="#2975aa",
        linewidth=2,
    )

    major_positions = [x_min, (x_min + x_max) // 2, x_max]
    for mp in major_positions:
        gene_ax.vlines(
            x=mp,
            ymin=chrom_line_y - 0.08,
            ymax=chrom_line_y,
            colors="black",
            linewidth=1,
        )
        if mp == x_min:
            label = f"{chrom_label}  {mp:,}"
        else:
            label = f"{mp:,}"
        gene_ax.text(
            mp,
            chrom_line_y - 0.10,
            label,
            ha="center",
            va="top",
            fontsize=12,
        )

    for _, row in gene_df.iterrows():
        start = row["start"]
        end   = row["end"]
        gene_name = row["gene"]
        strand    = row["strand"]
        gene_width = end - start

        # Vertical placement per gene (customized)
        if gene_name in ["POM121C", "HIP1"]:
            exon_center = 0.1
        elif gene_name in ["SPDYE5", "NSUN5P1", "PMS2P3"]:
            exon_center = 0.3
        else:
            continue

        exon_height = 0.1
        rect = patches.Rectangle(
            (start, exon_center - exon_height / 2),
            gene_width,
            exon_height,
            linewidth=0,
            facecolor="lightgray",
        )
        gene_ax.add_patch(rect)

        arrow_spacing  = 3500
        arrow_positions = np.arange(start + arrow_spacing / 1.5, end, arrow_spacing)
        arrow_text = ">" if strand == "+" else "<"
        for pos in arrow_positions:
            gene_ax.text(
                pos,
                exon_center,
                arrow_text,
                ha="center",
                va="center",
                fontsize=17,
                color="white",
            )

        if gene_name == "PMS2P3":
            gene_ax.text(
                end + 17000,
                exon_center - 0.01,
                gene_name,
                ha="right",
                va="center",
                fontsize=12,
            )
        else:
            gene_ax.text(
                start - 1500,
                exon_center - 0.01,
                gene_name,
                ha="right",
                va="center",
                fontsize=12,
            )
    gene_ax.axis("off")
    positions = all_variants_df["Position"].astype(int).values

    for i, pos in enumerate(positions):
        if pos < x_min or pos > x_max:
            continue

        bottom_data = np.array([pos, chrom_line_y])
        bottom_display = gene_ax.transData.transform(bottom_data)
        bottom_fig = fig.transFigure.inverted().transform(bottom_display)

        top_x = 0.029 + i / 47
        top_y = 0.225

        line = matplotlib.lines.Line2D(
            [bottom_fig[0], top_x],
            [bottom_fig[1], top_y],
            transform=fig.transFigure,
            color="gray",
            linewidth=1,
            linestyle="-",
        )
        fig.add_artist(line)

    png_out = f"{out_prefix}.png"
    pdf_out = f"{out_prefix}.pdf"

    plt.savefig(png_out, dpi=600)
    plt.savefig(pdf_out, format="pdf", dpi=400)


if __name__ == "__main__":
    main()
