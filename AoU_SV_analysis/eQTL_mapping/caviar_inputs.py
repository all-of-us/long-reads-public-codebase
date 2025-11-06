#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd


def compute_ld_r2_matrix(geno_df):
    X = geno_df.to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        corr = np.corrcoef(X)
    corr = np.nan_to_num(corr, nan=0.0)
    r2 = corr ** 2
    np.fill_diagonal(r2, 1.0)
    return r2


def run(args):
    mage_zscores = pd.read_csv(
        args.sv_z,
        sep="\t",
        header=None,
        names=["VariantID", "GeneID", "Zscore", "pscore", "pos"]
    )
    kgp_zscores = pd.read_csv(
        args.snp_z,
        sep="\t",
        header=None,
        compression="gzip",
        names=["VariantID", "GeneID", "Zscore", "pscore", "pos"]
    )

    mage_genotypes = pd.read_csv(args.sv_gt, delim_whitespace=True)
    kgp_genotypes = pd.read_csv(args.snp_gt, delim_whitespace=True, compression="gzip")

    if mage_genotypes.columns[0] != "VariantID":
        mage_genotypes = mage_genotypes.rename(columns={mage_genotypes.columns[0]: "VariantID"})
    if kgp_genotypes.columns[0] != "VariantID":
        kgp_genotypes = kgp_genotypes.rename(columns={kgp_genotypes.columns[0]: "VariantID"})

    mage_genotypes.iloc[:, 1:] = mage_genotypes.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    kgp_genotypes.iloc[:, 1:] = kgp_genotypes.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

    all_genes = sorted(
        set(mage_zscores["GeneID"].dropna()) |
        set(kgp_zscores["GeneID"].dropna())
    )

    os.makedirs(args.out_dir, exist_ok=True)

    for gene in all_genes:
        gene_folder = os.path.join(args.out_dir, gene)
        os.makedirs(gene_folder, exist_ok=True)

        mage_gene = mage_zscores.loc[mage_zscores["GeneID"] == gene, ["VariantID", "Zscore", "pscore"]]
        kgp_gene = kgp_zscores.loc[kgp_zscores["GeneID"] == gene, ["VariantID", "Zscore", "pscore"]]

        if kgp_gene.empty and mage_gene.empty:
            continue

        if len(kgp_gene) > args.top_kgp_snps:
            kgp_gene = kgp_gene.nsmallest(args.top_kgp_snps, "pscore")

        kgp_ids = list(kgp_gene["VariantID"])
        mage_ids = [vid for vid in mage_gene["VariantID"] if vid not in set(kgp_ids)]
        ordered_ids = kgp_ids + mage_ids

        if len(ordered_ids) == 0:
            continue

        mage_gene_gt = mage_genotypes[mage_genotypes["VariantID"].isin(mage_gene["VariantID"])]
        kgp_gene_gt = kgp_genotypes[kgp_genotypes["VariantID"].isin(kgp_gene["VariantID"])]

        gt_combined = pd.concat([kgp_gene_gt, mage_gene_gt], axis=0, ignore_index=True)
        gt_combined = gt_combined.drop_duplicates(subset=["VariantID"], keep="first")
        gt_combined = gt_combined.set_index("VariantID")
        gt_combined = gt_combined.reindex(ordered_ids)

        z_map = pd.concat(
            [
                kgp_gene[["VariantID", "Zscore"]],
                mage_gene[["VariantID", "Zscore"]],
            ],
            axis=0,
        ).drop_duplicates("VariantID", keep="first").set_index("VariantID")

        z_in_order = pd.DataFrame(
            {
                "VariantID": ordered_ids,
                "Zscore": z_map.reindex(ordered_ids)["Zscore"].astype(float).values,
            }
        )

        ld_r2 = compute_ld_r2_matrix(gt_combined)

        zscore_file = os.path.join(gene_folder, "variant.zscore")
        ld_file = os.path.join(gene_folder, "variant.ld")
        order_file = os.path.join(gene_folder, "variant.list")

        z_in_order.to_csv(zscore_file, sep=" ", index=False, header=False)
        np.savetxt(ld_file, ld_r2, fmt="%.6f")
        pd.Series(ordered_ids).to_csv(order_file, index=False, header=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-z",dest="sv_zscore",default="eQTL_SV_zscore.txt")
    parser.add_argument("--sv-gt",dest="sv_genotype",default="SV_genotypes.txt")
    parser.add_argument("--snp-z",dest="snp_zscore",default="eQTL_SNP_zscore.txt.gz")
    parser.add_argument("--snp-gt", dest="snp_genotype", default="SNP_genotypes.txt.gz")
    parser.add_argument("--out-dir", dest="out_dir", default="caviar/")
    parser.add_argument("--top-kgp-snps", dest="top_kgp_snps", type=int,default=1000)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
