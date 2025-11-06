#!/usr/bin/env python3
import argparse
import csv
import numpy as np
import pandas as pd


def calculate_ld(sv_gt, snp_gt):
    sv_centered = sv_gt - np.mean(sv_gt)
    snp_centered = snp_gt - np.mean(snp_gt)
    covariance = np.sum(sv_centered * snp_centered)
    var_sv = np.sum(sv_centered ** 2)
    var_snp = np.sum(snp_centered ** 2)
    if var_sv == 0 or var_snp == 0:
        return np.nan
    return (covariance ** 2) / (var_sv * var_snp)


def load_genotypes(path, index_col):
    df = pd.read_csv(path, sep="\t")
    df.set_index(index_col, inplace=True)
    return df


def load_gwas_pairs(path):
    pairs = set()
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            sv_id = str(parts[0])
            snp_id = parts[1].split("|")[0]
            pairs.add((sv_id, snp_id))
    return pairs


def compute_ld_dict(sv_genotypes, snp_genotypes, gwas_pairs):
    ld_dict = {}
    for sv_id, snp_id in gwas_pairs:
        if sv_id in sv_genotypes.index and snp_id in snp_genotypes.index:
            sv_gt = sv_genotypes.loc[sv_id].values.astype(float)
            snp_gt = snp_genotypes.loc[snp_id].values.astype(float)
            ld = calculate_ld(sv_gt, snp_gt)
            if not np.isnan(ld):
                ld_dict[(sv_id, snp_id)] = ld
        else:
            print(f"Missing genotype for pair: {sv_id}\t{snp_id}")
    return ld_dict


def filter_pairs_by_ld(pairs_file, ld_dict, out_file):
    with open(pairs_file, "r") as in_f, open(out_file, "w", newline="") as out_f:
        reader = csv.reader(in_f, delimiter="\t")
        writer = csv.writer(out_f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            sv_id = row[0]
            snp_id = row[1].split("|")[0]
            ld_value = ld_dict.get((sv_id, snp_id), np.nan)
            if not np.isnan(ld_value) and ld_value >= 0.5:
                writer.writerow(row + [f"{ld_value}"])


def run(args):
    sv_genotypes = load_genotypes(args.sv_gt, index_col="VariantID")
    snp_genotypes = load_genotypes(args.snp_gt, index_col="Variant_ID")
    common_samples = sorted(sv_genotypes.columns.intersection(snp_genotypes.columns))
    sv_genotypes = sv_genotypes[common_samples]
    snp_genotypes = snp_genotypes[common_samples]
    gwas_pairs = load_gwas_pairs(args.pairs)
    ld_dict = compute_ld_dict(sv_genotypes, snp_genotypes, gwas_pairs)
    filter_pairs_by_ld(args.pairs, ld_dict, args.out_pairs_ld)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sv-gt", dest="sv_gt", default="SV_genotype.tsv")
    parser.add_argument("--snp-gt", dest="snp_gt", default="SNP_genotype.tsv.gz")
    parser.add_argument("--pairs", dest="pairs", default="SVs_near_GWAS_SNPs_100kb.txt")
    parser.add_argument("--out-pairs-ld", dest="out_pairs_ld", default="1KGP_GWAS_IDs_LD_100K.txt")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
