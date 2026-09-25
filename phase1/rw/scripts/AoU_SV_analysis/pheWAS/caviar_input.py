r"""Write CAVIAR inputs for each SV-condition association: z-scores and an LD (r) matrix
for the SV and the top SNPs at the same locus. README.md section 5, step 7.

Inputs:
  --sv-pheWAS-file  tab-separated with columns SV, Disease, Z_disc (the SV's z-score)
  --snp-dir         per-locus SNP association results named
                    <SV>_<Disease>.<Disease>.glm.firth[.gz] (PLINK2 Firth output with
                    ID, P, A1, ALT, Z_STAT columns). These come from a separate SNP
                    association run at each locus, not from association_testing.sh.
  --ld-dir          <SV>_<Disease>.vcor files from compute_ld_r.sh
  --top-n           keep at most this many SNPs per locus, by smallest p-value (default 1000)
SNP z-scores are oriented to the ALT allele.
Outputs, in --out-dir/<SV>_<Disease>/:
  variant.zscore, variant.list  the SV first, then the SNPs
  variant.ld                    r matrix in the same order (pairs without LD set to 0)
Run run_caviar.sh afterwards.

Usage:
  python caviar_input.py --sv-pheWAS-file sv_disease.txt --snp-dir snp_assoc/ \
      --ld-dir ld/ --top-n 1000 --out-dir caviar_phewas/
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd


def find_file(patterns):
    matches = []

    for pattern in patterns:
        matches.extend(glob.glob(pattern))

    matches = sorted(set(matches))
    return matches[0] if matches else None


def find_glm_file(snp_dir, sv, disease_name):
    patterns = [
        os.path.join(snp_dir, f"{sv}_{disease_name}.{disease_name}.glm.firth"),
        os.path.join(snp_dir, f"{sv}_{disease_name}.{disease_name}.glm.firth.gz"),
        os.path.join(snp_dir, f"{sv}_{disease_name}*.glm.firth"),
        os.path.join(snp_dir, f"{sv}_{disease_name}*.glm.firth.gz"),
    ]

    return find_file(patterns)


def find_ld_file(ld_dir, sv, disease_name):
    patterns = [
        os.path.join(ld_dir, f"{sv}_{disease_name}.vcor"),
        os.path.join(ld_dir, f"{sv}_{disease_name}*.vcor"),
    ]

    return find_file(patterns)


def read_glm_top_snps(glm_file, top_n):
    df = pd.read_csv(
        glm_file,
        sep=r"\s+",
        engine="python",
        dtype=str,
        compression="infer"
    )

    required_cols = {"ID", "P", "A1", "ALT", "Z_STAT"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns in {glm_file}: {missing}")

    if "ERRCODE" in df.columns:
        df = df[df["ERRCODE"] == "."].copy()

    df["P"] = pd.to_numeric(df["P"], errors="coerce")
    df["Z_STAT"] = pd.to_numeric(df["Z_STAT"], errors="coerce")

    df = df[df["P"].notna()].copy()
    df = df[df["Z_STAT"].notna()].copy()

    df["zscore"] = df["Z_STAT"]

    # Orient SNP z-scores to the ALT allele (PLINK reports them for A1).
    flip_mask = df["A1"].astype(str) != df["ALT"].astype(str)
    df.loc[flip_mask, "zscore"] = -df.loc[flip_mask, "zscore"]

    df = df.sort_values("P", ascending=True).copy()

    if len(df) > top_n:
        df = df.head(top_n).copy()

    return df[["ID", "zscore", "P", "A1", "ALT", "Z_STAT"]].copy()


def build_ld_matrix(ld_file, variants):
    n = len(variants)
    variant_to_idx = {v: i for i, v in enumerate(variants)}

    mat = np.zeros((n, n), dtype=float)
    np.fill_diagonal(mat, 1.0)

    ld_df = pd.read_csv(
        ld_file,
        sep=r"\s+",
        engine="python",
        dtype=str
    )

    required_cols = {"ID_A", "ID_B", "UNPHASED_R"}
    missing = required_cols - set(ld_df.columns)

    if missing:
        raise ValueError(f"Missing columns in {ld_file}: {missing}")

    ld_df["UNPHASED_R"] = pd.to_numeric(
        ld_df["UNPHASED_R"],
        errors="coerce"
    )

    n_pairs_used = 0

    for _, row in ld_df.iterrows():
        var_a = row["ID_A"]
        var_b = row["ID_B"]
        r = row["UNPHASED_R"]

        if pd.isna(r):
            continue

        if var_a in variant_to_idx and var_b in variant_to_idx:
            i = variant_to_idx[var_a]
            j = variant_to_idx[var_b]

            mat[i, j] = float(r)
            mat[j, i] = float(r)

            n_pairs_used += 1

    return mat, n_pairs_used


def run(args):
    os.makedirs(args.out_dir, exist_ok=True)

    sv_df = pd.read_csv(args.sv_pheWAS_file, sep="\t", dtype=str)

    required_cols = {"SV", "Disease", "Z_disc"}
    missing = required_cols - set(sv_df.columns)

    if missing:
        raise ValueError(f"Missing columns in {args.sv_pheWAS_file}: {missing}")

    sv_df["Z_disc"] = pd.to_numeric(sv_df["Z_disc"], errors="coerce")

    n_done = 0
    n_skip = 0

    for _, row in sv_df.iterrows():
        sv = str(row["SV"]).strip()
        disease_name = str(row["Disease"]).strip()
        sv_z = row["Z_disc"]

        if not sv or not disease_name:
            n_skip += 1
            continue

        if pd.isna(sv_z):
            n_skip += 1
            continue

        prefix = f"{sv}_{disease_name}"
        out_dir = os.path.join(args.out_dir, prefix)
        os.makedirs(out_dir, exist_ok=True)

        glm_file = find_glm_file(args.snp_dir, sv, disease_name)
        ld_file = find_ld_file(args.ld_dir, sv, disease_name)

        if glm_file is None:
            n_skip += 1
            continue

        if ld_file is None:
            n_skip += 1
            continue

        try:
            snp_df = read_glm_top_snps(glm_file, top_n=args.top_n)
        except Exception as e:
            n_skip += 1
            continue

        snp_df = snp_df[snp_df["ID"] != sv].copy()

        if snp_df.empty:
            n_skip += 1
            continue

        # The SV is always the first variant.
        variants = [sv] + snp_df["ID"].tolist()
        zscores = [sv_z] + snp_df["zscore"].tolist()

        zscore_out = os.path.join(out_dir, "variant.zscore")

        zscore_df = pd.DataFrame({
            "variant": variants,
            "zscore": zscores
        })

        zscore_df.to_csv(
            zscore_out,
            sep=" ",
            header=False,
            index=False
        )

        variant_list_out = os.path.join(out_dir, "variant.list")

        pd.Series(variants).to_csv(
            variant_list_out,
            index=False,
            header=False
        )

        try:
            ld_mat, n_ld_pairs_used = build_ld_matrix(ld_file, variants)
        except Exception as e:
            n_skip += 1
            continue

        if ld_mat.shape[0] != len(variants) or ld_mat.shape[1] != len(variants):
            raise ValueError(
                f"LD matrix shape {ld_mat.shape} does not match "
                f"number of variants {len(variants)} for {sv} / {disease_name}"
            )

        expected_pairs = len(variants) * (len(variants) - 1) // 2

        if n_ld_pairs_used < expected_pairs:
            print(
                f"[WARN] {sv} / {disease_name}: LD pairs used = {n_ld_pairs_used}, "
                f"expected = {expected_pairs}. Missing pairs were set to 0."
            )

        ld_out = os.path.join(out_dir, "variant.ld")

        np.savetxt(
            ld_out,
            ld_mat,
            fmt="%.10g",
            delimiter="\t"
        )

        print(
            f"[DONE] {sv} / {disease_name}: "
            f"{len(variants)} variants written to {out_dir}"
        )

        n_done += 1


def main():
    parser = argparse.ArgumentParser(
        description="Prepare CAVIAR input files for PheWAS fine-mapping."
    )

    parser.add_argument(
        "--sv-pheWAS-file",
        default="SV_pheWAS_summary.txt"
    )

    parser.add_argument(
        "--snp-dir",
        default="snp_disease_plink",
    )

    parser.add_argument(
        "--ld-dir",
        default="LD_r"
    )

    parser.add_argument(
        "--out-dir",
        default="caviar_input"
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=1000
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
