#!/usr/bin/env python3
r"""Build the PLINK2 covariate file for the PheWAS: FID, IID, age, sex, PC1..PCn.

README.md section 5, step 5.

Inputs:
  --sample-file  sample IDs, first whitespace-separated column
  --pc-file      tab-separated with columns research_id and pca_features (a list of
                 principal-component values per sample)
  --info-file    CSV with columns person_id, date_of_birth, sex_at_birth
                 (for example AoU_samples_info.csv from sample_info.py)
Age is computed at --reference-date (default 2025-10-01). Sex is coded 1 = male,
2 = female, and unknown as 0 or NA (--unknown-sex).
Output:
  --out-covar    tab-separated; missing values written as NA

Usage:
  python covar.py --sample-file sample_IDs.txt --pc-file pcs.txt --info-file sample_info.tsv \
      --out-covar covar.txt --reference-date 2020-01-01 --n-pcs 10
"""

import argparse
import ast
import re
import pandas as pd


SEX_MAP = {
    "male": 1,
    "female": 2,
    "m": 1,
    "f": 2,
    "1": 1,
    "1.0": 1,
    "2": 2,
    "2.0": 2,
}


def normalize_sample_id(x):
    return str(x).strip()


def load_samples(path):
    samples = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # If the file has one column, use it as IID.
            # If it has multiple columns, still use the first column by default.
            iid = line.split()[0]
            samples.append(normalize_sample_id(iid))

    sample_df = pd.DataFrame({"IID": samples})

    n_dup = sample_df["IID"].duplicated().sum()
    if n_dup > 0:
        print(f"[Warning] Found {n_dup} duplicated sample IDs in sample file. Keeping first occurrence.")
        sample_df = sample_df.drop_duplicates("IID", keep="first")

    return sample_df


def parse_pc_list(pc_str):
    if pd.isna(pc_str):
        return []

    s = str(pc_str).strip()
    if not s:
        return []

    try:
        vals = ast.literal_eval(s)
        if isinstance(vals, (list, tuple)):
            return list(vals)
    except Exception:
        pass

    nums = re.findall(r"[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", s)
    return nums


def extract_pc_columns(pc_str, n_pc):
    vals = parse_pc_list(pc_str)

    out = {}
    for i in range(n_pc):
        key = f"PC{i + 1}"
        out[key] = vals[i] if i < len(vals) else pd.NA

    return out


def load_pc_file(path, n_pc):
    pc_df = pd.read_csv(path, sep="\t", dtype=str)

    required = {"research_id", "pca_features"}
    missing = required - set(pc_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    pc_df = pc_df[["research_id", "pca_features"]].copy()
    pc_df["IID"] = pc_df["research_id"].apply(normalize_sample_id)

    n_dup = pc_df["IID"].duplicated().sum()
    if n_dup > 0:
        print(f"[Warning] Found {n_dup} duplicated sample IDs in PC file. Keeping first occurrence.")
        pc_df = pc_df.drop_duplicates("IID", keep="first")

    pc_expanded = pc_df["pca_features"].apply(
        lambda x: pd.Series(extract_pc_columns(x, n_pc=n_pc))
    )

    pc_df = pd.concat([pc_df[["IID"]], pc_expanded], axis=1)

    for i in range(1, n_pc + 1):
        pc_df[f"PC{i}"] = pd.to_numeric(pc_df[f"PC{i}"], errors="coerce")

    return pc_df


def compute_age(dob_series, ref_date):
    dob = pd.to_datetime(dob_series, errors="coerce", utc=True)

    age = ref_date.year - dob.dt.year

    before_birthday = (
        (ref_date.month < dob.dt.month)
        | ((ref_date.month == dob.dt.month) & (ref_date.day < dob.dt.day))
    )

    age = age - before_birthday.astype(float)

    return age


def normalize_sex_value(v, unknown_sex_value="0"):
    if pd.isna(v):
        return pd.NA if unknown_sex_value == "NA" else 0

    s = str(v).strip().lower()

    if s in SEX_MAP:
        return SEX_MAP[s]

    return pd.NA if unknown_sex_value == "NA" else 0


def load_info_file(path, ref_date, unknown_sex_value):
    info_df = pd.read_csv(path, dtype=str)

    required = {"person_id", "date_of_birth", "sex_at_birth"}
    missing = required - set(info_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {missing}")

    info_df = info_df[["person_id", "date_of_birth", "sex_at_birth"]].copy()
    info_df["IID"] = info_df["person_id"].apply(normalize_sample_id)

    n_dup = info_df["IID"].duplicated().sum()
    if n_dup > 0:
        print(f"[Warning] Found {n_dup} duplicated sample IDs in info file. Keeping first occurrence.")
        info_df = info_df.drop_duplicates("IID", keep="first")

    info_df["age"] = compute_age(info_df["date_of_birth"], ref_date)
    info_df["sex"] = info_df["sex_at_birth"].apply(
        lambda x: normalize_sex_value(x, unknown_sex_value=unknown_sex_value)
    )

    info_df["age"] = pd.to_numeric(info_df["age"], errors="coerce")
    info_df["sex"] = pd.to_numeric(info_df["sex"], errors="coerce")

    return info_df[["IID", "age", "sex"]]


def run(args):
    ref_date = pd.Timestamp(args.reference_date)

    sample_df = load_samples(args.sample_file)
    print(f"Loaded {len(sample_df)} target samples")

    pc_df = load_pc_file(args.pc_file, n_pc=args.n_pcs)
    print(f"Loaded PC info for {len(pc_df)} samples")

    info_df = load_info_file(
        args.info_file,
        ref_date=ref_date,
        unknown_sex_value=args.unknown_sex,
    )

    # One row per target sample; missing demographics or PCs stay NA.
    covar_df = sample_df.merge(info_df, on="IID", how="left")
    covar_df = covar_df.merge(pc_df, on="IID", how="left")

    covar_df.insert(0, "FID", args.fid_value)

    pc_cols = [f"PC{i}" for i in range(1, args.n_pcs + 1)]
    final_cols = ["FID", "IID", "age", "sex"] + pc_cols
    covar_df = covar_df[final_cols].copy()

    covar_df["age"] = pd.to_numeric(covar_df["age"], errors="coerce")
    covar_df["sex"] = pd.to_numeric(covar_df["sex"], errors="coerce")

    for pc in pc_cols:
        covar_df[pc] = pd.to_numeric(covar_df[pc], errors="coerce")

    covar_df.to_csv(args.out_covar, sep="\t", index=False, na_rep="NA")

def main():
    parser = argparse.ArgumentParser(
        description="Prepare covariate file for PLINK Firth regression."
    )

    parser.add_argument(
        "--sample-file",
        default="sample_IDs",
    )

    parser.add_argument(
        "--pc-file",
        default="ancestry_preds.tsv",
    )

    parser.add_argument(
        "--info-file",
        default="All_EHR_samples_info.csv",
    )

    parser.add_argument(
        "--out-covar",
        default="covar.txt"
    )

    parser.add_argument(
        "--reference-date",
        default="2025-10-01"
    )

    parser.add_argument(
        "--n-pcs",
        type=int,
        default=10
    )

    parser.add_argument(
        "--fid-value",
        default="0"
    )

    parser.add_argument(
        "--unknown-sex",
        choices=["0", "NA"],
        default="0"
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
