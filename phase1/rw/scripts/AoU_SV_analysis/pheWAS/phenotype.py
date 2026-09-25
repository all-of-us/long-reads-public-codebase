r"""Build the PLINK2 binary phenotype file for the PheWAS.

One column per selected condition: case (default 2) if the sample has any EHR record
of that concept, otherwise control (default 1). README.md section 5, step 4.

Inputs:
  --disease-file  condition concept IDs, one per line (these become the column names)
  --sample-file   sample IDs, first whitespace-separated column
  --ehr-file      CSV with columns person_id and concept_id
                  (for example AoU_samples_EHR_disorder.csv from extract_AoU_condition.py)
Output:
  --out-pheno     tab-separated: FID, IID, one column per condition

Usage:
  python phenotype.py --disease-file disease_concepts.txt --sample-file sample_IDs.txt \
      --ehr-file conditions.tsv --out-pheno pheno.txt --fid-value 0 --case-value 2 --control-value 1
"""
import argparse
import pandas as pd


def normalize_sample_id(x):
    return str(x).strip()


def clean_pheno_name(x):
    return str(x).strip().replace(" ", "_")


def load_selected_diseases(path):
    diseases = []

    with open(path) as f:
        for line in f:
            disease = line.strip()
            if disease:
                diseases.append(disease)

    diseases = list(dict.fromkeys(diseases))

    return diseases


def load_samples(path, fid_value="0"):
    samples = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # If the sample file has multiple columns, use the first column as IID.
            iid = line.split()[0]
            samples.append(normalize_sample_id(iid))

    sample_df = pd.DataFrame({
        "FID": fid_value,
        "IID": samples
    })

    n_dup = sample_df["IID"].duplicated().sum()
    if n_dup > 0:
        print(f"[Warning] Found {n_dup} duplicated sample IDs. Keeping first occurrence.")
        sample_df = sample_df.drop_duplicates("IID", keep="first")

    return sample_df


def run(args):
    selected_diseases = load_selected_diseases(args.disease_file)
    selected_disease_set = set(selected_diseases)

    print(f"Loaded {len(selected_diseases)} selected diseases")

    sample_df = load_samples(args.sample_file, fid_value=args.fid_value)
    print(f"Loaded {len(sample_df)} target samples")

    ehr_df = pd.read_csv(args.ehr_file, dtype=str)

    required_cols = {"person_id", "concept_id"}
    missing = required_cols - set(ehr_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.ehr_file}: {missing}")

    ehr_df = ehr_df[["person_id", "concept_id"]].copy()
    ehr_df["person_id"] = ehr_df["person_id"].apply(normalize_sample_id)
    ehr_df["concept_id"] = ehr_df["concept_id"].astype(str).str.strip()

    # Keep records of the selected conditions in the target samples.
    ehr_df = ehr_df[ehr_df["concept_id"].isin(selected_disease_set)].copy()
    print(f"Selected disease records retained: {len(ehr_df)}")

    target_samples = set(sample_df["IID"])
    ehr_df = ehr_df[ehr_df["person_id"].isin(target_samples)].copy()
    print(f"Selected disease records in target samples: {len(ehr_df)}")

    ehr_df["phenotype"] = ehr_df["concept_id"].apply(clean_pheno_name)

    ehr_df = ehr_df.drop_duplicates(subset=["person_id", "phenotype"]).copy()

    selected_pheno_cols = [clean_pheno_name(x) for x in selected_diseases]
    selected_pheno_cols = list(dict.fromkeys(selected_pheno_cols))

    if ehr_df.empty:
        case_table = pd.DataFrame(index=[], columns=selected_pheno_cols)
    else:
        case_table = (
            ehr_df.assign(value=args.case_value)
            .pivot(index="person_id", columns="phenotype", values="value")
        )

    pheno_df = sample_df.merge(
        case_table,
        left_on="IID",
        right_index=True,
        how="left"
    )

    for col in selected_pheno_cols:
        if col not in pheno_df.columns:
            pheno_df[col] = pd.NA

    final_cols = ["FID", "IID"] + selected_pheno_cols
    pheno_df = pheno_df[final_cols].copy()

    # Samples with no record of a condition are controls.
    pheno_df[selected_pheno_cols] = (
        pheno_df[selected_pheno_cols]
        .fillna(args.control_value)
        .astype(int)
    )

    pheno_df.to_csv(args.out_pheno, sep="\t", index=False)

    case_counts = (pheno_df[selected_pheno_cols] == args.case_value).sum(axis=0)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare PLINK binary phenotype file."
    )

    parser.add_argument(
        "--disease-file",
        default="disease_IDs.txt"
    )

    parser.add_argument(
        "--sample-file",
        default="sample_IDs.txt"
    )

    parser.add_argument(
        "--ehr-file",
        default="samples_EHRs.csv"
    )

    parser.add_argument(
        "--out-pheno",
        default="pheno.txt"
    )

    parser.add_argument(
        "--fid-value",
        default="0"
    )

    parser.add_argument(
        "--case-value",
        type=int,
        default=2
    )

    parser.add_argument(
        "--control-value",
        type=int,
        default=1
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
