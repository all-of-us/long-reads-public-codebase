"""Extract standard SNOMED disorder diagnoses for the study participants from the
All of Us EHR. README.md section 5, step 1; Methods, "Phenome-wide SV analysis".

Runs in the Researcher Workbench: queries the condition_occurrence table of the CDR
named by the WORKSPACE_CDR environment variable through BigQuery.
Inputs (fixed names, current directory):
  sample_IDs.txt                    All of Us person_id values, one per line
  athena_ohdsi_SNOMED/CONCEPT.csv   OHDSI Athena vocabulary export (tab-separated)
Output:
  AoU_samples_EHR_disorder.csv  one row per condition record: person_id, concept_id,
                                concept_name, concept_class_id, standard_concept.
                                Used as --ehr-file for phenotype.py.

Usage:
  python extract_AoU_condition.py
"""
import os
import pandas as pd
from google.cloud import bigquery

sample_ids_file = "sample_IDs.txt"
vocabulary_file = "athena_ohdsi_SNOMED/CONCEPT.csv"
output_file = "AoU_samples_EHR_disorder.csv"
CDR = os.environ["WORKSPACE_CDR"]
client = bigquery.Client()

with open(sample_ids_file, "r") as f:
    sample_ids = [
        int(line.strip())
        for line in f
        if line.strip()
    ]

query = f"""
SELECT
    co.person_id,
    co.condition_concept_id AS concept_id
FROM `{CDR}.condition_occurrence` AS co
WHERE co.person_id IN UNNEST(@sample_ids)
  AND co.condition_concept_id IS NOT NULL
  AND co.condition_concept_id != 0
"""

job_config = bigquery.QueryJobConfig(
    query_parameters=[
        bigquery.ArrayQueryParameter(
            "sample_ids",
            "INT64",
            sample_ids
        )
    ]
)

df_conditions = client.query(query, job_config=job_config).to_dataframe(
    create_bqstorage_client=("BIGQUERY_STORAGE_API_ENABLED" in os.environ)
)


df_vocab = pd.read_csv(
    vocabulary_file,
    sep="\t",
    usecols=[
        "concept_id",
        "concept_name",
        "concept_class_id",
        "standard_concept"
    ],
    low_memory=False
)

# Make concept_id type consistent
df_conditions["concept_id"] = pd.to_numeric(
    df_conditions["concept_id"],
    errors="coerce"
).astype("Int64")

df_vocab["concept_id"] = pd.to_numeric(
    df_vocab["concept_id"],
    errors="coerce"
).astype("Int64")

df_merged = pd.merge(
    df_conditions,
    df_vocab,
    on="concept_id",
    how="left"
)


# Keep standard (S) SNOMED concepts of class Disorder.
df_merged = df_merged[
    (df_merged["concept_class_id"] == "Disorder") &
    (df_merged["standard_concept"] == "S")
].copy()

df_merged.to_csv(output_file, index=False)
