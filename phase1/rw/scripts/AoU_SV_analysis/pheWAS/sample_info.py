#!/usr/bin/env python3
"""Download demographics for All of Us participants with EHR data: gender, date of
birth, race, ethnicity, and sex at birth. README.md section 5, step 2.

Runs in the Researcher Workbench, against the CDR named by the WORKSPACE_CDR
environment variable (pandas.read_gbq). Covers every participant with EHR data;
covar.py later keeps only the study samples.
Output:
  AoU_samples_info.csv  used as --info-file for covar.py

Usage:
  python sample_info.py
"""

import os
import pandas as pd

cdr = os.environ["WORKSPACE_CDR"]

sample_info_sql = f"""
SELECT
    person.person_id,
    person.gender_concept_id,
    p_gender_concept.concept_name AS gender,
    person.birth_datetime AS date_of_birth,
    person.race_concept_id,
    p_race_concept.concept_name AS race,
    person.ethnicity_concept_id,
    p_ethnicity_concept.concept_name AS ethnicity,
    person.sex_at_birth_concept_id,
    p_sex_at_birth_concept.concept_name AS sex_at_birth
FROM `{cdr}.person` person
LEFT JOIN `{cdr}.concept` p_gender_concept
    ON person.gender_concept_id = p_gender_concept.concept_id
LEFT JOIN `{cdr}.concept` p_race_concept
    ON person.race_concept_id = p_race_concept.concept_id
LEFT JOIN `{cdr}.concept` p_ethnicity_concept
    ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id
LEFT JOIN `{cdr}.concept` p_sex_at_birth_concept
    ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id
WHERE person.person_id IN (
    SELECT DISTINCT person_id
    FROM `{cdr}.cb_search_person`
    WHERE has_ehr_data = 1
)
"""

sample_info_df = pd.read_gbq(
    sample_info_sql,
    dialect="standard",
    use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    progress_bar_type="tqdm_notebook"
)

sample_info_df.to_csv("AoU_samples_info.csv", index=False)
