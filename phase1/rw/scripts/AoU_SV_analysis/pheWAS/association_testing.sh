#!/usr/bin/env bash

# PheWAS: PLINK2 Firth logistic regression of every SV against each condition.
# README.md section 5, step 6.
#
# Inputs:
#   AoU_SV.pgen/.pvar/.psam  SV genotypes (plink2 --make-pgen of the SV VCF)
#   pheno.txt                from phenotype.py
#   covar.txt                from covar.py; covariates age, sex, PC1-PC10
#                            (variance-standardized). For a sex-specific condition,
#                            remove sex from --covar-name below.
#   disease.txt              phenotype column names to test, one per line
# Output:
#   <disease>.<disease>.glm.firth.gz for each condition
# Runs up to MAX_JOBS conditions in parallel, THREADS threads each. Requires plink2 and
# bash 4.3 or later (wait -n).
#
# Usage: bash association_testing.sh

DISEASE_FILE="disease.txt"
MAX_JOBS=10
THREADS=4

# Test one condition; PLINK2 writes <disease>.<disease>.glm.firth.
run_one() {
    local disease="$1"

    echo "[START] ${disease}"

    plink2 \
        --pfile AoU_SV \
        --pheno pheno.txt \
        --pheno-name "${disease}" \
        --covar covar.txt \
        --covar-name age sex PC1-PC10 \
        --covar-variance-standardize \
        --glm firth hide-covar \
        --threads "${THREADS}" \
        --out "${disease}"

    if [[ -f "${disease}.${disease}.glm.firth" ]]; then
        gzip -f "${disease}.${disease}.glm.firth"
        echo "[DONE] ${disease}"
    else
        echo "[WARN] ${disease}.${disease}.glm.firth not found for ${disease}" >&2
    fi
}

# Launch one job per condition, keeping at most MAX_JOBS running.
job_count=0

while IFS= read -r disease || [[ -n "$disease" ]]; do
    [[ -z "$disease" ]] && continue
    run_one "$disease" &
    ((job_count+=1))
    if (( job_count >= MAX_JOBS )); then
        wait -n
        ((job_count-=1))
    fi
done < "$DISEASE_FILE"
wait
echo "All jobs finished."
