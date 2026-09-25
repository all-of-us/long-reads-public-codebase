#!/usr/bin/env bash

# For each SV-condition pair, compute pairwise LD (signed r) among the SV and its
# candidate SNPs, for CAVIAR fine-mapping. README.md section 5, step 7.
#
# Inputs:
#   variant_disease_pairs.tsv.gz  tab-separated: SV ID, condition, then one or more SNP IDs
#   AoU_SV_SNPs.pgen/.pvar/.psam  PLINK2 fileset containing both the SVs and the SNPs
# Output:
#   LD_r/<SV>_<condition>.vcor    plink2 --r-unphased output (existing files are skipped)
# Requires plink2.
#
# Usage: bash compute_ld_r.sh

PAIR_FILE="variant_disease_pairs.tsv.gz"
PFILE_PREFIX="AoU_SV_SNPs"

OUTDIR="LD_r"
THREADS=6

mkdir -p "${OUTDIR}"

TMP_PAIR_LIST=$(mktemp)
TMP_SNP_LIST=$(mktemp)
TMP_VAR_LIST=$(mktemp)

trap 'rm -f "${TMP_PAIR_LIST}" "${TMP_SNP_LIST}" "${TMP_VAR_LIST}"' EXIT

# Unique (SV, condition) pairs.
gzip -cd "${PAIR_FILE}" | awk -F'\t' '
NF >= 3 {
    print $1 "\t" $2
}
' | sort -u > "${TMP_PAIR_LIST}"

while IFS=$'\t' read -r SV disease; do

    echo "[INFO] SV=${SV}; disease=${disease}"

    gzip -cd "${PAIR_FILE}" | awk -F'\t' -v sv="${SV}" -v dis="${disease}" '
    NF >= 3 && $1 == sv && $2 == dis {
        for (i = 3; i <= NF; i++) {
            if ($i != "") print $i
        }
    }
    ' | sort -u > "${TMP_SNP_LIST}"

    if [[ ! -s "${TMP_SNP_LIST}" ]]; then
        echo "[WARN] No SNPs found for ${SV} / ${disease}" >&2
        continue
    fi

    {
        echo "${SV}"
        cat "${TMP_SNP_LIST}"
    } | sort -u > "${TMP_VAR_LIST}"

    PREFIX="${OUTDIR}/${SV}_${disease}"

    if [[ -f "${PREFIX}.vcor" ]]; then
        echo "[SKIP] LD already exists: ${PREFIX}.vcor"
        continue
    fi

    # Subset to the SV and its SNPs, then compute all pairwise r.
    plink2 \
        --pfile "${PFILE_PREFIX}" \
        --extract "${TMP_VAR_LIST}" \
        --make-pgen \
        --threads "${THREADS}" \
        --out "${PREFIX}.subset"

    plink2 \
        --pfile "${PREFIX}.subset" \
        --r-unphased ref-based \
        --ld-window-kb 999999999 \
        --ld-window-r2 0 \
        --threads "${THREADS}" \
        --out "${PREFIX}"

    rm -f "${PREFIX}.subset.pgen" \
          "${PREFIX}.subset.pvar" \
          "${PREFIX}.subset.psam"

done < "${TMP_PAIR_LIST}"
