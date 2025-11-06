#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

def normalize_sample_id(s: str) -> str:
    return str(s).strip()

def first_token(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip()
    if not x:
        return np.nan
    return x.split()[0]

def convert_gt_cell(gt):
    if pd.isna(gt):
        return np.nan
    g0 = first_token(gt)
    if g0 in ("0/0", "0|0"):
        return 0.0
    elif g0 in ("0/1", "1/0", "1|0", "0|1"):
        return 1.0
    elif g0 in ("1/1", "1|1"):
        return 2.0
    else:
        return np.nan

def ensure_numeric_df(df: pd.DataFrame) -> pd.DataFrame:
    return df.apply(pd.to_numeric, errors='coerce')

def autodetect_read(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine='python')

def normalize_sex_value(v):
    if pd.isna(v):
        return np.nan
    s = str(v).strip()
    if s in ('1', '1.0'):
        return 1.0
    if s in ('2', '2.0'):
        return 2.0
    try:
        x = float(s)
        if x in (1.0, 2.0):
            return x
    except Exception:
        pass
    return np.nan

def meets_criteria(row: pd.Series) -> bool:
    counts = row.value_counts(dropna=True)
    return all(counts.get(x, 0) >= 37 for x in (0.0, 1.0, 2.0))

def run_eqtl(geno_file, expr_file, covar_file, pairs_file, out_pair_results, out_gene_bh):
    genotypes_raw = pd.read_csv(geno_file, index_col=0, dtype=str)
    geno_cols_norm = [normalize_sample_id(c) for c in genotypes_raw.columns]
    genotypes_raw.columns = geno_cols_norm

    expression_df = pd.read_csv(expr_file, index_col=0, dtype=str)
    expr_cols_norm = [normalize_sample_id(c) for c in expression_df.columns]
    _, idx_first = np.unique(expr_cols_norm, return_index=True)
    expression_df = expression_df.iloc[:, sorted(idx_first)]
    expr_cols_norm = [normalize_sample_id(c) for c in expression_df.columns]
    expression_df.columns = expr_cols_norm
    expression_df = expression_df.applymap(first_token)
    expression_df = ensure_numeric_df(expression_df)

    covar_wide = pd.read_csv(covar_file, sep='\t', dtype=str)
    covariates_df = covar_wide.set_index(covar_wide.columns[0]).T
    covariates_df.index = [normalize_sample_id(i) for i in covariates_df.index]
    covariates_df = covariates_df[['sex']].copy()
    covariates_df['sex'] = covariates_df['sex'].map(normalize_sex_value)
    covariates_df = covariates_df.dropna(axis=0, how='any')
    covariates_df = ensure_numeric_df(covariates_df)
    var_mask_global = covariates_df.var(axis=0, skipna=True) > 0
    covariates_df = covariates_df.loc[:, var_mask_global]

    geno_order = list(genotypes_raw.columns)
    samples_all = [s for s in geno_order if (s in expression_df.columns) and (s in covariates_df.index)]
    if len(samples_all) == 0:
        samples_all = [s for s in geno_order if (s in expression_df.columns)]
        covariates_df = pd.DataFrame(index=samples_all)
    genotypes_raw = genotypes_raw.loc[:, samples_all]
    expression_df = expression_df.loc[:, samples_all]
    covariates_df = covariates_df.reindex(samples_all)

    genotypes_df = genotypes_raw.applymap(convert_gt_cell)
    eligible_sv_mask = genotypes_df.apply(meets_criteria, axis=1)
    varied_genotypes_df = genotypes_df.loc[eligible_sv_mask]
    filtered_genes = expression_df.loc[expression_df.nunique(axis=1) >= 3]

    pairs_df = autodetect_read(pairs_file)
    cols_lower = {c.lower(): c for c in pairs_df.columns}
    gene_col = cols_lower['gene']
    sv_col = cols_lower['sv']

    results = []
    n_skipped_svdn = 0
    n_singular = 0
    candidate_covars = covariates_df.columns.tolist()

    for _, r in pairs_df.iterrows():
        gene = str(r[gene_col])
        sv = str(r[sv_col])
        if (gene not in filtered_genes.index) or (sv not in varied_genotypes_df.index):
            continue
        y = filtered_genes.loc[gene].rename('expr')
        gt = genotypes_df.loc[sv].rename('GT')
        Xc = covariates_df.copy()
        df = pd.concat([y, gt, Xc], axis=1).dropna()
        gt_counts = df['GT'].value_counts().reindex([0.0, 1.0, 2.0]).fillna(0)
        if int(gt_counts.min()) < 37:
            n_skipped_svdn += 1
            continue
        covar_cols_present = [c for c in candidate_covars if c in df.columns]
        covar_cols_pair = [c for c in covar_cols_present if df[c].nunique(dropna=True) > 1]
        X_fit = df[['GT'] + covar_cols_pair]
        y_fit = df['expr'].astype(float)
        X_fit = sm.add_constant(X_fit, has_constant='add')
        if X_fit.shape[0] <= X_fit.shape[1]:
            n_singular += 1
            continue
        try:
            model = sm.OLS(y_fit, X_fit).fit()
            if 'GT' in model.pvalues and not np.isnan(model.pvalues['GT']):
                results.append({
                    'gene': gene,
                    'SV': sv,
                    'p_value': float(model.pvalues['GT']),
                    'coef': float(model.params.get('GT', np.nan))
                })
        except Exception:
            n_singular += 1
            continue

    results_df = pd.DataFrame(results, columns=['gene', 'SV', 'p_value', 'coef'])
    results_df.to_csv(out_pair_results, index=False)

    if results_df.empty:
        pd.DataFrame(columns=['gene', 'SV', 'p_value', 'coef']).to_csv(out_pair_results, index=False)
        pd.DataFrame(columns=['gene', 'SV', 'bonferroni_corrected_pval', 'bh_corrected_pval']).to_csv(out_gene_bh, index=False)
        return

    agg = results_df.groupby('gene').agg(min_p=('p_value', 'min'), n_tests=('p_value', 'size'))
    idx_min = results_df.groupby('gene')['p_value'].idxmin()
    sv_at_min = results_df.loc[idx_min, ['gene', 'SV']].set_index('gene')
    agg = agg.join(sv_at_min, how='left')
    agg['bonferroni_corrected_pval'] = (agg['min_p'] * agg['n_tests']).clip(upper=1.0)

    corrected_pvals_df = agg[['SV', 'bonferroni_corrected_pval']].reset_index()
    _, qvals, _, _ = multipletests(corrected_pvals_df['bonferroni_corrected_pval'].values, alpha=0.05, method='fdr_bh')
    corrected_pvals_df['bh_corrected_pval'] = qvals
    corrected_pvals_df.to_csv(out_gene_bh, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--geno", dest="geno_file", default='genotypes.csv')
    parser.add_argument("--expr", dest="expr_file", default='TMM_expression.csv')
    parser.add_argument("--covar", dest="covar_file", default='covariates.csv')
    parser.add_argument("--pairs", dest="pairs_file", default='SV_gene_1Mb')
    parser.add_argument("--out-pairs", dest="out_pair_results", default='eQTL_result.csv')
    parser.add_argument("--out-bh", dest="out_gene_bh", default='eQTL_result.bh.csv')
    args = parser.parse_args()
    run_eqtl(args.geno_file, args.expr_file, args.covar_file, args.pairs_file, args.out_pair_results, args.out_gene_bh)

if __name__ == "__main__":
    main()
