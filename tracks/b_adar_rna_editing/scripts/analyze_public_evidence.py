#!/usr/bin/env python3
"""Summarize public ADAR/KRAS expression evidence and draw figures (reads output of fetch_public_evidence.py)."""
import argparse, json, re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CANDIDATES = ['ASPC1_PANCREAS', 'PANC1_PANCREAS', 'SW1990_PANCREAS', 'HPAFII_PANCREAS', 'BXPC3_PANCREAS']
IFN = ['ISG15', 'IFIT1', 'MX1', 'OAS1', 'STAT1']


def g12d(s):
    return 'G12D' in str(s).split(';')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    data, out = Path(a.data), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {}

    ccle = pd.read_csv(data / 'ccle_expression_kras.csv')
    det = pd.read_csv(data / 'kras_mutation_details.csv')
    ccle['is_pancreas'] = ccle.lineage.eq('PANCREAS')
    ccle['G12D'] = ccle.KRAS_mutations.map(g12d)
    ccle['ADAR_pct_all_lines'] = ccle.ADAR.rank(pct=True) * 100
    ccle['ADARB1_pct_all_lines'] = ccle.ADARB1.rank(pct=True) * 100
    ccle['IFN_score'] = ccle[IFN].mean(axis=1)
    d = det[det.study == 'ccle_broad_2025'].copy()
    d = d[d.proteinChange.eq('G12D')]
    d['G12D_VAF'] = d.t_alt_count / (d.t_alt_count + d.t_ref_count)
    ccle = ccle.merge(d[['sample', 'G12D_VAF']], on='sample', how='left')

    panc = ccle[ccle.is_pancreas]
    cand = ccle[ccle['sample'].isin(CANDIDATES)].set_index('sample').loc[CANDIDATES].reset_index()
    cand['label'] = cand['sample'].str.replace('_PANCREAS', '', regex=False)
    cols = ['label', 'KRAS_mutations', 'G12D_VAF', 'ADAR', 'ADAR_pct_all_lines', 'ADARB1', 'ADARB1_pct_all_lines', 'ADARB2', 'KRAS', 'IFN_score']
    cand[cols].round(3).to_csv(out / 'candidate_cell_lines.csv', index=False)

    g12d_p = panc[panc.G12D]
    other_p = panc[~panc.G12D]
    res['ccle_n_lines'] = int(len(ccle))
    res['ccle_n_pancreas'] = int(len(panc))
    res['ccle_n_pancreas_G12D'] = int(len(g12d_p))
    res['ADAR_log2tpm1_median_all_lines'] = float(ccle.ADAR.median())
    res['ADAR_log2tpm1_median_pancreas'] = float(panc.ADAR.median())
    res['ADAR_log2tpm1_median_pancreas_G12D'] = float(g12d_p.ADAR.median())
    res['ADAR_log2tpm1_min_pancreas'] = float(panc.ADAR.min())
    res['ADAR_pancreas_lines_below_all_line_median'] = int((panc.ADAR < ccle.ADAR.median()).sum())
    res['ADARB1_log2tpm1_median_pancreas'] = float(panc.ADARB1.median())
    res['ADARB2_log2tpm1_median_pancreas'] = float(panc.ADARB2.median())
    u = stats.mannwhitneyu(g12d_p.ADAR, other_p.ADAR)
    res['ADAR_G12D_vs_other_pancreas_MWU_p'] = float(u.pvalue)
    rho, p = stats.spearmanr(panc.ADAR, panc.IFN_score)
    res['ADAR_vs_IFN_score_spearman_pancreas'] = [float(rho), float(p)]
    res['ADAR_vs_KRAS_expr_spearman_pancreas'] = list(map(float, stats.spearmanr(panc.ADAR, panc.KRAS)))
    res['kras_mutation_counts_pancreas'] = panc.KRAS_mutations.value_counts().to_dict()
    panc.sort_values('ADAR')[['sample', 'KRAS_mutations', 'G12D_VAF', 'ADAR', 'ADAR_pct_all_lines', 'ADARB1', 'ADARB2', 'IFN_score']].round(3).to_csv(out / 'ccle_pancreas_all_lines.csv', index=False)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    order = ccle.groupby('lineage').ADAR.median().sort_values().index
    sel = [l for l in order if (ccle.lineage == l).sum() >= 15]
    ax[0].boxplot([ccle[ccle.lineage == l].ADAR for l in sel], orientation="horizontal", tick_labels=[l.replace('_', ' ')[:28] for l in sel], showfliers=False)
    ax[0].axvline(ccle.ADAR.median(), color='grey', ls=':')
    pos = sel.index('PANCREAS') + 1
    ax[0].get_yticklabels()[pos - 1].set_color('crimson')
    ax[0].set_xlabel('ADAR mRNA, log2(TPM+1), CCLE 2025')
    ax[0].set_title('ADAR1 across lineages (>=15 lines)')
    ax[0].tick_params(axis='y', labelsize=6)
    x = np.arange(len(cand))
    ax[1].bar(x, cand.ADAR, color=['#c0392b' if 'G12D' in str(k) else '#7f8c8d' for k in cand.KRAS_mutations], label='ADAR (p110+p150)')
    ax[1].bar(x, cand.ADARB1, color='#2980b9', alpha=0.6, width=0.4, label='ADARB1 (ADAR2)')
    ax[1].set_xticks(x)
    ax[1].set_xticklabels([f"{l}\n{k if k != 'none_called' else 'KRAS WT'}" for l, k in zip(cand.label, cand.KRAS_mutations)], fontsize=8)
    ax[1].set_ylabel('log2(TPM+1)')
    ax[1].set_title('Candidate validation lines')
    ax[1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out / 'fig1_ccle_adar.png', dpi=200)
    plt.close()

    tcga = pd.read_csv(data / 'tcga_paad_gdc_expression_kras.csv')
    tcga['KRAS_class'] = tcga.KRAS_mutations.map(lambda s: 'none_called' if s == 'none_called' else ('G12D' if g12d(s) else ('G12V' if 'G12V' in s else ('G12R' if 'G12R' in s else 'other_KRAS'))))
    res['tcga_n'] = int(len(tcga))
    res['tcga_KRAS_class_counts'] = tcga.KRAS_class.value_counts().to_dict()
    res['tcga_ADAR_tpm_median'] = float(tcga.ADAR.median())
    res['tcga_ADAR_tpm_IQR'] = [float(tcga.ADAR.quantile(.25)), float(tcga.ADAR.quantile(.75))]
    res['tcga_ADARB1_tpm_median'] = float(tcga.ADARB1.median())
    tg = tcga[tcga.KRAS_class == 'G12D']
    to = tcga[tcga.KRAS_class != 'G12D']
    res['tcga_ADAR_G12D_median'] = float(tg.ADAR.median())
    res['tcga_ADAR_nonG12D_median'] = float(to.ADAR.median())
    res['tcga_ADAR_G12D_vs_other_MWU_p'] = float(stats.mannwhitneyu(tg.ADAR, to.ADAR).pvalue)
    res['tcga_G12D_ADAR_tpm_min'] = float(tg.ADAR.min())
    res['tcga_G12D_ADAR_tpm_q10'] = float(tg.ADAR.quantile(.1))
    res['tcga_ADAR_vs_IFN_spearman'] = list(map(float, stats.spearmanr(tcga.ADAR, tcga[IFN].mean(axis=1))))
    ccle_tpm = np.power(2.0, g12d_p.ADAR) - 1
    res['ccle_G12D_pancreas_ADAR_tpm_median'] = float(ccle_tpm.median())

    cp = pd.read_csv(data / 'cptac_pdac_protein.csv')
    cm = pd.read_csv(data / 'cptac_pdac_mrna_expression_kras.csv')
    cpm = cm.rename(columns={'KRAS_mutations': 'KRAS_mutations_cptac'})
    cell = pd.to_numeric(cp.neoplastic_cellularity, errors='coerce')
    ok = cp.ADAR.notna() & cell.notna()
    res['cptac_n_protein'] = int(cp.ADAR.notna().sum())
    if ok.sum() > 10:
        res['cptac_ADAR_protein_vs_cellularity_spearman'] = list(map(float, stats.spearmanr(cp.ADAR[ok], cell[ok])))
    res['cptac_ADAR_protein_present_fraction'] = float(cp.ADAR.notna().mean())
    m = cp.merge(cm[['sample', 'KRAS_mutations']], on='sample', how='left')
    if m.ADAR.notna().sum() > 10:
        mg = m[m.KRAS_mutations.map(g12d) & m.ADAR.notna()]
        mo = m[~m.KRAS_mutations.map(g12d) & m.ADAR.notna()]
        res['cptac_ADAR_protein_G12D_n'] = int(len(mg))
        res['cptac_ADAR_protein_G12D_vs_other_p'] = float(stats.mannwhitneyu(mg.ADAR, mo.ADAR).pvalue)

    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    ax[0].hist(tcga.ADAR, bins=30, color='#7f8c8d')
    ax[0].axvline(tcga.ADAR.median(), color='crimson')
    ax[0].set_xlabel('ADAR TPM (TCGA-PAAD GDC)')
    ax[0].set_ylabel('tumors')
    ax[0].set_title(f'ADAR1 in {len(tcga)} PDAC tumors')
    groups = ['G12D', 'G12V', 'G12R', 'other_KRAS', 'none_called']
    ax[1].boxplot([tcga[tcga.KRAS_class == k].ADAR for k in groups], tick_labels=[f'{k}\n(n={int((tcga.KRAS_class == k).sum())})' for k in groups], showfliers=False)
    ax[1].set_ylabel('ADAR TPM')
    ax[1].set_title('By KRAS genotype')
    if ok.sum() > 10:
        ax[2].scatter(cell[ok], cp.ADAR[ok], s=12)
        ax[2].set_xlabel('Neoplastic cellularity (%)')
        ax[2].set_ylabel('ADAR1 protein (log2 TMT ratio)')
        ax[2].set_title('CPTAC-PDAC: ADAR1 protein vs purity')
    plt.tight_layout()
    plt.savefig(out / 'fig2_tcga_cptac_adar.png', dpi=200)
    plt.close()

    json.dump(res, open(out / 'public_evidence_summary.json', 'w'), indent=2, default=str)
    print(json.dumps(res, indent=2, default=str))
    print(cand[cols].round(2).to_string())


if __name__ == '__main__':
    main()
