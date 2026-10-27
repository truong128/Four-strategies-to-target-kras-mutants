#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config(); t=root/'results'/'tables'; f=root/'results'/'figures'
r=pd.read_csv(t/'ranked_hits.csv'); p=pd.read_csv(t/'ligand_properties.csv')
merged=r.merge(p,left_on='ligand',right_on='id',how='left')
if (t/'admet_ai.csv').exists():
    a=pd.read_csv(t/'admet_ai.csv')
    if 'id' in a.columns:
        merged=merged.merge(a.drop(columns=['smiles'],errors='ignore').drop_duplicates('id'),on='id',how='left',suffixes=('','_admet'))
    elif 'smiles' in a.columns:
        merged=merged.merge(a.drop_duplicates('smiles'),on='smiles',how='left',suffixes=('','_admet'))
merged.to_csv(t/'final_candidates.csv',index=False)
# plot selectivity vs resistance
x=merged.WT_selectivity_margin; y=merged.resistance_worst_penalty
plt.figure(figsize=(8,6)); plt.scatter(x,y,s=35)
for _,q in merged.head(15).iterrows():
    if np.isfinite(q.WT_selectivity_margin) and np.isfinite(q.resistance_worst_penalty): plt.text(q.WT_selectivity_margin,q.resistance_worst_penalty,str(q.ligand),fontsize=7)
plt.axvline(cfg['ranking']['min_selectivity_margin_kcal'],ls='--'); plt.axhline(cfg['ranking']['max_resistance_penalty_kcal'],ls='--'); plt.xlabel('WT selectivity margin (kcal/mol; larger is better)'); plt.ylabel('Worst resistance penalty (kcal/mol; smaller is better)'); plt.title('KRAS G12D selectivity vs resistance robustness (Vina heuristic)'); plt.tight_layout(); plt.savefig(f/'selectivity_vs_resistance.png',dpi=200); plt.close()
# markdown summary
cols=['ligand','G12D_mean_affinity','WT_affinity','WT_selectivity_margin','resistance_worst_penalty','composite_rank_score','QED']
cols=[c for c in cols if c in merged.columns]
md=['# KRAS G12D E2E run summary','', '## Important interpretation', '',
    '- Docking scores are ranking heuristics, not measured binding free energies.',
    '- Modeled secondary mutants are stress-test models; experimental structures are preferred when listed in config/pdb_panel.tsv.',
    '- This package covers the noncovalent MRTX1133-like Switch-II-pocket track. It does NOT model cyclophilin-A ternary/covalent RAS(ON) chemistry.',
    '- A lead must still pass biochemical KRAS G12D potency, WT/paralog selectivity, pERK, resistant-cell, organoid/TME, PK and in-vivo validation.', '',
    '## Top computational candidates','',merged[cols].head(50).to_markdown(index=False)]
(root/'results'/'SUMMARY.md').write_text('\n'.join(md))
print('Wrote',root/'results'/'SUMMARY.md')
