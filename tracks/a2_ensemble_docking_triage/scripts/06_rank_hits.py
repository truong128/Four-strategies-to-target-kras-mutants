#!/usr/bin/env python3
import csv
import numpy as np, pandas as pd
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config(); panel=pd.read_csv(root/'config'/'pdb_panel.tsv',sep='\t'); dock=pd.read_csv(root/'results'/'tables'/'docking_long.csv')
d=dock.pivot(index='ligand',columns='pdb_id',values='affinity')
primary=[cfg['ranking']['primary_gdp'],cfg['ranking']['primary_active']]
exp_res=panel.query("role=='resistance' and source_type=='experimental'")['pdb_id'].tolist()
all_res=panel.query("role=='resistance'")['pdb_id'].tolist()
wt=cfg['ranking']['wt']
out=[]
for lig,row in d.iterrows():
    prim=np.nanmean([row.get(x,np.nan) for x in primary]); wt_e=row.get(wt,np.nan)
    exp=np.array([row.get(x,np.nan) for x in exp_res],float); allr=np.array([row.get(x,np.nan) for x in all_res],float)
    exp_mean=np.nanmean(exp) if np.isfinite(exp).any() else np.nan
    all_mean=np.nanmean(allr) if np.isfinite(allr).any() else np.nan
    worst=np.nanmax(allr) if np.isfinite(allr).any() else np.nan  # least favorable (less negative)
    select=(wt_e-prim) if np.isfinite(wt_e) and np.isfinite(prim) else np.nan
    penalty=(all_mean-prim) if np.isfinite(all_mean) and np.isfinite(prim) else np.nan
    worst_pen=(worst-prim) if np.isfinite(worst) and np.isfinite(prim) else np.nan
    out.append({'ligand':lig,'G12D_mean_affinity':prim,'WT_affinity':wt_e,'WT_selectivity_margin':select,
                'resistance_exp_mean':exp_mean,'resistance_all_mean':all_mean,'resistance_mean_penalty':penalty,
                'resistance_worst_penalty':worst_pen})
r=pd.DataFrame(out)
# rank-based composite avoids pretending Vina kcal/mol are quantitative free energies.
def pct_good(s, lower_better=True):
    rank=s.rank(method='average',ascending=lower_better,pct=True)
    return 1-rank+1/len(s)
r['s_affinity']=pct_good(r['G12D_mean_affinity'],True)
r['s_selectivity']=pct_good(-r['WT_selectivity_margin'],True) # larger margin is better
r['s_robustness']=pct_good(r['resistance_worst_penalty'],True)
w=cfg['ranking']['weights']; r['composite_rank_score']=w['g12d_affinity']*r.s_affinity+w['wt_selectivity']*r.s_selectivity+w['resistance_robustness']*r.s_robustness
r['pass_selectivity']=r.WT_selectivity_margin>=cfg['ranking']['min_selectivity_margin_kcal']
r['pass_resistance']=r.resistance_worst_penalty<=cfg['ranking']['max_resistance_penalty_kcal']
r=r.sort_values(['pass_selectivity','pass_resistance','composite_rank_score'],ascending=[False,False,False])
r.to_csv(root/'results'/'tables'/'ranked_hits.csv',index=False)
print(r.head(25).to_string(index=False))
