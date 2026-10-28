#!/usr/bin/env python
import argparse,pandas as pd,numpy as np
from common import load_yaml,read_any_table
ap=argparse.ArgumentParser()
for x in ["base","admet","synthesis","config","ranked","pareto"]:ap.add_argument("--"+x,required=True)
for x in ["boltz","md","mmpbsa","openfe","resistance"]:ap.add_argument("--"+x)
a=ap.parse_args();C=load_yaml(a.config)["ranking"];b=read_any_table(a.base);ad=read_any_table(a.admet);sy=read_any_table(a.synthesis)
keys=[k for k in ["id","smiles"] if k in b and k in ad];df=b.merge(ad,on=keys,how="left",suffixes=("","_admet"))
keys=[k for k in ["id","smiles"] if k in df and k in sy];df=df.merge(sy,on=keys,how="left",suffixes=("","_synth"))
for path in [a.boltz, a.md, a.mmpbsa, a.openfe, a.resistance]:
    if not path: continue
    ex=read_any_table(path); keys=[k for k in ["id","smiles"] if k in df and k in ex]
    if keys: df=df.merge(ex,on=keys,how="left",suffixes=("","_extra"))
ac=[c for c in ad.columns if c not in {"id","smiles","ml_probability","docking_score","dock_status","admet_filters_applied"}]
df["admet_pass_fraction"]=df[ac].notna().mean(axis=1) if ac else np.nan
dirs={"ml_probability":"max","docking_score":"min","gnina_cnn_score":"max","admet_pass_fraction":"max","G12D_mean_affinity":"min","WT_selectivity_margin":"max","resistance_worst_penalty":"min","boltz_affinity_probability":"max","boltz_affinity_pred_value":"min","boltz_confidence_score":"max","md_stability":"max","mmgbsa":"min","openfe_abfe_dg_kcal":"min","synthetic_accessibility":"max","retrosynthesis_score":"max"}
score=np.zeros(len(df));tot=0
for c,w in C["weights"].items():
    if c not in df:continue
    v=pd.to_numeric(df[c],errors="coerce")
    if v.notna().sum()<2:continue
    r=v.rank(pct=True,ascending=(dirs.get(c)=="min"),na_option="bottom");score+=float(w)*r.fillna(1).values;tot+=float(w)
df["weighted_rank_score"]=score/max(tot,1e-9);df=df.sort_values("weighted_rank_score")
objs=[]
for c,d in dirs.items():
    if c in df and pd.to_numeric(df[c],errors="coerce").notna().sum()>1:
        x=pd.to_numeric(df[c],errors="coerce").values.astype(float);x=-x if d=="max" else x
        fin=np.isfinite(x);worst=np.nanmax(x[fin])+abs(np.nanmax(x[fin]))+1;x[~fin]=worst;objs.append(x)
pareto=np.ones(len(df),dtype=bool)
if objs:
    A=np.vstack(objs).T
    for i in range(len(A)):
        if np.any(np.all(A<=A[i],axis=1)&np.any(A<A[i],axis=1)):pareto[i]=False
df["pareto"]=pareto;df.to_csv(a.ranked,index=False);df[df.pareto].to_csv(a.pareto,index=False)
