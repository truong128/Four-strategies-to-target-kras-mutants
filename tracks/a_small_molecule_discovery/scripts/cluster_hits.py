#!/usr/bin/env python
import argparse,pandas as pd
from rdkit import Chem,DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from common import load_yaml,read_any_table
ap=argparse.ArgumentParser();ap.add_argument("--input",required=True);ap.add_argument("--config",required=True);ap.add_argument("--output",required=True)
a=ap.parse_args();C=load_yaml(a.config)["clustering"];df=read_any_table(a.input);n=int(C["max_input"])
if len(df)>n:df=df.nlargest(n,"ml_probability") if "ml_probability" in df else df.head(n)
fps=[];ix=[]
for i,s in zip(df.index,df.smiles):
    m=Chem.MolFromSmiles(str(s))
    if m:fps.append(AllChem.GetMorganFingerprintAsBitVect(m,int(C["fingerprint_radius"]),nBits=int(C["fingerprint_nbits"])));ix.append(i)
df=df.loc[ix].reset_index(drop=True)
if not len(df):df.to_csv(a.output,index=False);raise SystemExit
D=[]
for i in range(1,len(fps)):D.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])])
clusters=Butina.ClusterData(D,len(fps),1-float(C["tanimoto_cutoff"]),isDistData=True)
O=[]
for j,cl in enumerate(clusters):
    s=df.iloc[list(cl)].copy()
    if "ml_probability" in s:s=s.sort_values("ml_probability",ascending=False)
    elif "docking_score" in s:s=s.sort_values("docking_score")
    s=s.head(int(C["keep_per_cluster"]));s["cluster_id"]=j;O.append(s)
pd.concat(O,ignore_index=True).to_csv(a.output,index=False)
