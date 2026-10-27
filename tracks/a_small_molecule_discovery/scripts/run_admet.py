#!/usr/bin/env python
import argparse,subprocess,shutil,os,pandas as pd
from common import load_yaml,read_any_table
ap=argparse.ArgumentParser()
for x in ["input","config","predictions","filtered"]:ap.add_argument("--"+x,required=True)
a=ap.parse_args();C=load_yaml(a.config)["admet"];exe=C.get("exe","admet_predict")
if shutil.which(exe) is None:raise SystemExit("ADMET-AI not found; install admet-ai")
df=read_any_table(a.input)
n=int(C.get("keep_top_n",0))
if n and len(df)>n:
    df=df.nlargest(n,"ml_probability") if "ml_probability" in df.columns else df.head(n)
tmp="_admet_input.csv";df[["smiles"]].to_csv(tmp,index=False)
q=subprocess.run([exe,"--data_path",tmp,"--save_path",a.predictions,"--smiles_column","smiles"])
if q.returncode:raise SystemExit(q.returncode)
p=pd.read_csv(a.predictions)
if len(p)==len(df):
    for c in df.columns:
        if c not in p:p[c]=df[c].values
p.to_csv(a.predictions,index=False)
aliases={"max_hERG_probability":["hERG","herg","hERG_probability","hERG_inhibition"],
"max_AMES_probability":["AMES","Ames","ames","AMES_probability","mutagenicity"],
"max_hepatotoxicity_probability":["Hepatotoxicity","hepatotoxicity","DILI","dili","hepatotoxicity_probability"]}
keep=pd.Series(True,index=p.index);applied=[]
for rule,v in C.get("hard_filters",{}).items():
    col=next((x for x in aliases.get(rule,[]) if x in p.columns),None)
    if col:
        z=pd.to_numeric(p[col],errors="coerce");keep&=z.isna()|(z<=float(v));applied.append(f"{col}<={v}")
out=p[keep].copy();out["admet_filters_applied"]=";".join(applied)
n=int(C.get("keep_top_n",0))
if n and len(out)>n:
    out=out.nlargest(n,"ml_probability") if "ml_probability" in out else out.head(n)
out.to_csv(a.filtered,index=False)
try:os.remove(tmp)
except OSError:pass
