#!/usr/bin/env python
import argparse,json
from common import load_yaml,iter_smiles_ids,standardize_smiles,basic_descriptors,passes_library_filter
ap=argparse.ArgumentParser()
ap.add_argument("--input",required=True); ap.add_argument("--output",required=True)
ap.add_argument("--config",required=True); ap.add_argument("--qc",required=True)
a=ap.parse_args(); c=load_yaml(a.config)["library"]
st={"read":0,"valid":0,"kept":0,"invalid":0,"filtered":0}; seen=set()
with open(a.output,"w") as fo:
    for smi, mid in iter_smiles_ids(a.input):
        st["read"]+=1; cs,m=standardize_smiles(smi)
        if m is None: st["invalid"]+=1; continue
        st["valid"]+=1; d=basic_descriptors(m)
        if not passes_library_filter(d, c): st["filtered"]+=1; continue
        if (cs,mid) in seen: continue
        seen.add((cs,mid)); st["kept"]+=1; fo.write(f"{cs}\t{mid}\n")
json.dump(st,open(a.qc,"w"),indent=2)
