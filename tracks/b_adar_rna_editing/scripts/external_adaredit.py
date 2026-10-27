#!/usr/bin/env python3
"""Auditable bridge to *externally generated* AdarEdit predictions.
Does not misinterpret Alu-native-site editability as modified-guide editing yield.
Never assumes undocumented checkpoint inference APIs.
"""
import argparse,csv,hashlib,json,re
from pathlib import Path

def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 a=argparse.ArgumentParser();a.add_argument('--guides',required=True);a.add_argument('--predictions',required=True);a.add_argument('--out',required=True);a.add_argument('--model-name',default='AdarEdit');a.add_argument('--model-checkpoint',required=True);a.add_argument('--model-context',required=True);a.add_argument('--prediction-type',required=True,choices=['native_site_editability_probability','native_site_editability_label']);a.add_argument('--model-version',required=True);args=a.parse_args()
 gs=read(args.guides);ps=read(args.predictions)
 if not ps or not gs:raise SystemExit('Missing guide/prediction rows')
 gmap={r['guide_id']:r for r in gs}
 if len(gmap)!=len(gs):raise SystemExit('Duplicate guide IDs')
 required={'guide_id','model_output'}
 if not required.issubset(ps[0]):raise SystemExit('predictions require guide_id,model_output. These must originate from a verified inference execution, not generated placeholders.')
 seen=set();rows=[]
 for p in ps:
  gid=p['guide_id'];x=p['model_output']
  if gid in seen or gid not in gmap:raise SystemExit('Duplicate/unknown guide_id: '+gid)
  seen.add(gid)
  if args.prediction_type.endswith('probability'):
   try:val=float(x)
   except:raise SystemExit('Invalid probability: '+x)
   if not 0<=val<=1:raise SystemExit('Out-of-range probability')
  else:
   if x not in ('0','1'):raise SystemExit('Label must be 0 or 1')
  rows.append({'guide_id':gid,'native_site_model_output':x,'model_context':args.model_context,'model_checkpoint':args.model_checkpoint,'model_version':args.model_version,'status':'OUT_OF_DOMAIN_FOR_CHEMICALLY_MODIFIED_RESTORE2_GUIDES','predicted_RESTORE2_editing_fraction':'','clinical_interpretation':'NOT_VALIDATED'})
 if len(seen)!=len(gs):raise SystemExit(f'Missing predictions for {len(gs)-len(seen)} guides; refusing partial import')
 out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 with (out/'11_adaredit_baseline.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 metadata={'model':args.model_name,'model_context':args.model_context,'model_version':args.model_version,'model_checkpoint':args.model_checkpoint,'prediction_type':args.prediction_type,'external_prediction_file_sha256':sha(args.predictions),'guide_file_sha256':sha(args.guides),'number_guides':len(gs),'warning':'AdarEdit checkpoints trained primarily on endogenous Alu-derived RNA-editing substrates. It is not a validated editing-efficiency model for KRAS G12D RESTORE 2.0 modified guides. Scores must never be substituted for measured editing fraction.'}
 (out/'11_adaredit_provenance.json').write_text(json.dumps(metadata,indent=2))
 print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()
