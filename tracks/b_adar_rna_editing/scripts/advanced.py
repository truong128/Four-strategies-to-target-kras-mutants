#!/usr/bin/env python3
"""Steps 7-12. Explicitly distinguishes measured, computed, and unavailable evidence."""
import argparse,csv,json,subprocess,shutil,re,hashlib,math,random
from pathlib import Path
from collections import defaultdict

def readcsv(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def output(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def fasta(path):
 name=None;seq=[]
 with open(path) as f:
  for line in f:
   line=line.strip()
   if line.startswith('>'):
    if name is not None:yield name,''.join(seq).upper().replace('U','T');seq=[]
    name=line[1:].split()[0]
   elif line:seq.append(line)
  if name is not None:yield name,''.join(seq).upper().replace('U','T')
def revcomp(seq):return seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
def run(args):
 d=Path(args.input);out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 guides=readcsv(d/'02_guides.csv'); bystanders=readcsv(d/'03_bystander_sites.csv')
 if not guides:raise ValueError('No guides generated. Check input variant and reference transcript.')
 # Step 7: per-position chemical recipe; these are design HYPOTHESES, not published RESTORE2 recipes.
 chem=[]; sites=[]
 for g in guides:
  L=len(g['guide_rna_5to3']);center=int(g['guide_index_0based'])
  for recipe in ('UNMODIFIED_PO','TERMINAL_2OME_PS_EXPLORATORY'):
   cid=g['guide_id']+'__'+recipe
   chem.append(dict(construct_id=cid,guide_id=g['guide_id'],recipe=recipe,verified_published_pattern='NO',parameterized_for_MD='NO',status='DESIGN_HYPOTHESIS' if recipe!='UNMODIFIED_PO' else 'CONTROL'))
   for i,base in enumerate(g['guide_rna_5to3']):
    terminal= i<3 or i>=L-3
    sugar='2OME' if recipe!='UNMODIFIED_PO' and terminal else 'RNA_2OH'
    # phosphorothioate link defined on 3-prime outgoing bond, terminal PS only
    linkage='NONE_3PRIME' if i==L-1 else ('PS_STEREO_UNSPECIFIED' if recipe!='UNMODIFIED_PO' and (i<3 or i>=L-4) else 'PO')
    sites.append(dict(construct_id=cid,guide_position_1based=i+1,base=base,sugar=sugar,link_to_next_3prime=linkage,is_opposite_target_A=str(i==center),MD_ready='NO' if recipe!='UNMODIFIED_PO' else 'NEEDS_3D_INPUT'))
 output(out/'07_chemistry_constructs.csv',chem,list(chem[0]))
 output(out/'07_chemistry_positions.csv',sites,list(sites[0]))
 # Step 8: correct RNA/RNA duplex topology via RNAcofold, not RNAfold single strand
 folds=[];exe=shutil.which('RNAcofold')
 for g in guides:
  entry=dict(guide_id=g['guide_id'],status='UNAVAILABLE_RNAcofold',mfe_kcal_mol='',structure='',warning='Unmodified-RNA approximation; modified chemical pattern not modeled')
  if exe:
   duplex=g['target_rna_5to3']+'&'+g['guide_rna_5to3']
   p=subprocess.run([exe,'--noPS'],input=duplex+'\n',text=True,capture_output=True,timeout=45)
   if p.returncode==0:
    lines=p.stdout.strip().splitlines()
    pattern=r'^(\S+\&\S+)\s+\(\s*([-+\d.]+)\)'
    match=next((re.search(pattern,line) for line in lines if re.search(pattern,line)),None)
    if match:entry.update(status='OK_UNMODIFIED_ONLY',structure=match.group(1),mfe_kcal_mol=match.group(2))
    else:entry['status']='PARSE_FAILED'
   else:entry['status']='FAILED_'+p.stderr[:70]
  folds.append(entry)
 output(out/'08_cofold.csv',folds,list(folds[0]))
 # Step 9: transcriptome-wide complementarity enumeration (NOT biological off-target clearance).
 hits=[]; offstatus='NOT_RUN_NO_TRANSCRIPTOME'
 if args.transcriptome:
  if not Path(args.transcriptome).exists():raise ValueError('transcriptome FASTA does not exist')
  offstatus='FULL_LENGTH_SEQUENCE_SCAN_COMPLETED'
  sequences=list(fasta(args.transcriptome))
  for g in guides:
   # exact full-length guide complement, with up to N mismatches; strand orientation explicit
   target=revcomp(g['guide_rna_5to3'].replace('U','T'))
   for name,seq in sequences:
    for strand,t in [('+',seq)]:
     if len(t)<len(target):continue
     # Pigeonhole indexing: for <=k mismatches at least one of k+1 disjoint
     # segments must be exactly matched. For short sequences, fall back to scan.
     k=args.max_mismatches
     if k<0:raise ValueError('max-mismatches must be >=0')
     candidates=set()
     if len(target)>=8*(k+1):
      segments=k+1
      for seg in range(segments):
       lo=seg*len(target)//segments;hi=(seg+1)*len(target)//segments
       seed=target[lo:hi]
       where=t.find(seed)
       while where!=-1:
        start=where-lo
        if 0<=start<=len(t)-len(target):candidates.add(start)
        where=t.find(seed,where+1)
     else:candidates=range(len(t)-len(target)+1)
     for i in sorted(candidates):
      mismatches=sum(a!=b for a,b in zip(target,t[i:i+len(target)]))
      if mismatches<=k:
       hits.append(dict(guide_id=g['guide_id'],transcript_id=name,strand=strand,start_1based=i+1,mismatches=mismatches,hit_sequence=t[i:i+len(target)],status='SEQUENCE_COMPLEMENTARITY_ONLY'))
   if len(hits)>args.max_hits:raise RuntimeError('Too many off-target hits; restrict inputs or increase --max-hits; partial output intentionally suppressed')
 output(out/'09_transcriptome_hits.csv',hits,['guide_id','transcript_id','strand','start_1based','mismatches','hit_sequence','status'])
 (out/'09_offtarget_limits.json').write_text(json.dumps(dict(status=offstatus,maximum_mismatches=args.max_mismatches,transcriptome_fasta=args.transcriptome,limitations=['Scanning full-length near-complementarity, not bulges/partial binding','Transcript identifiers may be redundant','No calibrated ADAR editing prediction','Sequence matches do not demonstrate transcriptome-wide safety','Only sequences present in supplied FASTA evaluated','Sense transcript sequences only; genomic reverse-strand loci must be represented as sense transcripts in FASTA']),indent=2))
 # Step 10: MD execution gate. Accept verified input manifest; never generate fake modified forcefields.
 md=[]; manifest={}
 if args.md_manifest:
  manifest={r['construct_id']:r for r in readcsv(args.md_manifest)}
 for c in chem:
  m=manifest.get(c['construct_id'],{})
  file_names=('structure_pdb','topology_prmtop','coordinates_inpcrd')
  present=all(m.get(k) and Path(m[k]).is_file() for k in file_names)
  metadata=all(m.get(k) for k in ('chemistry_verified','stereochemistry_verified','parameterization_validated','template_validated'))
  qualified=present and metadata and all(m[k].strip().upper()=='YES' for k in ('chemistry_verified','stereochemistry_verified','parameterization_validated','template_validated'))
  md.append(dict(construct_id=c['construct_id'],status='READY_FOR_EXTERNAL_MD' if qualified else 'BLOCKED_MISSING_VALIDATED_INPUTS',structure_pdb=m.get('structure_pdb',''),topology_prmtop=m.get('topology_prmtop',''),coordinates_inpcrd=m.get('coordinates_inpcrd',''),evidence_of_run='NOT_RUN'))
 output(out/'10_md_readiness.csv',md,list(md[0]))
 # Step 11: cross-dataset experimental model: no measurements -> strictly abstain.
 preds=[];mlstatus='ABSTAIN_NO_LABELED_EXPERIMENTAL_DATA'
 if args.measured_editing:
  try:
   import numpy as np
   from sklearn.ensemble import RandomForestRegressor
   from sklearn.model_selection import GroupKFold, cross_val_predict
   from sklearn.metrics import mean_absolute_error
  except ImportError:raise RuntimeError('ML option requires numpy and scikit-learn: pip install -r requirements-ml.txt')
  training=readcsv(args.measured_editing)
  req=['guide_rna_5to3','target_rna_5to3','editing_fraction','study_or_target_group']
  if any(not all(k in r for k in req) for r in training):raise ValueError('Measured editing CSV required columns: '+','.join(req))
  def features(g):
   s=g['guide_rna_5to3'].replace('T','U');t=g['target_rna_5to3'].replace('T','U')
   return [len(s),s.count('G')/len(s),s.count('C')/len(s),s.count('U')/len(s),t.count('A')/len(t),t.count('U')/len(t)]
  X=np.array([features(r) for r in training]);y=np.array([float(r['editing_fraction']) for r in training]);groups=[r['study_or_target_group'] for r in training]
  if len(set(groups))<3 or len(y)<30 or np.any((y<0)|(y>1)):
   mlstatus='ABSTAIN_INSUFFICIENT_INDEPENDENT_GROUPS_OR_BAD_LABELS'
  else:
   model=RandomForestRegressor(n_estimators=300,min_samples_leaf=3,random_state=42,n_jobs=1)
   cv=GroupKFold(n_splits=min(5,len(set(groups))))
   cross=cross_val_predict(model,X,y,cv=cv,groups=groups)
   mae=float(mean_absolute_error(y,cross));baseline=float(mean_absolute_error(y,np.full_like(y,np.mean(y))))
   if mae >= baseline:
    mlstatus=f'ABSTAIN_MODEL_NOT_BETTER_THAN_BASELINE:MAE={mae:.3f}:baseline={baseline:.3f}'
   else:
    mlstatus=f'EXPLORATORY_MODEL_NOT_VALIDATED_FOR_KRAS:groupCV_MAE={mae:.3f}'
    model.fit(X,y)
    for g in guides:preds.append(dict(guide_id=g['guide_id'],predicted_editing_fraction=round(float(model.predict([features(g)])[0]),4),status=mlstatus))
 output(out/'11_editing_ml.csv',preds,['guide_id','predicted_editing_fraction','status'])
 (out/'11_model_status.json').write_text(json.dumps(dict(status=mlstatus,requirement='Independent KRAS external test set required before clinically meaningful claims'),indent=2))
 # Step 12: no therapeutic efficacy ranking without measured evidence; generate candidate triage only.
 bcounts=defaultdict(int)
 for b in bystanders:
  if b['is_target']=='NO':bcounts[b['guide_id']]+=1
 hits_count=defaultdict(int)
 for h in hits:hits_count[h['guide_id']]+=1
 final=[]
 for g in guides:
  ready=next(x['status'] for x in md if x['construct_id']==g['guide_id']+'__UNMODIFIED_PO')
  final.append(dict(guide_id=g['guide_id'],variant_id=g['variant_id'],potential_bystander_A_count=bcounts[g['guide_id']],transcriptome_sequence_hits=hits_count[g['guide_id']] if args.transcriptome else 'NOT_ASSESSED',WT_allele_specificity='UNPROVEN',editing_efficiency='UNPROVEN',MD_status=ready,therapeutic_priority='NOT_RANKABLE_WITHOUT_EXPERIMENTAL_EVIDENCE',interpretation='Engineering triage only; not a therapeutic ranking'))
 output(out/'12_candidate_evidence.csv',final,list(final[0]))
 summary={'guides':len(guides),'chemistry_constructs':len(chem),'RNAcofold_available':bool(exe),'offtarget_status':offstatus,'offtarget_hits':len(hits),'MD_ready':sum(r['status']=='READY_FOR_EXTERNAL_MD' for r in md),'ML_status':mlstatus,'therapeutic_ranking':'NOT_PERFORMED_UNVALIDATED','warning':'This pipeline neither establishes allele-specific ADAR editing nor PDAC treatment efficacy'}
 (out/'STEPS_07_12_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',default='results/run');p.add_argument('--out',default='results/run');p.add_argument('--transcriptome');p.add_argument('--max-mismatches',type=int,default=2);p.add_argument('--max-hits',type=int,default=200000);p.add_argument('--md-manifest');p.add_argument('--measured-editing');run(p.parse_args())
