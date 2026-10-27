#!/usr/bin/env python3
"""Auditable RESTORE 2.0 reference ingestion: never infer chemistry from prose.
Accept verified, explicitly normalized rows exported from publisher supplementary tables.
Pure standard library; untouched source files are SHA256 logged.
"""
import argparse,csv,hashlib,json,re,sys,datetime
from pathlib import Path

VALID_SUGARS={'RNA_2OH','2OME','2F','DNA_2H','2MOE','LNA','OTHER'}
VALID_BONDS={'PO','PS_RANDOM','PS_RP','PS_SP','OTHER','NONE'}
REQUIRED=['construct_id','target_id','guide_5to3','guide_position_1based','base','sugar','link_to_next_3prime','source_table','source_row','evidence_level']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load_csv(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def save(p,rows,cols):
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def validate(rows):
 if not rows:raise ValueError('No reference records supplied')
 if any(k not in rows[0] for k in REQUIRED):raise ValueError('Missing normalized fields: '+','.join(k for k in REQUIRED if k not in rows[0]))
 groups={}
 for ix,r in enumerate(rows,2):
  cid=r['construct_id'];seq=r['guide_5to3'].upper().replace('T','U');p=int(r['guide_position_1based']);base=r['base'].upper().replace('T','U')
  if not re.fullmatch('[ACGU]+',seq):raise ValueError(f'row {ix}: invalid RNA sequence')
  if not 1<=p<=len(seq) or seq[p-1]!=base:raise ValueError(f'row {ix}: base-position mismatch')
  if r['sugar'] not in VALID_SUGARS or r['link_to_next_3prime'] not in VALID_BONDS:raise ValueError(f'row {ix}: unknown chemistry, use OTHER with supporting record')
  if p==len(seq) and r['link_to_next_3prime']!='NONE':raise ValueError(f'row {ix}: terminal outgoing bond must be NONE')
  if p<len(seq) and r['link_to_next_3prime']=='NONE':raise ValueError(f'row {ix}: non-terminal bond missing')
  if not r['source_table'] or not r['source_row']:raise ValueError(f'row {ix}: source provenance required')
  if r['evidence_level'] not in ('PUBLISHED_EXACT','AUTHOR_CONFIRMED','EXPLORATORY') :raise ValueError(f'row {ix}: invalid evidence level')
  k=groups.setdefault(cid,{'seq':seq,'p':set(),'target':r['target_id'],'level':set()})
  if k['seq']!=seq or k['target']!=r['target_id'] or p in k['p']:raise ValueError(f'row {ix}: duplicate/conflicting positions')
  k['p'].add(p);k['level'].add(r['evidence_level'])
 for cid,g in groups.items():
  if g['p']!=set(range(1,len(g['seq'])+1)):raise ValueError(f'{cid}: positions are incomplete')
  if len(g['level'])!=1:raise ValueError(f'{cid}: mixed evidence levels')
 return groups

def guide_mapping(rows,guides):
 # Transfer source *chemical positions* to a novel guide, never transfer original sequence.
 library={}
 for r in rows:library.setdefault(r['construct_id'],[]).append(r)
 out=[];reject=[]
 for g in guides:
  seq=g['guide_rna_5to3'];gi=int(g['guide_index_0based']);L=len(seq)
  for cid,source in library.items():
   source=sorted(source,key=lambda x:int(x['guide_position_1based']))
   # Editing position is required explicitly; it cannot be inferred from an orphan C alone.
   targetpositions=[i for i,s in enumerate(source) if str(s.get('is_orphan_C','')).strip().upper()=='YES']
   if len(source)!=L or len(targetpositions)!=1 or targetpositions[0]!=gi:
    reject.append({'guide_id':g['guide_id'],'template_id':cid,'reason':'LEN_OR_EDIT_SITE_INCOMPATIBLE_OR_UNANNOTATED'});continue
   targetid=g['guide_id']+'__transfer__'+cid
   for i,r in enumerate(source):out.append({'construct_id':targetid,'guide_id':g['guide_id'],'template_id':cid,'position_1based':i+1,'base':seq[i],'sugar':r['sugar'],'link_to_next_3prime':r['link_to_next_3prime'],'status':'TEMPLATE_TRANSFER_NOT_EXPERIMENTALLY_VALIDATED','source_table':r['source_table'],'source_row':r['source_row']})
 return out,reject

def main():
 ap=argparse.ArgumentParser(description='Validate published chemistry reference; exact data must be manually normalized from publisher supplement')
 ap.add_argument('--reference-csv',required=True);ap.add_argument('--guides-csv');ap.add_argument('--out',default='results/restore2');ap.add_argument('--publisher-supplement');ap.add_argument('--publisher-source-data')
 a=ap.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 rows=load_csv(a.reference_csv);groups=validate(rows)
 report={'paper_doi':'10.1038/s41467-025-64434-7','reference_sha256':sha(a.reference_csv),'constructs':len(groups),'records':len(rows),'verified_experimentally_on_KRAS':False,'source_files':{},'notes':['PUBLISHED_EXACT is provenance status of source recipe only; transferred constructs are UNVALIDATED on KRAS','No editing efficiency is derived from chemistry alone','Unknown PS stereo-random ensembles require explicit molecular modeling','Source Data numerical labels must be matched by construct ID and assay conditions']}
 for field in ('publisher_supplement','publisher_source_data'):
  fn=getattr(a,field)
  if fn:
   if not Path(fn).exists():raise FileNotFoundError(fn)
   report['source_files'][field]={'path':str(fn),'sha256':sha(fn)}
 if a.guides_csv:
  trans,rej=guide_mapping(rows,load_csv(a.guides_csv))
  save(out/'restore2_transferred_chemistry.csv',trans,['construct_id','guide_id','template_id','position_1based','base','sugar','link_to_next_3prime','status','source_table','source_row'])
  save(out/'restore2_rejections.csv',rej,['guide_id','template_id','reason'])
  report['transferred_rows']=len(trans);report['rejected_pairs']=len(rej)
 (out/'reference_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
