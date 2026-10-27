#!/usr/bin/env python3
"""Transfer only exact positional chemistries to KRAS candidates, never assert efficacy."""
import csv,argparse,json,collections,re
from pathlib import Path
def rd(p):
 with open(p,newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def wr(p,rows,cols):
 with open(p,'w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--reference-dir',default='reference/publisher');p.add_argument('--guides',default='results/run/02_guides.csv');p.add_argument('--out',default='results/run/restore2_published');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 ref=rd(Path(a.reference_dir)/'publisher_constructs.csv');rs=rd(Path(a.reference_dir)/'publisher_residues.csv');guides=rd(a.guides)
 by=collections.defaultdict(list)
 for r in rs:by[r['construct_key']].append(r)
 output=[];counts=collections.Counter()
 for g in guides:
  seq=g['guide_rna_5to3'];gi=int(g['guide_index_0based']);L=len(seq)
  for c in ref:
   if c['status']!='EXACT_STANDARD_ALPHABET':counts['SOURCE_NEEDS_MANUAL_REVIEW']+=1;continue
   if int(c['length_reported'])!=L:counts['LENGTH_MISMATCH']+=1;continue
   sym=c['symmetry'];mt=re.fullmatch(r'(\d+)-1-(\d+)',sym)
   if not mt or int(mt[1])!=gi:counts['ORPHAN_POSITION_MISMATCH']+=1;continue
   source=sorted(by[c['construct_key']],key=lambda x:int(x['position_1based']))
   if source[gi]['base']!='C' or seq[gi]!='C':counts['ORPHAN_C_NOT_MATCHED']+=1;continue
   if any(r['sugar'].endswith(('SPECIAL_BASE_I','SPECIAL_BASE_X','SPECIAL_BASE_B','SPECIAL_BASE_Z')) or r['link_to_next_3prime'].startswith('TERMINAL') for r in source):counts['SPECIAL_BASE_STRUCTURE_REQUIRED']+=1;continue
   tid=g['guide_id']+'__'+c['construct_key'].replace(' ','_').replace(':','_')
   for i,r in enumerate(source):
    output.append({'transfer_id':tid,'guide_id':g['guide_id'],'template_source_key':c['construct_key'],'source_target':c['target'],'source_row':c['source_row'],'position_1based':i+1,'kr as_base'.replace(' ',''):seq[i],'sugar':r['sugar'],'link_to_next_3prime':r['link_to_next_3prime'],'status':'PUBLISHED_CHEMISTRY_POSITION_TRANSFER_UNVALIDATED_ON_KRAS'})
   counts['TRANSFERRED_CONSTRUCTS']+=1
 wr(out/'transferred_per_residue.csv',output,['transfer_id','guide_id','template_source_key','source_target','source_row','position_1based','kras_base','sugar','link_to_next_3prime','status'])
 report={'n_guides':len(guides),'n_construct_templates':len(ref),'n_transferred_designs':counts['TRANSFERRED_CONSTRUCTS'],'n_residue_rows':len(output),'rejection_reasons':dict(counts),'caveat':'A source-verified chemistry recipe is not validated ADAR edit efficiency in KRAS; PS_RP/SP stereochemistry unspecified; no MD structures generated.'}
 (out/'transfer_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
