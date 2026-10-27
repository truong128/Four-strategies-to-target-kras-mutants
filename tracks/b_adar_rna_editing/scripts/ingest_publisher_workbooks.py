#!/usr/bin/env python3
"""RESTORE 2.0 direct publisher-workbook ingestion; conservative, provenance-first.
Requires artifact_tool; re-import raw workbook, never guess ambiguous conjugate chemistry.
"""
import argparse,csv,hashlib,json,re
from pathlib import Path
from artifact_tool import Blob,SpreadsheetFile
P=re.compile(r'(m|f|d|o|l)?([ACGUTIXBZ])([*]?)',re.I)
MAP={'m':'2OME','f':'2F','d':'DNA_2H','o':'2MOE','l':'LNA','':'RNA_2OH'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,rows,headers):
 with Path(p).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=headers,extrasaction='ignore');w.writeheader();w.writerows(rows)
def main():
 a=argparse.ArgumentParser();a.add_argument('--supplement',required=True);a.add_argument('--source-data',required=True);a.add_argument('--out',default='results/publisher_ingest');args=a.parse_args()
 out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 sup=SpreadsheetFile.import_xlsx(Blob.load(args.supplement)); ss=sup.worksheets.get_item_at(0)
 raw=ss.get_range('A1:E300').values
 constructs=[];residues=[];issues=[];byid={};bykey={}
 for excelrow,r in enumerate(raw,1):
  if not isinstance(r[1],str) or not isinstance(r[2],(int,float)) or not isinstance(r[4],str):continue
  target,ver,L,sym,annot=r[:5];target=str(target or '').strip();ver=ver.strip();s=re.sub(r'\s+','',annot)
  chunks=list(P.finditer(s));extra=P.sub('',s)
  base=''.join(m.group(2).upper().replace('T','U') for m in chunks)
  flags=[]
  if len(chunks)!=int(L):flags.append('PARSED_LENGTH_MISMATCH')
  if any(m.group(2).upper().replace('T','U') not in 'ACGU' for m in chunks):flags.append('SPECIAL_NUCLEOBASE_PARAMETERIZATION_REQUIRED')
  if extra:flags.append('SPECIAL_ANNOTATION_OR_UNPARSED:'+extra[:80])
  if not re.fullmatch(r'\d+-1-\d+',str(sym)):flags.append('UNEXPECTED_SYMMETRY')
  else:
   aa,bb=map(int,sym.split('-1-'))
   if aa+1+bb!=int(L):flags.append('SYMMETRY_LENGTH_MISMATCH')
  cid=f'{target}::{ver}::row{excelrow}'
  status='EXACT_STANDARD_ALPHABET' if not flags else 'MANUAL_REVIEW_REQUIRED'
  constructs.append({'construct_key':cid,'target':target,'on_version':ver,'length_reported':L,'length_parsed':len(chunks),'symmetry':sym,'sequence_5to3':base,'publisher_sequence_raw':annot,'source_row':excelrow,'status':status,'issues':';'.join(flags)})
  byid.setdefault('v'+ver.lstrip('v'),set()).add(target)
  bykey[(target,ver)]=cid
  if flags:issues.append({'construct_key':cid,'issues':';'.join(flags)})
  # Even manual review records retain parsed residues for diagnostic work; NOT MD eligible.
  for pos,m in enumerate(chunks,1):
   pfx,nuc,star=m.group(1) or '',m.group(2).upper().replace('T','U'),m.group(3)
   chem=MAP[pfx.lower()]
   if nuc not in 'ACGU':chem+=':SPECIAL_BASE_'+nuc
   link='PS_RANDOM' if star else 'PO'
   if pos==len(chunks):link='NONE' if not star else 'TERMINAL_PS_MARKER_NEEDS_REVIEW'
   residues.append({'construct_key':cid,'target':target,'on_version':ver,'position_1based':pos,'base':nuc,'sugar':chem,'link_to_next_3prime':link,'status':status,'source_row':excelrow})
 write(out/'publisher_constructs.csv',constructs,list(constructs[0]))
 write(out/'publisher_residues.csv',residues,list(residues[0]))
 write(out/'publisher_parser_issues.csv',issues,['construct_key','issues'])
 src=SpreadsheetFile.import_xlsx(Blob.load(args.source_data))
 assay=[];sheets=[]
 # Source workbook is figure-shaped; use only rows with header "ON name" and numeric replicates, preserve context.
 for ii in range(len(src.worksheets.items)):
  sheet=src.worksheets.get_item_at(ii); vals=sheet.get_range('A1:Q200').values
  sheets.append(sheet.name)
  context=''
  for i,row in enumerate(vals):
   if isinstance(row[0],str) and ('Figure' in row[0]):context=' | '.join(str(x) for x in row[:3] if x is not None)
   for k,x in enumerate(row):
    if not isinstance(x,str) or 'ON name' not in x:continue
    for j in range(i+1,min(i+90,len(vals))):
     r=vals[j];name=r[k]
     if not isinstance(name,str) or not re.match(r'^v?\d+[.]\d+$',name.strip()):
      if j>i+2 and all(x is None for x in r):break
      continue
     nums=[z for z in r[k+1:k+7] if isinstance(z,(float,int)) and not isinstance(z,bool)]
     if not nums:continue
     v='v'+name.strip().lstrip('v');targets=sorted(byid.get(v,[]))
     assay.append({'source_sheet':sheet.name,'context':context,'excel_row':j+1,'on_version':v,'targets_with_same_version':';'.join(targets),'n_possible_targets':len(targets),'numeric_values_raw_json':json.dumps(nums),'status':'AMBIGUOUS_TARGET_DO_NOT_JOIN' if len(targets)!=1 else 'MATCH_CANDIDATE_REQUIRES_ASSAY_CONTEXT'})
 write(out/'source_data_assay_index.csv',assay,['source_sheet','context','excel_row','on_version','targets_with_same_version','n_possible_targets','numeric_values_raw_json','status'])
 report={'article_doi':'10.1038/s41467-025-64434-7','source_sha256':{'Supplementary Data 1':sha(args.supplement),'Source Data':sha(args.source_data)},'n_constructs':len(constructs),'n_exact_standard_records':sum(x['status']=='EXACT_STANDARD_ALPHABET' for x in constructs),'n_need_manual_review':len(issues),'n_residues':len(residues),'source_data_sheets':sheets,'n_source_index_rows':len(assay),'n_ambiguous_source_rows':sum(x['n_possible_targets']!=1 for x in assay),'strong_warning':'Source figure tables do not provide a universal assay-independent editing rate; NEVER train ML using uncurated numeric arrays; chemical patterns transferred to KRAS are UNVALIDATED. PDB/Amber parameterization NOT produced.'}
 (out/'publisher_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
