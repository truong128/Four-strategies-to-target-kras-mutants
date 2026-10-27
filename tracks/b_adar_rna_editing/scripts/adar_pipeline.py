#!/usr/bin/env python3
"""KRAS ADAR E2E: reproducible computational prioritization, not an editing efficacy predictor."""
import argparse,csv,json,math,os,re,subprocess,sys,time,urllib.request
from pathlib import Path
from collections import defaultdict
COMP=str.maketrans('ACGTU','TGCAA')
CODONS={}
for aas,letters in [('FF','TTTTTC'),('LL','TTATTGCTTCTCCTACTG'),('II','ATTATCAT A'.replace(' ','')),('M','ATG')]: pass
GENETIC={'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def writecsv(path,rows,cols):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def fetch_cds(tx,cache):
 cache.parent.mkdir(parents=True,exist_ok=True)
 if cache.exists():return ''.join(x.strip() for x in cache.read_text().splitlines() if not x.startswith('>')).upper().replace('U','T')
 url=f'https://rest.ensembl.org/sequence/id/{tx}?type=cds'
 req=urllib.request.Request(url,headers={'Content-Type':'text/plain','User-Agent':'KRAS-ADAR-academic/1.0'})
 try:
  with urllib.request.urlopen(req,timeout=30) as response:s=response.read().decode().strip().upper()
 except Exception as e:raise RuntimeError(f'Ensembl fetch failed: {e}; use --cds-fasta with verified canonical KRAS CDS')
 if not re.fullmatch('[ACGT]+',s):raise ValueError('CDS contains unexpected bases')
 cache.write_text(f'>{tx} CDS Ensembl REST\n'+s+'\n');return s
def loadfasta(path):
 return ''.join(line.strip() for line in Path(path).read_text().splitlines() if line and not line.startswith('>')).upper().replace('U','T')
def translate(cds):return ''.join(GENETIC.get(cds[i:i+3],'X') for i in range(0,len(cds)-2,3))
def variants(file):
 with open(file,newline='') as f:return list(csv.DictReader(f))
def edit_variant(wt,changes):
 s=list(wt)
 for c in changes:
  m=re.fullmatch(r'c\.(\d+)([ACGT])>([ACGT])',c)
  if not m:raise ValueError('Only SNV HGVS c.NRef>Alt supported: '+c)
  pos,ref,alt=int(m[1]),m[2],m[3]
  if not 1<=pos<=len(s) or wt[pos-1]!=ref:raise ValueError(f'Reference mismatch {c}; observed {wt[pos-1] if pos<=len(s) else "OUT_OF_RANGE"}')
  s[pos-1]=alt
 return ''.join(s)
def gc(s):return round(100*(s.count('G')+s.count('C'))/len(s),2)
def fold(seq):
 exe=__import__('shutil').which('RNAfold')
 if not exe:return ('NOT_RUN','')
 p=subprocess.run([exe,'--noPS'],input=seq.replace('T','U')+'\n',text=True,capture_output=True,timeout=30)
 if p.returncode:return ('ERROR',p.stderr[:200])
 lines=p.stdout.strip().splitlines()
 return ('OK',lines[-1] if lines else '')
def run(args):
 out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
 cds=loadfasta(args.cds_fasta) if args.cds_fasta else fetch_cds(args.transcript,out/'cache'/(args.transcript+'.fasta'))
 if len(cds)%3 or cds[:3]!='ATG' or cds[33:36] not in ('GGT','GGC','GGA','GGG'):
  raise ValueError(f'Failed CDS sanity checks; codon12={cds[33:36]}, length={len(cds)}; check transcript/isoform')
 panel=variants(args.variants)
 ann=[];guides=[];bystanders=[];chem=[];conflicts=[]
 for row in panel:
  name=row['variant_id'].strip(); alterations=[x.strip() for x in row['hgvs_c'].split(';') if x.strip()]
  mutant=edit_variant(cds,alterations)
  proteinWT=translate(cds);proteinM=translate(mutant)
  positions=[int(re.match(r'c\.(\d+)',x)[1])-1 for x in alterations]
  editable=[]
  for pos in positions:
   if mutant[pos]!='A':continue
   corrected=mutant[:pos]+'G'+mutant[pos+1:]
   idx=pos//3
   if GENETIC[corrected[idx*3:idx*3+3]]==GENETIC[cds[idx*3:idx*3+3]]:
    editable.append(pos)
  ann.append({'variant_id':name,'hgvs_c':row['hgvs_c'],'protein_changes':','.join(f'{proteinWT[i]}{i+1}{proteinM[i]}' for i in range(len(proteinWT)) if proteinWT[i]!=proteinM[i]),'edited_cds_positions_1based':';'.join(str(i+1) for i in editable),'directly_correctable':'YES' if editable else 'NO','other_uncorrected_changes_after_single_edit':'VERIFY' if len(alterations)>1 else 'NONE','transcript':args.transcript})
  for pos in editable:
   for L in args.lengths:
    # Include both generic offsets and published RESTORE2 asymmetric orphan-C locations.
    # guide index is 5'-side of orphan C; target offset is reverse-indexed.
    guide_indices={L-1-x for x in (L//2, L//2-5, L//2+5)}
    guide_indices.update(idx for idx in (29,24,25,26) if 0<=idx<L)
    for offset in sorted({L-1-idx for idx in guide_indices}):
     start=pos-offset
     if start<0 or start+L>len(cds):continue
     fragment=mutant[start:start+L];wfrag=cds[start:start+L]
     guide=list(rc(fragment));gi=L-1-offset
     if guide[gi]!='T':raise ValueError('Unexpected guide base at target')
     guide[gi]='C';g=''.join(guide)
     gid=f'{name}_c{pos+1}_L{L}_O{offset}'
     candidate={'guide_id':gid,'variant_id':name,'target_cds_position':pos+1,'guide_length':L,'target_cds_start_1based':start+1,'target_rna_5to3':fragment.replace('T','U'),'guide_rna_5to3':g.replace('T','U'),'target_index_0based':offset,'guide_index_0based':gi,'gc_percent':gc(g),'wt_local_mismatches':sum(a!=b for a,b in zip(fragment,wfrag)),'wt_guide_target_base':wfrag[offset],'mutant_guide_target_base':fragment[offset],'sequence_basis':'Ensembl CDS or user FASTA; 5-prime/3-prime antisense','wt_site_pairing':'G-C canonical vs mutant A-C mismatch; WT guide binding possible' if wfrag[offset]=='G' and fragment[offset]=='A' else 'INDEPENDENT_VALIDATION_REQUIRED'}
     guides.append(candidate)
     for i,b in enumerate(fragment):
      if b=='A':bystanders.append({'guide_id':gid,'cds_position':start+i+1,'relative_to_target':i-offset,'is_target':'YES' if i==offset else 'NO','context_5prime':fragment[i-1] if i else 'EDGE','context_3prime':fragment[i+1] if i<L-1 else 'EDGE','guide_counterbase':g[L-1-i]})
     for label in ('UNMODIFIED_REFERENCE','EXPLORATORY_TERMINAL_2OMe_2F_PS'):
      if label=='UNMODIFIED_REFERENCE':notes='All ribose OH, all PO; control only'
      else:notes='Illustrative terminal modifications, NOT exact RESTORE 2.0 pattern; no efficacy claim; central window unmodified pending chemistry mapping'
      chem.append({'guide_id':gid,'chemistry_id':label,'recipe_status':'REFERENCE' if label.startswith('UNMOD') else 'EXPLORATORY_NOT_VALIDATED','notes':notes,'can_run_modified_MD':'NO: requires exact structure, PS stereochemistry and parameters' if not label.startswith('UNMOD') else 'CONDITIONAL: only after 3D coordinates and physical validation'})
     if fragment!=wfrag:
      conflicts.append({'guide_id':gid,'mutant_target':fragment,'wt_target':wfrag,'interpretation':'Guide can also hybridize WT; allele-specific editing NOT established'})
 writecsv(out/'01_variants.csv',ann,['variant_id','hgvs_c','protein_changes','edited_cds_positions_1based','directly_correctable','other_uncorrected_changes_after_single_edit','transcript'])
 gc_cols=['guide_id','variant_id','target_cds_position','guide_length','target_cds_start_1based','target_rna_5to3','guide_rna_5to3','target_index_0based','guide_index_0based','gc_percent','wt_local_mismatches','wt_guide_target_base','mutant_guide_target_base','sequence_basis','wt_site_pairing']
 writecsv(out/'02_guides.csv',guides,gc_cols)
 writecsv(out/'03_bystander_sites.csv',bystanders,['guide_id','cds_position','relative_to_target','is_target','context_5prime','context_3prime','guide_counterbase'])
 writecsv(out/'04_chemistry_manifest.csv',chem,['guide_id','chemistry_id','recipe_status','notes','can_run_modified_MD'])
 writecsv(out/'05_wt_comparison.csv',conflicts,['guide_id','mutant_target','wt_target','interpretation'])
 folds=[]
 if args.fold:
  for g in guides:
   status,val=fold(g['guide_rna_5to3']);folds.append({'guide_id':g['guide_id'],'RNAfold_status':status,'structure_energy':val})
 writecsv(out/'06_secondary_structure.csv',folds,['guide_id','RNAfold_status','structure_energy'])
 # Never silently fabricate structural predictions or transcriptome-scale off-target calculations.
 offstatus='NOT_RUN'
 if args.transcriptome:
  offstatus='INPUT_REGISTERED_NOT_SCORED: transcriptome aligner & editing-context model required; see README'
 (out/'07_offtarget_status.json').write_text(json.dumps({'status':offstatus,'transcriptome':args.transcriptome or None,'note':'No off-target clearance claim is produced.'},indent=2))
 (out/'08_md_input_manifest.csv').write_text('guide_id,structure_file,forcefield_params,stereochemistry,ad ar_complex_template,readiness\n'.replace('ad ar','adar')+''.join(f'{g["guide_id"]},MISSING,MISSING,MISSING,MISSING,BLOCKED\n' for g in guides))
 manifest={'timestamp_utc':__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),'transcript':args.transcript,'cds_source':args.cds_fasta or 'Ensembl REST cached','variants':len(panel),'guides':len(guides),'correctable_variants':sum(a['directly_correctable']=='YES' for a in ann),'status':'COMPUTATIONAL_CANDIDATE_GENERATION_ONLY','not_completed':['Transcriptome off-target alignment','Experimental calibrated editing efficiency prediction','Exact RESTORE2 chemistry transfer','Molecular dynamics','ADAR catalytic dynamics','Cellular delivery','Tumor efficacy']}
 (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2))
 report=['# KRAS ADAR computational report','',f'Variants: {len(panel)}; designed guides: {len(guides)}; directly correctable: {manifest["correctable_variants"]}.','','## Blocking gates','- Confirm Ensembl transcript version/MANE and clinical variant HGVS before experimental work.','- WT allele can hybridize: WT discrimination has NOT been proven.','- Bystander sites tabulated, NOT experimentally predicted editing rates.','- Off-target clearance needs transcriptome-scale assessment; not executed automatically.','- Chemical pattern is an exploratory placeholder; published RESTORE2 pattern MUST be transcribed position by position.','- MD blocked until validated 3D structures, chemistry and force-field parameters supplied.','- Cellular editing and KRAS signaling must be measured experimentally.','', '## Variant results']
 report += [f'- {a["variant_id"]}: {a["hgvs_c"]}; {a["protein_changes"]}; editable CDS bases {a["edited_cds_positions_1based"] or "none"}' for a in ann]
 (out/'SUMMARY.md').write_text('\n'.join(report)+'\n')
 print(json.dumps(manifest,indent=2))
def main():
 p=argparse.ArgumentParser();p.add_argument('--variants',default='config/variants.csv');p.add_argument('--transcript',default='ENST00000256078');p.add_argument('--cds-fasta');p.add_argument('--out',default='results/run');p.add_argument('--lengths',type=int,nargs='+',default=[33,45,59]);p.add_argument('--fold',action='store_true');p.add_argument('--transcriptome');a=p.parse_args()
 try:run(a)
 except Exception as e:print('FATAL:',e,file=sys.stderr);sys.exit(2)
if __name__=='__main__':main()
