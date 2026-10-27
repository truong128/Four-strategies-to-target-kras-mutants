#!/usr/bin/env python3
"""Evidence-gated full-workflow audit. Never labels hypothetical results validated."""
import argparse, csv, hashlib, json, shutil, subprocess, sys
from pathlib import Path
from datetime import datetime, timezone

TOOLS={
 'RNAfold':('viennarna','RNA secondary structure, unmodified approximation','https://github.com/ViennaRNA/ViennaRNA'),
 'RNAcofold':('viennarna','target/guide cofold, unmodified approximation','https://github.com/ViennaRNA/ViennaRNA'),
 'bowtie2':('bowtie2','optional alignment; not an editing predictor','https://github.com/BenLangmead/bowtie2'),
 'minimap2':('minimap2','optional transcript/genome mapping','https://github.com/lh3/minimap2'),
 'samtools':('samtools','BAM operations','https://github.com/samtools/samtools'),
 'bcftools':('bcftools','VCF processing','https://github.com/samtools/bcftools'),
 'vep':('ensembl-vep','variant annotation with human GRCh38 cache','https://github.com/Ensembl/ensembl-vep'),
 'RNAhybrid':('viennarna-optional','RNA hybrid screening; distinct from editing probability','https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid'),
 'tleap':('ambertools','Amber topology setup; custom chemistry parameters essential','https://ambermd.org/AmberTools.php'),
 'cpptraj':('ambertools','MD trajectory analysis','https://ambermd.org/AmberTools.php'),
 'pmemd.cuda':('amber-gpu','production GPU MD, separate licensed installation','https://ambermd.org/'),
 'sander':('ambertools','alternative MD engine','https://ambermd.org/AmberTools.php'),
 'java':('java','JACUSA2 runtime for measured RNA editing','https://github.com/dieterich-lab/JACUSA2'),
}
PROJECTS={
 'AdarEdit':{'url':'https://github.com/Scientific-Computing-Lab/AdarEdit','role':'natural Alu-ADAR editability baseline only, not RESTORE efficacy'},
 'EditPredict':{'url':'https://github.com/wjd198605/EditPredict','role':'A-to-I baseline, out of domain for modified guides'},
 'JACUSA2':{'url':'https://github.com/dieterich-lab/JACUSA2','role':'RNA-seq editing site calls from BAM; needs matched controls and DNA/variant filters'},
 'REDItools2':{'url':'https://github.com/BioinfoUNIBA/REDItools2','role':'independent RNA editing calling and controls'},
 'ViennaRNA':{'url':'https://github.com/ViennaRNA/ViennaRNA','role':'unmodified RNA secondary structure reference'},
 'EnsemblVEP':{'url':'https://github.com/Ensembl/ensembl-vep','role':'somatic VCF gene/transcript consequences'},
}
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for x in iter(lambda:f.read(1024*1024),b''):h.update(x)
 return h.hexdigest()
def readrows(p):
 if not p.exists():return []
 with open(p,newline='') as f:return list(csv.DictReader(f))
def audit(base, out):
 base=Path(base);out=Path(out);out.mkdir(parents=True,exist_ok=True)
 tool_rows=[]
 for exe,(package,scope,url) in TOOLS.items():
  found=shutil.which(exe)
  tool_rows.append({'tool':exe,'installed':bool(found),'executable':found,'purpose':scope,'source':url})
 with open(out/'tool_inventory.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(tool_rows[0]));w.writeheader();w.writerows(tool_rows)
 specs=[
 ('01','KRAS CDS & HGVS','01_variants.csv','CODED','Reference transcript version, allele phase, tissue VCF still require patient evidence'),
 ('02','Guide design','02_guides.csv','CODED','Allele-specificity not established'),
 ('03','Bystanders','03_bystander_sites.csv','ENUMERATED_ONLY','Not measured bystander editing'),
 ('04','Published chemistry transfer','restore2_published/transferred_per_residue.csv','REFERENCE_DERIVED','Patterns transferred; original target assay not KRAS and PS stereo-random'),
 ('05','RNA secondary structure','08_cofold.csv','CONDITIONAL','Unmodified thermodynamic approximation'),
 ('06','Transcriptome binding','09_transcriptome_hits.csv','CONDITIONAL','Need expressed human transcriptome and expression evidence'),
 ('07','ADAR1/2 biological context','not_implemented','NOT_AVAILABLE','PDAC ADAR expression, localization and IFN status required'),
 ('08','Protein-RNA complex','not_implemented','NOT_AVAILABLE','Experimental template, residue mapping, modified force field'),
 ('09','Production MD','10_md_readiness.csv','GATED','No claim of run unless validated topology and executed trajectories'),
 ('10','Editing ML','11_editing_ml.csv','GATED','No independently validated RESTORE-on-KRAS efficiency estimator'),
 ('11','Measured transcriptome editing','not_implemented','EXTERNAL_EXPERIMENT','JACUSA2/REDItools2 need treated/control RNA-seq plus DNA evidence'),
 ('12','Therapeutic ranking','12_candidate_evidence.csv','NOT_VALIDATED','Do not substitute sequence/MD scores for tumor efficacy'),
 ]
 rows=[]
 for num,label,path,code_status,note in specs:
  p=base/path
  rows.append({'stage':num,'stage_name':label,'code_status':code_status,'output_exists':p.exists(),'output_rows':len(readrows(p)) if p.suffix=='.csv' and p.exists() else '', 'limitation':note})
 with open(out/'stage_matrix.csv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 assets={}
 for k,p in {'supplement':'source_material/41467_2025_64434_MOESM3_ESM.xlsx','source_data':'source_material/41467_2025_64434_MOESM6_ESM.xlsx','constructs':'reference/publisher/publisher_constructs.csv','residues':'reference/publisher/publisher_residues.csv'}.items():
  pp=Path(p);assets[k]={'exists':pp.exists(),'sha256':digest(pp) if pp.exists() else None}
 manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'status':'RESEARCH_PIPELINE_AUDIT_NOT_A_THERAPEUTIC_VALIDATION','assets':assets,'projects':PROJECTS,'required_to_claim_full_E2E':['Patient-specific VCF/DNA and matched tumor RNA','ADAR protein expression in PDAC model','Full transcriptome/WT off-target study with controls','Exact verified chemically modified structures and force-field validation','ADAR-RNA complex templates and independent MD replicas','Measured KRAS G12D editing with allele discrimination','Functional KRAS-GTP/pERK and viability','Delivery and toxicity in PDAC models']}
 (out/'audit.json').write_text(json.dumps(manifest,indent=2))
 (out/'AUDIT_REPORT.md').write_text('# KRAS-ADAR Pipeline Audit\n\n**No therapeutic efficacy claim is supported by sequence screening or MD alone.**\n\n'+ '\n'.join(f'- **{r["stage"]} {r["stage_name"]}**: {r["code_status"]}; output={r["output_exists"]}. {r["limitation"]}' for r in rows)+'\n\n## Additional verified tools\n'+ '\n'.join(f'- [{k}]({v["url"]}): {v["role"]}' for k,v in PROJECTS.items())+'\n')
 print(json.dumps({'tool_inventory':len(tool_rows),'stages':len(rows),'report':str(out/'AUDIT_REPORT.md')},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--results',default='results/run');p.add_argument('--out',default='results/audit');a=p.parse_args();audit(a.results,a.out)
