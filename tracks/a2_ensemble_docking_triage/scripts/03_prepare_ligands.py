#!/usr/bin/env python3
from pathlib import Path
import csv, re, subprocess, sys
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas as pd
from common import ensure_dirs
root=ensure_dirs(); inp=root/'inputs'/'library.csv'; sdfdir=root/'work'/'ligands'/'sdf'; pqdir=root/'work'/'ligands'/'pdbqt'
rows=[]
if inp.exists():
    rows=list(csv.DictReader(open(inp)))
for old in list(pqdir.glob('*.pdbqt'))+[f for f in sdfdir.glob('*.sdf') if not f.name.endswith('_ideal.sdf')]:
    old.unlink()

def to_pdbqt(sdf,out):
    out.unlink(missing_ok=True)
    subprocess.run(['obabel',str(sdf),'-O',str(out),'-p','7.4','--partialcharge','gasteiger'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    ok=out.exists() and any(ln.startswith(('ATOM','HETATM')) for ln in out.read_text().splitlines())
    if not ok: out.unlink(missing_ok=True)
    return ok

# Always include crystallographic MRTX1133 positive control from RCSB CCD 6IC.
controls=[('MRTX1133_6IC',sdfdir/'6IC_ideal.sdf')]
meta=[]
for name,sdf in controls:
    out=pqdir/f'{name}.pdbqt'
    if not to_pdbqt(sdf,out): sys.exit('Could not convert the 6IC positive control to PDBQT')
    meta.append({'id':name,'smiles':'','source':'RCSB_6IC','MW':None,'cLogP':None,'TPSA':None,'QED':None})

uncharger=rdMolStandardize.Uncharger()
seen={m['id'] for m in meta}
for r in rows:
    raw_name=(r.get('id') or '').strip(); smi=(r.get('smiles') or '').strip()
    if not raw_name or not smi: continue
    name=re.sub(r'[^A-Za-z0-9_.-]','_',raw_name).lstrip('.')
    if not name or name in seen:
        print('SKIP duplicate or reserved id',raw_name,file=sys.stderr); continue
    mol=Chem.MolFromSmiles(smi)
    if mol is None:
        print('SKIP invalid SMILES',name,file=sys.stderr); continue
    # largest fragment + neutralize where possible
    frags=Chem.GetMolFrags(mol,asMols=True,sanitizeFrags=True); mol=max(frags,key=lambda m:m.GetNumHeavyAtoms()); mol=uncharger.uncharge(mol)
    can=Chem.MolToSmiles(mol,isomericSmiles=True)
    mol=Chem.AddHs(mol)
    ps=AllChem.ETKDGv3(); ps.randomSeed=0xC0FFEE
    if AllChem.EmbedMolecule(mol,ps)!=0:
        print('SKIP embed failed',name,file=sys.stderr); continue
    try: AllChem.MMFFOptimizeMolecule(mol,maxIters=500)
    except Exception: pass
    sdf=sdfdir/f'{name}.sdf'; w=Chem.SDWriter(str(sdf)); w.write(mol); w.close()
    out=pqdir/f'{name}.pdbqt'
    if not to_pdbqt(sdf,out):
        print('SKIP PDBQT conversion failed',name,file=sys.stderr); sdf.unlink(missing_ok=True); continue
    seen.add(name)
    base=Chem.RemoveHs(mol)
    meta.append({'id':name,'smiles':can,'source':'user','MW':Descriptors.MolWt(base),'cLogP':Crippen.MolLogP(base),'TPSA':rdMolDescriptors.CalcTPSA(base),'QED':QED.qed(base)})
if len(meta)==len(controls):
    print('WARNING: no valid user compounds in inputs/library.csv; only the MRTX1133 control will be docked',file=sys.stderr)
pd.DataFrame(meta).to_csv(root/'results'/'tables'/'ligand_properties.csv',index=False)
print('Prepared ligands:',len(meta))
