#!/usr/bin/env python3
"""Align experimental structures to 7RPZ, model selected secondary mutants from 7RPZ,
strip bound screening ligand while preserving GDP/GTP-analog/Mg, and compute docking center.
Modeled mutants are stress-test models, NOT experimental structures.
"""
from pathlib import Path
import csv, json, re
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Superimposer, Select
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config(); raw=root/'work'/'raw_pdb'; aligned=root/'work'/'aligned_pdb'
parser=PDBParser(QUIET=True)
ref=parser.get_structure('ref',str(raw/'7RPZ.pdb'))

# get CA atom mapping by residue number for robust alignment

def ca_map(struct, chain='A'):
    d={}
    model=next(struct.get_models())
    if chain not in model: chain=list(model.child_dict)[0]
    for res in model[chain]:
        if res.id[0]==' ' and 'CA' in res:
            d[res.id[1]]=res['CA']
    return d
refmap=ca_map(ref,'A')

def align_pdb(src,dst):
    s=parser.get_structure('x',str(src)); m=ca_map(s,'A')
    common=sorted(set(refmap)&set(m))
    common=[i for i in common if 5<=i<=165]
    sup=Superimposer(); sup.set_atoms([refmap[i] for i in common],[m[i] for i in common])
    sup.apply(list(s.get_atoms()))
    io=PDBIO(); io.set_structure(s); io.save(str(dst))
    return float(sup.rms), len(common)

rows=list(csv.DictReader(open(root/'config'/'pdb_panel.tsv'),delimiter='\t'))
metrics=[]
for r in rows:
    pid=r['pdb_id']
    if pid.startswith('MODEL_'): continue
    rms,n=align_pdb(raw/f'{pid}.pdb',aligned/f'{pid}_aligned.pdb')
    metrics.append({'id':pid,'rmsd_to_7RPZ_CA':rms,'n_ca':n})

# model mutations from G12D template. Use PDBFixer on protein; then merge cofactors later.
model_specs={
'MODEL_V9E':'VAL-9-GLU','MODEL_V9W':'VAL-9-TRP','MODEL_V9Q':'VAL-9-GLN',
'MODEL_G13P':'GLY-13-PRO','MODEL_T58Y':'THR-58-TYR','MODEL_R68G':'ARG-68-GLY',
'MODEL_Y96W':'TYR-96-TRP','MODEL_H95Q':'HIS-95-GLN','MODEL_H95L':'HIS-95-LEU'}
# cofactor lines from 7RPZ: GDP/MG only, no MRTX1133 (6IC)
ref_lines=(raw/'7RPZ.pdb').read_text().splitlines()
cof=[ln for ln in ref_lines if ln.startswith('HETATM') and ln[17:20].strip() in {'GDP','MG'}]
for mid,mut in model_specs.items():
    fixer=PDBFixer(filename=str(raw/'7RPZ.pdb'))
    fixer.removeHeterogens(keepWater=False)
    fixer.applyMutations([mut],'A')
    fixer.findMissingResidues(); fixer.findMissingAtoms(); fixer.addMissingAtoms(); fixer.addMissingHydrogens(7.4)
    tmp=aligned/f'{mid}_protein.pdb'
    with open(tmp,'w') as f: PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
    lines=[ln for ln in tmp.read_text().splitlines() if not ln.startswith('END')]
    (aligned/f'{mid}_aligned.pdb').write_text('\n'.join(lines+cof+['END'])+'\n')
    tmp.unlink(missing_ok=True)

# find 6IC center in reference aligned structure
coords=[]
for ln in (aligned/'7RPZ_aligned.pdb').read_text().splitlines():
    if ln.startswith('HETATM') and ln[17:20].strip()=='6IC':
        elem=ln[76:78].strip().upper()
        if elem!='H': coords.append([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
if not coords: raise RuntimeError('Could not find 6IC in 7RPZ')
center=np.asarray(coords).mean(axis=0)
(root/'results'/'grid.json').write_text(json.dumps({'center_x':float(center[0]),'center_y':float(center[1]),'center_z':float(center[2]),**cfg['grid']},indent=2))

# receptor PDBs: strip water and organic inhibitors; preserve nucleotide and Mg.
keep_hets={'GDP','GCP','GNP','GTP','MG'}
nucleotides=keep_hets-{'MG'}
first_altloc={' ','A'}
def xyz(ln): return np.array([float(ln[30:38]),float(ln[38:46]),float(ln[46:54])])
for r in rows:
    pid=r['pdb_id']; src=aligned/f'{pid}_aligned.pdb'; dst=root/'work'/'receptors'/f'{pid}.pdb'
    lines=src.read_text().splitlines()
    nuc_xyz=[xyz(ln) for ln in lines if ln.startswith('HETATM') and ln[17:20].strip() in nucleotides
             and ln[21:22] in {'A',' '} and ln[16:17] in first_altloc]
    out=[]
    for ln in lines:
        if ln.startswith('ATOM') and ln[21:22]=='A' and ln[16:17] in first_altloc:
            out.append(ln[:16]+' '+ln[17:])
        elif ln.startswith('HETATM') and ln[17:20].strip() in keep_hets and ln[21:22] in {'A',' '} and ln[16:17] in first_altloc:
            # keep only the Mg2+ coordinated by the nucleotide (<=3.0 A), drop lattice/crystallization Mg
            if ln[17:20].strip()=='MG' and not any(np.linalg.norm(xyz(ln)-q)<=3.0 for q in nuc_xyz): continue
            out.append(ln[:16]+' '+ln[17:])
        elif ln.startswith(('TER','END')):
            out.append(ln)
    dst.write_text('\n'.join(out)+'\n')

import pandas as pd
pd.DataFrame(metrics).to_csv(root/'results'/'tables'/'alignment_metrics.csv',index=False)
print('Grid center:',center)
print('Built receptor panel:',len(rows),'structures')
