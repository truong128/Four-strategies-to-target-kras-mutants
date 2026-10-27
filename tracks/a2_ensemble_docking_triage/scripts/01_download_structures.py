#!/usr/bin/env python3
from pathlib import Path
import csv, requests
from common import ensure_dirs
root=ensure_dirs(); out=root/'work'/'raw_pdb'
rows=list(csv.DictReader(open(root/'config'/'pdb_panel.tsv'),delimiter='\t'))
ids=sorted({r['pdb_id'] for r in rows if not r['pdb_id'].startswith('MODEL_')})
for pdb in ids:
    url=f'https://files.rcsb.org/download/{pdb}.pdb'
    dest=out/f'{pdb}.pdb'
    if dest.exists() and dest.stat().st_size>1000:
        print('exists',pdb); continue
    print('download',pdb,url)
    r=requests.get(url,timeout=60); r.raise_for_status(); dest.write_bytes(r.content)
# MRTX1133 CCD ligand and GDP for reference/optional downstream use
for ccd in ['6IC','GDP']:
    url=f'https://files.rcsb.org/ligands/download/{ccd}_ideal.sdf'
    dest=root/'work'/'ligands'/'sdf'/f'{ccd}_ideal.sdf'
    r=requests.get(url,timeout=60); r.raise_for_status(); dest.write_bytes(r.content)
print('Downloaded:',', '.join(ids))
