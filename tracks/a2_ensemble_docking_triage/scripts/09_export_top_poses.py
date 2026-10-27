#!/usr/bin/env python3
"""Export top-ranked primary-target docked poses as PDB for visual inspection / later MD setup."""
from pathlib import Path
import subprocess, pandas as pd
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config(); r=pd.read_csv(root/'results'/'tables'/'ranked_hits.csv')
outdir=root/'results'/'top_pose_pdb'; outdir.mkdir(parents=True,exist_ok=True)
pid=cfg['ranking']['primary_gdp']
for lig in r.head(int(cfg['report']['top_n']))['ligand']:
    src=root/'work'/'docking'/pid/f'{lig}_out.pdbqt'
    if not src.exists(): continue
    dst=outdir/f'{lig}_{pid}.pdb'
    subprocess.run(['obabel',str(src),'-O',str(dst)],check=True,stdout=subprocess.DEVNULL)
print('Exported top poses to',outdir)
