#!/usr/bin/env python3
from pathlib import Path
import csv, json, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config(); grid=json.load(open(root/'results'/'grid.json'))
panel=list(csv.DictReader(open(root/'config'/'pdb_panel.tsv'),delimiter='\t'))
ligands=sorted((root/'work'/'ligands'/'pdbqt').glob('*.pdbqt'))
outdir=root/'work'/'docking'
missing=[p['pdb_id'] for p in panel if not (root/'work'/'receptors'/f"{p['pdb_id']}.pdbqt").exists()]
if missing: sys.exit(f'Missing receptor PDBQT files (run 04 first): {missing}')
if not ligands: sys.exit('No ligand PDBQT files found (run 03 first)')

def one(pid,lig):
    d=outdir/pid; d.mkdir(parents=True,exist_ok=True); stem=lig.stem
    out=d/f'{stem}_out.pdbqt'; log=d/f'{stem}.log'
    cmd=['vina','--receptor',str(root/'work'/'receptors'/f'{pid}.pdbqt'),'--ligand',str(lig),
         '--center_x',str(grid['center_x']),'--center_y',str(grid['center_y']),'--center_z',str(grid['center_z']),
         '--size_x',str(grid['size_x']),'--size_y',str(grid['size_y']),'--size_z',str(grid['size_z']),
         '--exhaustiveness',str(cfg['docking']['exhaustiveness']),'--num_modes',str(cfg['docking']['num_modes']),
         '--energy_range',str(cfg['docking']['energy_range']),'--cpu',str(cfg['docking']['cpu_per_job']),
         '--seed',str(cfg['docking']['seed']),'--out',str(out)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    log.write_text(p.stdout+'\n'+p.stderr)
    if p.returncode!=0: return {'pdb_id':pid,'ligand':stem,'affinity':None,'status':'FAIL'}
    m=re.search(r'^\s*1\s+(-?\d+\.\d+)',p.stdout,re.M)
    return {'pdb_id':pid,'ligand':stem,'affinity':float(m.group(1)) if m else None,'status':'OK' if m else 'PARSE_FAIL'}

jobs=[]
with ThreadPoolExecutor(max_workers=int(cfg['docking']['parallel_jobs'])) as ex:
    for p in panel:
        for lig in ligands: jobs.append(ex.submit(one,p['pdb_id'],lig))
    rows=[]
    for i,f in enumerate(as_completed(jobs),1):
        r=f.result(); rows.append(r); print(f'[{i}/{len(jobs)}]',r['pdb_id'],r['ligand'],r['affinity'])
res=pd.DataFrame(rows).sort_values(['pdb_id','ligand']).reset_index(drop=True)
res.to_csv(root/'results'/'tables'/'docking_long.csv',index=False)
n_ok=int((res.status=='OK').sum())
print(f'Docking finished: {n_ok}/{len(res)} OK, {len(res)-n_ok} failed (see work/docking/<pdb>/<ligand>.log)')
if n_ok==0: sys.exit('All docking jobs failed')
