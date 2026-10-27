#!/usr/bin/env python3
from pathlib import Path
import csv, subprocess, sys
from common import ensure_dirs
root=ensure_dirs(); rows=list(csv.DictReader(open(root/'config'/'pdb_panel.tsv'),delimiter='\t'))
failed=[]
for r in rows:
    pid=r['pdb_id']; src=root/'work'/'receptors'/f'{pid}.pdb'; out=root/'work'/'receptors'/f'{pid}.pdbqt'
    out.unlink(missing_ok=True)
    # OpenBabel is used for a portable, reproducible first-pass receptor preparation.
    # Final leads should be re-prepared/protonated with a higher-fidelity workflow and inspected manually.
    # No --partialcharge here: OpenBabel Gasteiger has no Mg parameters and aborts the whole conversion
    # (exit code 0, empty output); Vina does not use PDBQT charges.
    cmd=['obabel',str(src),'-O',str(out),'-xr','-p','7.4']
    subprocess.run(cmd,check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    if not out.exists() or sum(1 for ln in out.read_text().splitlines() if ln.startswith(('ATOM','HETATM')))==0:
        failed.append(pid)
if failed:
    sys.exit(f'Receptor PDBQT conversion produced no atoms for: {", ".join(failed)}')
print('Prepared receptor PDBQT panel:',len(rows))
