#!/usr/bin/env python3
import shutil, sys
from common import ensure_dirs
root=ensure_dirs()
required=['obabel','vina']
missing=[x for x in required if not shutil.which(x)]
mods=['numpy','pandas','yaml','requests','matplotlib','tabulate','Bio','rdkit','pdbfixer','openmm']
bad=[]
for m in mods:
    try: __import__(m)
    except Exception as e: bad.append((m,str(e)))
print('Project root:',root)
print('Executables:',{x:shutil.which(x) for x in required})
if missing or bad:
    print('\nMISSING DEPENDENCIES')
    if missing: print('Executables:',missing)
    for x,e in bad: print(x,e)
    print('\nCreate environment with: conda env create -f environment.yml && conda activate krasg12d_e2e')
    sys.exit(2)
print('Environment check: PASS')
