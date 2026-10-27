#!/usr/bin/env python3
import shutil, subprocess, pandas as pd
from common import ensure_dirs, load_config
root=ensure_dirs(); cfg=load_config()
(root/'results'/'tables'/'admet_ai.csv').unlink(missing_ok=True)
props=pd.read_csv(root/'results'/'tables'/'ligand_properties.csv')
sm=props[(props.source=='user') & props.smiles.notna() & (props.smiles!='')][['id','smiles']]
if not cfg['admet']['enabled'] or sm.empty:
    print('ADMET: skipped (disabled or no user ligands)'); raise SystemExit(0)
inp=root/'work'/'admet_input.csv'; out=root/'results'/'tables'/'admet_ai.csv'; sm.to_csv(inp,index=False)
if not shutil.which('admet_predict'):
    print('ADMET-AI CLI not found; install with pip install admet-ai'); raise SystemExit(0)
subprocess.run(['admet_predict','--data_path',str(inp),'--save_path',str(out),'--smiles_column','smiles'],check=True)
print('ADMET predictions:',out)
