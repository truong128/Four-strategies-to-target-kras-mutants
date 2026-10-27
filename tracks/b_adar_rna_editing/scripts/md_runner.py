#!/usr/bin/env python3
"""Optional short explicit-solvent MD for PREVALIDATED Amber topologies only.
Requires independently verified equilibrated solvated topology/coordinates.
Not a de novo MD setup, editing-efficiency model or catalytic simulation.
"""
import argparse,csv,json,subprocess,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--out',default='results/md');p.add_argument('--steps',type=int,default=500000);p.add_argument('--engine',choices=['pmemd.cuda','pmemd','sander'],default='sander');p.add_argument('--execute',action='store_true');args=p.parse_args()
if args.steps<100 or args.steps>50000000:raise ValueError('Steps outside safety range')
engine=shutil.which(args.engine)
if args.execute and not engine:raise RuntimeError(f'{args.engine} missing on PATH')
rows=list(csv.DictReader(open(args.manifest)));out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
report=[]
for r in rows:
 id=r['construct_id']; keys=['topology_prmtop','coordinates_inpcrd','structure_pdb']
 metadata=['chemistry_verified','stereochemistry_verified','parameterization_validated','template_validated','equilibrated_explicit_solvent_ions']
 valid=all(r.get(k) and Path(r[k]).exists() for k in keys) and all(r.get(k,'').upper()=='YES' for k in metadata)
 status='BLOCKED_NEEDS_VALIDATED_EQUILIBRATED_SOLVATED_AMBER_INPUT' if not valid else 'DRY_RUN_READY'
 if valid:
  folder=out/id;folder.mkdir(parents=True,exist_ok=True)
  # NPT restart, 2-fs timestep, hydrogen bond constraints; depends on periodic pre-equilibrated system
  inp=folder/'production.in'
  inp.write_text(f'''Pre-equilibrated explicit solvent MD; chemistry prevalidated\n &cntrl\n  imin=0,irest=1,ntx=5,ntb=2,ntp=1,\n  nstlim={args.steps},dt=0.002,ntc=2,ntf=2,\n  tempi=310.15,temp0=310.15,ntt=3,gamma_ln=2.0,\n  cut=9.0,ntpr=5000,ntwx=5000,ntwr=5000,\n  ioutfm=1,ntxo=2,ig=-1,\n /\n''')
  cmd=[args.engine,'-O','-i',str(inp.resolve()),'-o',str((folder/'production.out').resolve()),'-p',str(Path(r['topology_prmtop']).resolve()),'-c',str(Path(r['coordinates_inpcrd']).resolve()),'-r',str((folder/'production.rst7').resolve()),'-x',str((folder/'production.nc').resolve())]
  (folder/'command.json').write_text(json.dumps(cmd,indent=2))
  if args.execute:
   with (folder/'execution.log').open('w') as log:
    finished=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   status='COMPLETED_CHECK_TRAJECTORY_CONVERGENCE' if finished.returncode==0 else 'FAILED_'+str(finished.returncode)
  else:status='DRY_RUN_READY_COMMAND_WRITTEN'
 report.append(dict(construct_id=id,status=status))
with (out/'md_execution_report.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['construct_id','status']);w.writeheader();w.writerows(report)
print(json.dumps({'constructs':len(rows),'executed':args.execute,'results':report},indent=2))
