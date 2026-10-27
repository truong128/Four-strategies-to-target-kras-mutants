#!/usr/bin/env python
import argparse,os,re,subprocess,tempfile,shutil,pandas as pd
from common import load_yaml,read_any_table
ap=argparse.ArgumentParser()
for x in ["input","receptor","config","output"]:ap.add_argument("--"+x,required=True)
ap.add_argument("--center",nargs=3,type=float,required=True);ap.add_argument("--size",nargs=3,type=float,required=True)
a=ap.parse_args();C=load_yaml(a.config)["gnina"];exe=C.get("exe","gnina")
if shutil.which(exe) is None:raise SystemExit("GNINA not found")
df=read_any_table(a.input);n=int(C["max_input"])
if len(df)>n:df=df.nlargest(n,"ml_probability") if "ml_probability" in df else df.head(n)
rx=re.compile(r"CNNscore[:\s=]+(-?\d+(?:\.\d+)?)",re.I);ra=re.compile(r"CNNaffinity[:\s=]+(-?\d+(?:\.\d+)?)",re.I);O=[]
with tempfile.TemporaryDirectory() as td:
  for r in df.itertuples(index=False):
    sdf=os.path.join(td,"l.sdf");out=os.path.join(td,"o.sdf")
    if subprocess.run(["python",os.path.join(os.path.dirname(__file__),"smiles_to_sdf.py"),"--smiles",r.smiles,"--id",str(r.id),"--output",sdf]).returncode:continue
    q=subprocess.run([exe,"-r",a.receptor,"-l",sdf,"-o",out,"--center_x",str(a.center[0]),"--center_y",str(a.center[1]),"--center_z",str(a.center[2]),
      "--size_x",str(a.size[0]),"--size_y",str(a.size[1]),"--size_z",str(a.size[2]),"--cnn_scoring",C.get("cnn_scoring","rescore")],text=True,capture_output=True)
    t=(q.stdout or "")+(q.stderr or "");m=rx.search(t);z=ra.search(t)
    O.append({"id":r.id,"smiles":r.smiles,"gnina_cnn_score":float(m.group(1)) if m else float("nan"),"gnina_cnn_affinity":float(z.group(1)) if z else float("nan")})
pd.DataFrame(O).to_csv(a.output,index=False)
