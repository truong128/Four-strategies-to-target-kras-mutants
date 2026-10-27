#!/usr/bin/env python
"""Single-frame MM/GBSA for the same short list used by the MD stage."""
import argparse, os, re, shutil, subprocess, sys
import pandas as pd
from common import load_yaml, read_any_table

LEAP = """source leaprc.protein.ff14SB
source leaprc.gaff2
REC = loadpdb receptor.pdb
LIG = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
COM = combine {REC LIG}
saveamberparm REC receptor.prmtop receptor.inpcrd
saveamberparm LIG ligand.prmtop ligand.inpcrd
saveamberparm COM complex.prmtop complex.inpcrd
quit
"""

def amber_env(cfg):
    exe=cfg.get("exe", "MMPBSA.py")
    bindir=os.path.dirname(exe) if os.path.isabs(exe) else ""
    env=os.environ.copy()
    if bindir:
        env["PATH"]=bindir+os.pathsep+env.get("PATH", "")
    return env

def one(row, protein, cfg, work, env):
    tag=str(row.id).replace("/","_")
    d=os.path.abspath(os.path.join(work, tag)); os.makedirs(d, exist_ok=True)
    sdf=os.path.join(d, "lig.sdf")
    pose=getattr(row, "pose_sdf", None)
    if isinstance(pose, str) and os.path.isfile(pose):
        shutil.copy(pose, sdf)
    else:
        q=subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "smiles_to_sdf.py"),
                          "--smiles", str(row.smiles), "--id", tag, "--output", sdf], cwd=d)
        if q.returncode: return float("nan"), "embed_fail"
    ante=cfg.get("antechamber", "antechamber"); tleap=cfg.get("tleap", "tleap")
    q=subprocess.run([ante, "-i", sdf, "-fi", "sdf", "-o", "ligand.mol2", "-fo", "mol2", "-c", "bcc", "-at", "gaff2", "-s", "2", "-pf", "y"],
                     cwd=d, env=env, capture_output=True, text=True)
    if q.returncode: return float("nan"), "antechamber_fail"
    parmchk=os.path.join(os.path.dirname(ante), "parmchk2") if os.path.isabs(ante) else "parmchk2"
    subprocess.run([parmchk, "-i", "ligand.mol2", "-f", "mol2", "-o", "ligand.frcmod"], cwd=d, env=env, capture_output=True)
    subprocess.run(["cp", protein, os.path.join(d, "receptor.pdb")])
    open(os.path.join(d, "leap.in"), "w").write(LEAP)
    q=subprocess.run([tleap, "-f", "leap.in"], cwd=d, env=env, capture_output=True, text=True)
    if q.returncode or not os.path.isfile(os.path.join(d, "complex.prmtop")): return float("nan"), "tleap_fail"
    cpptraj=os.path.join(os.path.dirname(tleap), "cpptraj")
    open(os.path.join(d, "cpp.in"), "w").write("parm complex.prmtop\ntrajin complex.inpcrd\ntrajout complex.nc netcdf\nrun\nquit\n")
    subprocess.run([cpptraj, "-i", "cpp.in"], cwd=d, env=env, capture_output=True)
    inp=cfg.get("input_file", "mmpbsa.in")
    if not os.path.isabs(inp):
        inp=os.path.join(os.path.dirname(__file__), "..", inp)
    q=subprocess.run([cfg.get("exe", "MMPBSA.py"), "-O", "-i", os.path.abspath(inp),
                      "-cp", "complex.prmtop", "-rp", "receptor.prmtop", "-lp", "ligand.prmtop",
                      "-y", "complex.nc", "-o", "FINAL_MMPBSA.dat"], cwd=d, env=env, capture_output=True, text=True)
    text=open(os.path.join(d, "FINAL_MMPBSA.dat")).read() if os.path.isfile(os.path.join(d, "FINAL_MMPBSA.dat")) else (q.stdout or "")
    m=re.search(r"DELTA TOTAL\s+(-?\d+(?:\.\d+)?)", text)
    if not m: return float("nan"), "parse_fail"
    return float(m.group(1)), "ok"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--protein", required=True)
    ap.add_argument("--config", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--workdir", required=True)
    a=ap.parse_args(); C=load_yaml(a.config)["mmpbsa"]; env=amber_env(C)
    if not C.get("enabled", True):
        df=read_any_table(a.input)
        df["mmgbsa"]=float("nan"); df["mmgbsa_status"]="skipped: mmpbsa disabled"
        df[["id","smiles","mmgbsa","mmgbsa_status"]].to_csv(a.output, index=False)
        print(df["mmgbsa_status"].iloc[0]); return
    if not os.path.isfile(C.get("exe", "")) and not shutil_which(C.get("exe", "MMPBSA.py"), env):
        df=read_any_table(a.input)
        df["mmgbsa"]=float("nan"); df["mmgbsa_status"]="skipped: MMPBSA.py not found"
        df[["id","smiles","mmgbsa","mmgbsa_status"]].to_csv(a.output, index=False)
        print(df["mmgbsa_status"].iloc[0]); return
    df=read_any_table(a.input)
    n=int(C.get("max_input", 20))
    if "docking_score" in df.columns:
        df=df.assign(_s=pd.to_numeric(df["docking_score"], errors="coerce")).sort_values("_s").head(n)
    else:
        df=df.head(n)
    rows=[]
    for r in df.itertuples(index=False):
        val, status=one(r, os.path.abspath(a.protein), C, a.workdir, env)
        rows.append({"id": r.id, "smiles": r.smiles, "mmgbsa": val, "mmgbsa_status": status})
    pd.DataFrame(rows).to_csv(a.output, index=False)
    print(f"mmpbsa_ok {sum(1 for r in rows if r['mmgbsa_status']=='ok')}/{len(rows)}")

def shutil_which(exe, env):
    import shutil
    return shutil.which(exe, path=env.get("PATH"))

if __name__=="__main__":
    main()
