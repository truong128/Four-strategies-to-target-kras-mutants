#!/usr/bin/env python
"""Dock many ligands at once: one Vina process per worker, one core each."""
import argparse, os, re, shutil, subprocess, tempfile, multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from common import load_yaml, read_any_table

SCORE_RE=re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", re.M)

def _embed(smiles, mid, sdf, optimize):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m=Chem.MolFromSmiles(str(smiles))
    if m is None: return False
    m=Chem.AddHs(m); p=AllChem.ETKDGv3(); p.randomSeed=42
    if AllChem.EmbedMolecule(m, p)!=0: return False
    if optimize:
        try: AllChem.MMFFOptimizeMolecule(m, maxIters=50)
        except Exception:
            try: AllChem.UFFOptimizeMolecule(m, maxIters=50)
            except Exception: pass
    m.SetProp("_Name", str(mid))
    w=Chem.SDWriter(sdf); w.write(m); w.close()
    return True

def _dock_one(task):
    smiles, mid, extra, receptor, center, size, exh, modes, vina, meeko, optimize = task
    td=tempfile.mkdtemp(prefix="vina_")
    try:
        sdf=os.path.join(td, "lig.sdf"); pq=os.path.join(td, "lig.pdbqt"); outp=os.path.join(td, "out.pdbqt")
        if not _embed(smiles, mid, sdf, optimize):
            status="embed_fail"; score=float("nan")
        else:
            q=subprocess.run([meeko, "-i", sdf, "-o", pq], text=True, capture_output=True)
            if q.returncode:
                status="meeko_fail"; score=float("nan")
            else:
                cmd=[vina, "--receptor", receptor, "--ligand", pq,
                     "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
                     "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
                     "--exhaustiveness", str(exh), "--num_modes", str(modes), "--cpu", "1", "--out", outp]
                q=subprocess.run(cmd, text=True, capture_output=True)
                m=SCORE_RE.search((q.stdout or "")+"\n"+(q.stderr or ""))
                score=float(m.group(1)) if m else float("nan")
                status="ok" if m else "parse_fail"
        row={"smiles": smiles, "id": mid, "docking_score": score, "dock_status": status}
        row.update(extra)
        return row
    finally:
        shutil.rmtree(td, ignore_errors=True)

def main():
    ap=argparse.ArgumentParser()
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--smiles"); g.add_argument("--csv")
    ap.add_argument("--receptor", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--config", required=True); ap.add_argument("--center", nargs=3, type=float, required=True)
    ap.add_argument("--size", nargs=3, type=float, required=True); ap.add_argument("--redock", action="store_true")
    a=ap.parse_args(); c=load_yaml(a.config); dc=c["docking"]
    vina=dc.get("vina_exe", "vina"); meeko=dc.get("meeko_ligand_exe", "mk_prepare_ligand.py")
    for x in [vina, meeko]:
        if shutil.which(x) is None: raise SystemExit(f"Missing executable: {x}")
    df=read_any_table(a.smiles or a.csv)
    if "smiles" not in df or "id" not in df: raise SystemExit("Need smiles and id")
    if a.redock:
        n=int(c["redock"]["max_hits"])
        if len(df)>n:
            df=df.nlargest(n, "ml_probability") if "ml_probability" in df.columns else df.head(n)
    exh=dc.get("redock_exhaustiveness", 8) if a.redock else dc.get("exhaustiveness", 4)
    modes=int(dc.get("num_modes", 1)); optimize=bool(dc.get("embed_optimize", False))
    workers=int(dc.get("workers") or 16)
    workers=max(1, min(workers, len(df) or 1))
    extra_cols=[col for col in ["ml_probability", "p0", "p1", "deltaP"] if col in df.columns]
    tasks=[]
    for r in df.itertuples(index=False):
        extra={col: getattr(r, col) for col in extra_cols}
        tasks.append((r.smiles, r.id, extra, os.path.abspath(a.receptor), tuple(a.center), tuple(a.size),
                      exh, modes, vina, meeko, optimize))
    rows=[]
    if workers==1:
        rows=[_dock_one(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            for row in pool.map(_dock_one, tasks, chunksize=8):
                rows.append(row)
    pd.DataFrame(rows).to_csv(a.output, index=False)
    ok=sum(1 for r in rows if r["dock_status"]=="ok")
    print(f"docked {ok}/{len(rows)} workers {workers} exhaustiveness {exh}")

if __name__=="__main__":
    main()
