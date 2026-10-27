#!/usr/bin/env python
"""Run OpenMM for the top docked compounds and record ligand RMSD."""
import argparse, os, shutil, subprocess, sys
import numpy as np, pandas as pd

def reexec(python):
    if python and os.path.isfile(python) and os.path.realpath(python)!=os.path.realpath(sys.executable):
        os.execv(python, [python, os.path.abspath(__file__), *sys.argv[1:]])

def heavy_coords(mol):
    from rdkit import Chem
    conf=mol.GetConformer()
    return np.array([conf.GetAtomPosition(a.GetIdx()) for a in mol.GetAtoms() if a.GetAtomicNum()>1], dtype=float)

def pdb_ca(path):
    pts=[]
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip()=="CA":
            pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(pts, dtype=float)

def kabsch(P, Q):
    cP, cQ=P.mean(0), Q.mean(0)
    H=(Q-cQ).T@(P-cP)
    U,_,Vt=np.linalg.svd(H)
    R=Vt.T@U.T
    if np.linalg.det(R)<0:
        Vt=Vt.copy(); Vt[-1]*=-1; R=Vt.T@U.T
    return R, cP, cQ

def ligand_rmsd(start_sdf, final_pdb, protein_pdb):
    from rdkit import Chem
    start=Chem.MolFromMolFile(start_sdf, removeHs=False) if start_sdf.endswith(".sdf") else None
    if start is None: return float("nan")
    s=heavy_coords(start)
    ca0, ca1=pdb_ca(protein_pdb), pdb_ca(final_pdb)
    n=min(len(ca0), len(ca1))
    if n<10 or len(s)<3: return float("nan")
    R, cP, cQ=kabsch(ca0[:n], ca1[:n])
    fin=Chem.MolFromPDBFile(final_pdb, removeHs=True, sanitize=False)
    if fin is None: return float("nan")
    # Ligand is the non-polymer residue added before solvent. Use the last residue that is not water or ion.
    skip={"HOH","WAT","NA","CL","K","MG"}
    residues={}
    for atom in fin.GetAtoms():
        info=atom.GetPDBResidueInfo()
        if info is None: continue
        name=info.GetResidueName().strip()
        if name in skip: continue
        residues.setdefault((info.GetChainId(), info.GetResidueNumber(), name), []).append(atom.GetIdx())
    het=[(k,v) for k,v in residues.items() if k[2] not in {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}]
    if not het: return float("nan")
    idx=het[-1][1]
    conf=fin.GetConformer()
    q=np.array([conf.GetAtomPosition(i) for i in idx if fin.GetAtomWithIdx(i).GetAtomicNum()>1], dtype=float)
    m=min(len(s), len(q))
    if m<3: return float("nan")
    aligned=((q[:m]-cQ)@R.T)+cP
    return float(np.sqrt(((aligned-s[:m])**2).sum(1).mean()))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--protein", required=True)
    ap.add_argument("--config", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--workdir", required=True)
    a=ap.parse_args()
    import yaml
    cfg=yaml.safe_load(open(a.config)); C=cfg["md"]
    if not C.get("enabled", True):
        from common import read_any_table
        df=read_any_table(a.input)
        df["md_ligand_rmsd"]=float("nan"); df["md_stability"]=float("nan"); df["md_status"]="skipped: md disabled"
        df[["id","smiles","md_ligand_rmsd","md_stability","md_status"]].to_csv(a.output, index=False)
        print(df["md_status"].iloc[0]); return
    reexec(C.get("python"))
    from common import read_any_table
    df=read_any_table(a.input)
    n=int(C.get("max_short", 20))
    if "docking_score" in df.columns:
        df=df.assign(_s=pd.to_numeric(df["docking_score"], errors="coerce")).sort_values("_s").head(n)
    elif "ml_probability" in df.columns:
        df=df.nlargest(n, "ml_probability")
    else:
        df=df.head(n)
    os.makedirs(a.workdir, exist_ok=True)
    rows=[]
    here=os.path.dirname(os.path.abspath(__file__))
    for r in df.itertuples(index=False):
        tag=str(r.id).replace("/","_")
        d=os.path.join(a.workdir, tag); os.makedirs(d, exist_ok=True)
        sdf=os.path.join(d, "lig.sdf"); prefix=os.path.join(d, "md")
        pose=getattr(r, "pose_sdf", None)
        if isinstance(pose, str) and os.path.isfile(pose):
            shutil.copy(pose, sdf)
            embed_rc=0
        else:
            q=subprocess.run([sys.executable, os.path.join(here, "smiles_to_sdf.py"), "--smiles", str(r.smiles), "--id", tag, "--output", sdf])
            embed_rc=q.returncode
        if embed_rc:
            rows.append({"id": r.id, "smiles": r.smiles, "md_ligand_rmsd": float("nan"), "md_stability": float("nan"), "md_status": "embed_fail"})
            continue
        q=subprocess.run([sys.executable, os.path.join(here, "md_openmm.py"), "--protein", a.protein, "--ligand", sdf, "--ns", str(C.get("short_ns", 100.0)), "--output-prefix", prefix])
        final=prefix+"_final.pdb"
        if q.returncode or not os.path.isfile(final):
            rows.append({"id": r.id, "smiles": r.smiles, "md_ligand_rmsd": float("nan"), "md_stability": float("nan"), "md_status": "md_fail"})
            continue
        rmsd=ligand_rmsd(sdf, final, a.protein)
        stab=1.0/(1.0+rmsd) if np.isfinite(rmsd) else float("nan")
        rows.append({"id": r.id, "smiles": r.smiles, "md_ligand_rmsd": rmsd, "md_stability": stab, "md_status": "ok"})
    pd.DataFrame(rows).to_csv(a.output, index=False)
    print(f"md_finished {sum(1 for r in rows if r['md_status']=='ok')}/{len(rows)}")

if __name__=="__main__":
    main()
