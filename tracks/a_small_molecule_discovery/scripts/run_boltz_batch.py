#!/usr/bin/env python
"""Fold the clustered short list with Boltz-2 and record affinity plus confidence."""
import argparse, json, os, subprocess
from pathlib import Path
import pandas as pd
import yaml
from common import load_yaml, read_any_table

AA = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

def pdb_sequence(path):
    seq, last = [], None
    for line in open(path):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        key = (line[21], line[22:26])
        if key == last:
            continue
        last = key
        seq.append(AA.get(line[17:20].strip(), "X"))
    return "".join(seq)

def safe_name(raw, used):
    s = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(raw)).strip("_") or "cmpd"
    s = s[:80]
    base = s
    n = 2
    while s in used:
        s = f"{base}_{n}"
        n += 1
    used.add(s)
    return s

def write_skip(df, path, reason):
    out = df.copy()
    out["boltz_affinity_pred_value"] = float("nan")
    out["boltz_affinity_probability"] = float("nan")
    out["boltz_confidence_score"] = float("nan")
    out["boltz_status"] = reason
    out[["id", "smiles", "boltz_affinity_pred_value", "boltz_affinity_probability", "boltz_confidence_score", "boltz_status"]].to_csv(path, index=False)
    print(reason)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--protein-pdb")
    ap.add_argument("--protein-sequence")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()
    C = load_yaml(a.config)["boltz"]
    seq = a.protein_sequence or (pdb_sequence(a.protein_pdb) if a.protein_pdb else "")
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    df = read_any_table(a.input)
    n = int(C.get("max_input", 200))
    if len(df) > n:
        if "docking_score" in df.columns:
            df = df.assign(_s=pd.to_numeric(df["docking_score"], errors="coerce")).sort_values("_s").head(n)
        elif "ml_probability" in df.columns:
            df = df.nlargest(n, "ml_probability")
        else:
            df = df.head(n)
    csv_path = out / "boltz_scores.csv"
    if not C.get("enabled", True):
        write_skip(df, csv_path, "skipped: boltz disabled")
        return
    exe = C.get("exe", "boltz")
    if not seq:
        write_skip(df, csv_path, "skipped: protein sequence")
        return
    if not (os.path.isfile(exe) or shutil_which(exe)):
        write_skip(df, csv_path, f"skipped: {exe} not found")
        return
    input_dir = out / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    used = set()
    rows = []
    for r in df.itertuples(index=False):
        job = safe_name(r.id, used)
        doc = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": seq, "msa": "empty"}},
                {"ligand": {"id": "B", "smiles": str(r.smiles)}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
        (input_dir / f"{job}.yaml").write_text(yaml.safe_dump(doc, sort_keys=False))
        rows.append({"id": r.id, "job": job, "smiles": r.smiles})
    cmd = [
        exe, "predict", str(input_dir),
        "--out_dir", str(out / "pred"),
        "--model", str(C.get("model", "boltz2")),
        "--accelerator", str(C.get("accelerator", "gpu")),
        "--devices", str(int(C.get("devices", 1))),
        "--recycling_steps", str(int(C.get("recycling_steps", 3))),
        "--sampling_steps", str(int(C.get("sampling_steps", 200))),
        "--diffusion_samples", str(int(C.get("diffusion_samples", 1))),
        "--diffusion_samples_affinity", str(int(C.get("diffusion_samples_affinity", 5))),
        "--seed", str(int(C.get("seed", 42))),
        "--output_format", "mmcif",
        "--override",
    ]
    if not C.get("use_kernels", False):
        cmd.append("--no_kernels")
    if C.get("use_msa_server"):
        cmd.append("--use_msa_server")
    if C.get("cache"):
        cmd += ["--cache", str(C["cache"])]
    q = subprocess.run(cmd)
    pred = out / "pred"
    scored = []
    for row in rows:
        job = row["job"]
        aff = next(pred.rglob(f"affinity_{job}.json"), None) if pred.exists() else None
        conf = next(pred.rglob(f"confidence_{job}_model_0.json"), None) if pred.exists() else None
        aval = apv = cscore = float("nan")
        status = "ok"
        if aff is None and conf is None:
            status = "boltz_fail"
        if aff is not None:
            data = json.loads(aff.read_text())
            aval = data.get("affinity_pred_value", float("nan"))
            apv = data.get("affinity_probability_binary", float("nan"))
        if conf is not None:
            cscore = json.loads(conf.read_text()).get("confidence_score", float("nan"))
        scored.append({
            "id": row["id"], "smiles": row["smiles"],
            "boltz_affinity_pred_value": aval,
            "boltz_affinity_probability": apv,
            "boltz_confidence_score": cscore,
            "boltz_status": status,
        })
    pd.DataFrame(scored).to_csv(csv_path, index=False)
    ok = sum(1 for r in scored if r["boltz_status"] == "ok")
    print(f"boltz_ok {ok}/{len(scored)}")
    if q.returncode and ok == 0:
        raise SystemExit(q.returncode)

def shutil_which(exe):
    import shutil
    return shutil.which(exe)

if __name__ == "__main__":
    main()
