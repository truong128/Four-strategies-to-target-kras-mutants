#!/usr/bin/env python
"""Boltz-2 affinity of FDA-approved small molecules against KRAS G12D only."""
import json
import os
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import load_yaml

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
# Drug libraries are user supplied and are not part of this repository.
#   FDA_DRUG_PKL: pandas pickle with columns drug_id, drug_name, smile
#   FDA_DRUG_CSV: optional CSV with columns name, smiles
PKL = Path(os.environ.get("FDA_DRUG_PKL", "data/fda_drugs/approved_small_molecule_drugs.pkl"))
EXTRA = Path(os.environ.get("FDA_DRUG_CSV", "data/fda_drugs/fda_drugs_smiles.csv"))
WT_PDB = ROOT / "results" / "allele_dock" / "receptors" / "6MBU.pdb"


def chain_a(path):
    out, last = {}, None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != "A" or line[12:16].strip() != "CA":
            continue
        resi = int(line[22:26])
        if resi == last:
            continue
        last = resi
        out[resi] = AA3[line[17:20].strip()]
    seq = "".join(out[i] for i in range(1, 169))
    if seq[11] != "G":
        raise SystemExit("6MBU residue 12 is not glycine")
    return seq[:11] + "D" + seq[12:]


def parent_smiles(raw):
    mol = Chem.MolFromSmiles(str(raw))
    if mol is None:
        return None
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    organic = []
    for frag in frags:
        if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
            organic.append(frag)
    if not organic:
        return None
    frag = max(organic, key=lambda m: m.GetNumHeavyAtoms())
    if frag.GetNumHeavyAtoms() > 80 or Descriptors.MolWt(frag) > 1200:
        return None
    return Chem.MolToSmiles(frag, isomericSmiles=True)


def safe(raw):
    s = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(raw)).strip("_")
    return (s or "drug")[:80]


def load_drugs():
    if not PKL.exists():
        raise SystemExit(f"Set FDA_DRUG_PKL to your approved-drug table (missing: {PKL})")
    frame = pickle.load(open(PKL, "rb"))
    rows = []
    seen = set()
    for rec in frame.itertuples(index=False):
        smi = parent_smiles(rec.smile)
        if not smi or smi in seen:
            continue
        seen.add(smi)
        rows.append({"drug_id": str(rec.drug_id), "drug_name": str(rec.drug_name), "smiles": smi, "source": "approved_small_molecule_drugs_review"})
    if EXTRA.exists():
        extra = pd.read_csv(EXTRA)
        for rec in extra.itertuples(index=False):
            smi = parent_smiles(rec.smiles)
            if not smi or smi in seen:
                continue
            seen.add(smi)
            rows.append({"drug_id": safe(rec.name), "drug_name": str(rec.name), "smiles": smi, "source": "FDA_Drug_Data_SMILES_pipeline_ready"})
    return rows


def main():
    out = ROOT / "results" / "g12d_fda_boltz_msa"
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    msa_csv = ROOT / "results" / "kras_msa" / "G12D.csv"
    if not msa_csv.exists():
        raise SystemExit(f"Missing MSA {msa_csv}. Run scripts/fetch_kras_msa.py first.")
    seq = chain_a(WT_PDB)
    drugs = load_drugs()
    rows = []
    used = set()
    for drug in drugs:
        job = safe(drug["drug_id"])
        base = job
        n = 2
        while job in used:
            job = f"{base}_{n}"
            n += 1
        used.add(job)
        doc = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": seq, "msa": str(msa_csv.resolve())}},
                {"ligand": {"id": "B", "smiles": drug["smiles"]}},
            ],
            "properties": [{"affinity": {"binder": "B"}}],
        }
        (inputs / f"{job}.yaml").write_text(yaml.safe_dump(doc, sort_keys=False))
        rows.append({**drug, "job": job, "allele": "G12D"})
    pd.DataFrame(rows).to_csv(out / "jobs.csv", index=False)
    print(f"wrote {len(rows)} G12D Boltz-2 inputs", flush=True)
    cfg = load_yaml(ROOT / "config.yaml")["boltz"]
    exe = cfg.get("exe", "boltz")
    cmd = [
        exe, "predict", str(inputs),
        "--out_dir", str(out / "pred"),
        "--model", str(cfg.get("model", "boltz2")),
        "--accelerator", str(cfg.get("accelerator", "gpu")),
        "--devices", str(int(cfg.get("devices", 1))),
        "--recycling_steps", str(int(cfg.get("recycling_steps", 3))),
        "--sampling_steps", str(int(cfg.get("sampling_steps", 200))),
        "--diffusion_samples", str(int(cfg.get("diffusion_samples", 1))),
        "--diffusion_samples_affinity", str(int(cfg.get("diffusion_samples_affinity", 5))),
        "--seed", str(int(cfg.get("seed", 42))),
        "--output_format", "mmcif",
        "--override",
        "--no_kernels",
    ]
    print("exec", " ".join(cmd), flush=True)
    q = subprocess.run(cmd)
    pred = out / "pred"
    scored = []
    for row in rows:
        job = row["job"]
        aff = next(pred.rglob(f"affinity_{job}.json"), None) if pred.exists() else None
        conf = next(pred.rglob(f"confidence_{job}_model_0.json"), None) if pred.exists() else None
        aval = apv = cscore = float("nan")
        status = "ok" if aff is not None else "boltz_fail"
        if aff is not None:
            data = json.loads(aff.read_text())
            aval = data.get("affinity_pred_value", float("nan"))
            apv = data.get("affinity_probability_binary", float("nan"))
        if conf is not None:
            cscore = json.loads(conf.read_text()).get("confidence_score", float("nan"))
        scored.append({
            **row,
            "boltz_affinity_pred_value": aval,
            "boltz_affinity_probability": apv,
            "boltz_confidence_score": cscore,
            "boltz_status": status,
        })
    table = pd.DataFrame(scored).sort_values("boltz_affinity_probability", ascending=False)
    path = out / "g12d_fda_boltz.csv"
    table.to_csv(path, index=False)
    ok = int((table["boltz_status"] == "ok").sum())
    print(f"boltz_ok {ok}/{len(table)}", flush=True)
    print("wrote", path, flush=True)
    if q.returncode and ok == 0:
        raise SystemExit(q.returncode)


if __name__ == "__main__":
    main()
