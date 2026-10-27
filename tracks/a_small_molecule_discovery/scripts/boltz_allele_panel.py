#!/usr/bin/env python
"""Boltz-2 affinity of the final leads on 18 KRAS alleles.

Each allele is one substitution on wild-type 6MBU residues 1-168. Crystal
construct differences (extra N-terminal residues, C118S, missing loops) are
not copied into the sequence, so scores differ by the allele only.
Jobs = number of leads x 18 alleles.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from common import load_yaml, read_any_table

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
def chain_a_map(path):
    out, last = {}, None
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != "A" or line[12:16].strip() != "CA":
            continue
        resi = int(line[22:26])
        if resi == last:
            continue
        last = resi
        out[resi] = AA3[line[17:20].strip()]
    return out


def mutation_of(allele, mutation_field):
    if allele == "WT":
        return None
    if mutation_field:
        wt3, resi, mut3 = mutation_field.split("-")
        return int(resi), AA3[wt3], AA3[mut3]
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", allele)
    if not m:
        raise SystemExit(f"Cannot parse allele {allele}")
    return int(m.group(2)), m.group(1), m.group(3)


def apply_mutation(seq, allele, mutation_field):
    spec = mutation_of(allele, mutation_field)
    if spec is None:
        return seq
    resi, wt, mut = spec
    idx = resi - 1
    if seq[idx] != wt:
        raise SystemExit(f"{allele}: residue {resi} is {seq[idx]}, expected {wt}")
    return seq[:idx] + mut + seq[idx + 1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", default=str(ROOT / "config" / "kras_allele_panel.tsv"))
    ap.add_argument("--leads", default=str(ROOT / "results" / "cluster" / "diverse_hits.csv"))
    ap.add_argument("--wt-pdb", default=str(ROOT / "results" / "allele_dock" / "receptors" / "6MBU.pdb"))
    ap.add_argument("--config", default=str(ROOT / "config.yaml"))
    ap.add_argument("--output-dir", default=str(ROOT / "results" / "allele_boltz"))
    a = ap.parse_args()
    cfg = load_yaml(a.config)["boltz"]
    out = Path(a.output_dir)
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    wt_map = chain_a_map(a.wt_pdb)
    wt = "".join(wt_map[i] for i in range(1, 169))
    if len(wt) != 168 or wt[11] != "G":
        raise SystemExit("6MBU residues 1-168 are not the expected KRAS sequence")
    panel = pd.read_csv(a.panel, sep="\t")
    leads = read_any_table(a.leads)
    rows = []
    for rec in panel.itertuples(index=False):
        seq = apply_mutation(wt, rec.allele, str(rec.mutation) if pd.notna(rec.mutation) and str(rec.mutation) else "")
        for lead in leads.itertuples(index=False):
            job = f"{rec.allele}__{lead.id}"
            msa_csv = ROOT / "results" / "kras_msa" / f"{rec.allele}.csv"
            if not msa_csv.exists():
                raise SystemExit(f"Missing MSA {msa_csv}. Run scripts/fetch_kras_msa.py first.")
            doc = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq, "msa": str(msa_csv.resolve())}},
                    {"ligand": {"id": "B", "smiles": str(lead.smiles)}},
                ],
                "properties": [{"affinity": {"binder": "B"}}],
            }
            (inputs / f"{job}.yaml").write_text(yaml.safe_dump(doc, sort_keys=False))
            rows.append({
                "job": job, "id": lead.id, "smiles": lead.smiles,
                "allele": rec.allele, "pdb_id": rec.pdb_id, "source_type": rec.source_type,
            })
    manifest = out / "jobs.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    print(f"wrote {len(rows)} Boltz-2 inputs to {inputs}", flush=True)
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
    ]
    if not cfg.get("use_kernels", False):
        cmd.append("--no_kernels")
    if cfg.get("cache"):
        cmd += ["--cache", str(cfg["cache"])]
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
    long_df = pd.DataFrame(scored)
    g12d = long_df.loc[long_df["allele"] == "G12D", ["id", "boltz_affinity_probability"]].rename(
        columns={"boltz_affinity_probability": "G12D_affinity_probability"}
    )
    long_df = long_df.merge(g12d, on="id", how="left")
    long_df["delta_affinity_probability_vs_G12D"] = (
        long_df["boltz_affinity_probability"] - long_df["G12D_affinity_probability"]
    )
    long_path = out / "allele_boltz_long.csv"
    long_df.to_csv(long_path, index=False)
    wide = long_df.pivot(index="id", columns="allele", values="boltz_affinity_probability")
    wide.to_csv(out / "allele_boltz_wide.csv")
    ok = int((long_df["boltz_status"] == "ok").sum())
    print(f"boltz_ok {ok}/{len(long_df)}", flush=True)
    print("wrote", long_path, flush=True)
    if q.returncode and ok == 0:
        raise SystemExit(q.returncode)


if __name__ == "__main__":
    main()
