#!/usr/bin/env python
"""Boltz-2 affinity of the G12D-passing drugs on the other 17 KRAS alleles.

Input: results/g12d_fda_boltz_msa/g12d_passed_gt_0.5.csv (the rows of the G12D Boltz table that
you keep, for example boltz_affinity_probability > 0.5; needs drug_id, drug_name, smiles).
Uses the MSA files saved by fetch_kras_msa.py in results/kras_msa. G12D is not rerun.
Jobs = number of kept drugs x 17 alleles.
"""
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from boltz_allele_panel import apply_mutation, chain_a_map
from common import load_yaml

PASSED = ROOT / "results" / "g12d_fda_boltz_msa" / "g12d_passed_gt_0.5.csv"
PANEL = ROOT / "config" / "kras_allele_panel.tsv"
WT_PDB = ROOT / "results" / "allele_dock" / "receptors" / "6MBU.pdb"
OUT = ROOT / "results" / "fda_passed_17alleles_boltz"


def sequences():
    wt_map = chain_a_map(WT_PDB)
    wt = "".join(wt_map[i] for i in range(1, 169))
    if len(wt) != 168 or wt[11] != "G":
        raise SystemExit("6MBU residues 1-168 are not the expected KRAS sequence")
    panel = pd.read_csv(PANEL, sep="\t")
    out = []
    for rec in panel.itertuples(index=False):
        if rec.allele == "G12D":
            continue
        mutation = str(rec.mutation) if pd.notna(rec.mutation) and str(rec.mutation) else ""
        seq = apply_mutation(wt, rec.allele, mutation)
        msa = ROOT / "results" / "kras_msa" / f"{rec.allele}.csv"
        if not msa.exists():
            raise SystemExit(f"Missing MSA {msa}. Run scripts/fetch_kras_msa.py first.")
        out.append((rec.allele, rec.pdb_id, rec.source_type, seq, msa.resolve()))
    if len(out) != 17:
        raise SystemExit(f"Expected 17 non-G12D alleles, found {len(out)}")
    return out


def main():
    drugs = pd.read_csv(PASSED)
    if drugs["drug_id"].duplicated().any():
        raise SystemExit(f"Duplicate drug_id in {PASSED}")
    alleles = sequences()
    inputs = OUT / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    rows = []
    for allele, pdb_id, source_type, seq, msa in alleles:
        for drug in drugs.itertuples(index=False):
            job = f"{allele}__{drug.drug_id}"
            doc = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq, "msa": str(msa)}},
                    {"ligand": {"id": "B", "smiles": str(drug.smiles)}},
                ],
                "properties": [{"affinity": {"binder": "B"}}],
            }
            (inputs / f"{job}.yaml").write_text(yaml.safe_dump(doc, sort_keys=False))
            rows.append({
                "job": job,
                "drug_id": drug.drug_id,
                "drug_name": drug.drug_name,
                "smiles": drug.smiles,
                "allele": allele,
                "pdb_id": pdb_id,
                "source_type": source_type,
                "G12D_affinity_probability": float(drug.boltz_affinity_probability),
                "G12D_affinity_pred_value": float(drug.boltz_affinity_pred_value),
            })
    pd.DataFrame(rows).to_csv(OUT / "jobs.csv", index=False)
    print(f"wrote {len(rows)} Boltz-2 inputs to {inputs}", flush=True)

    cfg = load_yaml(ROOT / "config.yaml")["boltz"]
    exe = cfg.get("exe", "boltz")
    cmd = [
        exe, "predict", str(inputs),
        "--out_dir", str(OUT / "pred"),
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
    import subprocess
    q = subprocess.run(cmd)

    pred = OUT / "pred"
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
            "passed_gt_0.5": bool(apv > 0.5) if aff is not None else False,
        })
    long_df = pd.DataFrame(scored)
    long_df["delta_affinity_probability_vs_G12D"] = (
        long_df["boltz_affinity_probability"] - long_df["G12D_affinity_probability"]
    )
    long_path = OUT / "fda_17alleles_long.csv"
    long_df.to_csv(long_path, index=False)

    wide = long_df.pivot(index="drug_id", columns="allele", values="boltz_affinity_probability")
    g12d = drugs.set_index("drug_id")["boltz_affinity_probability"].rename("G12D")
    wide = wide.join(g12d)
    names = drugs.set_index("drug_id")["drug_name"]
    wide.insert(0, "drug_name", names)
    allele_cols = ["G12D"] + [a for a, *_ in alleles]
    wide = wide[["drug_name"] + allele_cols]
    wide_path = OUT / "fda_18alleles_wide.csv"
    wide.to_csv(wide_path)

    passed_rows = []
    for drug_id, block in long_df.groupby("drug_id", sort=False):
        ok_block = block[block["boltz_status"] == "ok"]
        hit = ok_block[ok_block["boltz_affinity_probability"] > 0.5]
        g12d_p = float(block["G12D_affinity_probability"].iloc[0])
        alleles_passed = ["G12D"] + hit["allele"].tolist()
        passed_rows.append({
            "drug_id": drug_id,
            "drug_name": block["drug_name"].iloc[0],
            "G12D_affinity_probability": g12d_p,
            "n_of_17_passed": int(len(hit)),
            "n_of_18_passed": int(len(hit) + (g12d_p > 0.5)),
            "alleles_passed": ",".join(alleles_passed),
        })
    summary = pd.DataFrame(passed_rows).sort_values(
        ["n_of_17_passed", "G12D_affinity_probability"], ascending=False
    )
    summary_path = OUT / "fda_passed_across_alleles.csv"
    summary.to_csv(summary_path, index=False)

    ok = int((long_df["boltz_status"] == "ok").sum())
    n_drugs = int((summary["n_of_17_passed"] > 0).sum())
    print(f"boltz_ok {ok}/{len(long_df)}", flush=True)
    print(f"drugs passing >0.5 on at least one of the 17 alleles: {n_drugs}/{len(summary)}", flush=True)
    print("wrote", long_path, flush=True)
    print("wrote", wide_path, flush=True)
    print("wrote", summary_path, flush=True)
    if q.returncode and ok == 0:
        raise SystemExit(q.returncode)


if __name__ == "__main__":
    main()
