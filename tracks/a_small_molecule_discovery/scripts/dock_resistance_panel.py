#!/usr/bin/env python
"""Dock the clustered shortlist across the KRAS resistance panel.

The billion-molecule screen stays on G12D. This step only re-docks a shortlist
into the G12D GDP state, the G12D GTP-like state, wild-type KRAS, and the
resistance structures. Scores are Vina rank heuristics.
"""
import argparse, json, os, shutil
from pathlib import Path
import numpy as np
import pandas as pd
from common import load_yaml, read_any_table
from dock_vina_batch import _dock_one
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor


def pct_good(series, lower_better=True):
    rank = series.rank(method="average", ascending=lower_better, pct=True)
    return 1 - rank + 1 / max(len(series), 1)


def write_empty(df, path, reason, long_path=""):
    out = df[["id", "smiles"]].copy() if {"id", "smiles"} <= set(df.columns) else df.copy()
    for col in [
        "G12D_mean_affinity", "WT_affinity", "WT_selectivity_margin",
        "resistance_exp_mean", "resistance_all_mean", "resistance_mean_penalty",
        "resistance_worst_penalty", "composite_rank_score",
        "pass_selectivity", "pass_resistance", "resistance_status",
    ]:
        out[col] = np.nan if col not in {"pass_selectivity", "pass_resistance", "resistance_status"} else False
    out["resistance_status"] = reason
    out.to_csv(path, index=False)
    if long_path:
        pd.DataFrame(columns=["id", "smiles", "pdb_id", "docking_score", "dock_status"]).to_csv(long_path, index=False)
    print(reason)


def rank_panel(long_df, manifest, ranking):
    receptors = manifest["receptors"]
    by_id = {r["pdb_id"]: r for r in receptors}
    primary = [ranking["primary_gdp"], ranking["primary_active"]]
    wt = ranking["wt"]
    exp_res = [r["pdb_id"] for r in receptors if r["role"] == "resistance" and r["source_type"] == "experimental"]
    all_res = [r["pdb_id"] for r in receptors if r["role"] == "resistance"]
    wide = long_df.pivot(index="id", columns="pdb_id", values="docking_score")
    smiles = long_df.drop_duplicates("id").set_index("id")["smiles"]
    rows = []
    for lig, row in wide.iterrows():
        prim_vals = np.array([row.get(x, np.nan) for x in primary], float)
        prim = np.nanmean(prim_vals) if np.isfinite(prim_vals).any() else np.nan
        wt_e = float(row.get(wt, np.nan))
        exp = np.array([row.get(x, np.nan) for x in exp_res], float)
        allr = np.array([row.get(x, np.nan) for x in all_res], float)
        exp_mean = np.nanmean(exp) if np.isfinite(exp).any() else np.nan
        all_mean = np.nanmean(allr) if np.isfinite(allr).any() else np.nan
        worst = np.nanmax(allr) if np.isfinite(allr).any() else np.nan
        select = (wt_e - prim) if np.isfinite(wt_e) and np.isfinite(prim) else np.nan
        penalty = (all_mean - prim) if np.isfinite(all_mean) and np.isfinite(prim) else np.nan
        worst_pen = (worst - prim) if np.isfinite(worst) and np.isfinite(prim) else np.nan
        item = {
            "id": lig,
            "smiles": smiles.get(lig, ""),
            "G12D_mean_affinity": prim,
            "WT_affinity": wt_e,
            "WT_selectivity_margin": select,
            "resistance_exp_mean": exp_mean,
            "resistance_all_mean": all_mean,
            "resistance_mean_penalty": penalty,
            "resistance_worst_penalty": worst_pen,
        }
        for pid in wide.columns:
            item[f"vina_{pid}"] = row.get(pid, np.nan)
            meta = by_id.get(pid, {})
            item[f"role_{pid}"] = meta.get("role", "")
        rows.append(item)
    ranked = pd.DataFrame(rows)
    if ranked.empty:
        return ranked
    ranked["s_affinity"] = pct_good(ranked["G12D_mean_affinity"], True)
    ranked["s_selectivity"] = pct_good(-ranked["WT_selectivity_margin"], True)
    ranked["s_robustness"] = pct_good(ranked["resistance_worst_penalty"], True)
    w = ranking["weights"]
    ranked["composite_rank_score"] = (
        float(w["g12d_affinity"]) * ranked["s_affinity"]
        + float(w["wt_selectivity"]) * ranked["s_selectivity"]
        + float(w["resistance_robustness"]) * ranked["s_robustness"]
    )
    ranked["pass_selectivity"] = ranked["WT_selectivity_margin"] >= float(ranking["min_selectivity_margin_kcal"])
    ranked["pass_resistance"] = ranked["resistance_worst_penalty"] <= float(ranking["max_resistance_penalty_kcal"])
    ranked["resistance_status"] = "ok"
    ranked = ranked.sort_values(
        ["pass_selectivity", "pass_resistance", "composite_rank_score"],
        ascending=[False, False, False],
    )
    return ranked


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--long", default="")
    a = ap.parse_args()
    cfg = load_yaml(a.config)
    rc = cfg.get("resistance", {})
    df = read_any_table(a.input)
    if not rc.get("enabled", True):
        write_empty(df, a.output, "skipped: resistance panel disabled", a.long)
        return
    manifest_path = Path(a.panel) / "panel_manifest.json"
    if not manifest_path.exists():
        write_empty(df, a.output, "skipped: panel manifest missing", a.long)
        return
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("disabled") or not manifest.get("receptors"):
        write_empty(df, a.output, "skipped: resistance panel disabled", a.long)
        return
    if "id" not in df.columns or "smiles" not in df.columns:
        raise SystemExit("Resistance docking needs id and smiles columns")
    n = int(rc.get("max_input", 200))
    if len(df) > n:
        if "docking_score" in df.columns:
            df = df.nsmallest(n, "docking_score")
        elif "ml_probability" in df.columns:
            df = df.nlargest(n, "ml_probability")
        else:
            df = df.head(n)
    dc = cfg["docking"]
    vina = dc.get("vina_exe", "vina")
    meeko = dc.get("meeko_ligand_exe", "mk_prepare_ligand.py")
    for exe in (vina, meeko):
        if shutil.which(exe) is None:
            write_empty(df, a.output, f"skipped: {exe} not found", a.long)
            return
    exh = int(rc.get("exhaustiveness", dc.get("redock_exhaustiveness", 8)))
    modes = int(dc.get("num_modes", 1))
    optimize = bool(dc.get("embed_optimize", False))
    workers = max(1, min(int(dc.get("workers") or 16), 16))
    center = tuple(manifest["center"])
    size = tuple(manifest["size"])
    tasks = []
    for receptor in manifest["receptors"]:
        pdbqt = str((Path(a.panel) / receptor["pdbqt"]).resolve())
        if not os.path.isfile(pdbqt):
            print("missing receptor", pdbqt)
            continue
        for rec in df.itertuples(index=False):
            tasks.append((rec.smiles, rec.id, {"pdb_id": receptor["pdb_id"]}, pdbqt, center, size,
                          exh, modes, vina, meeko, optimize))
    if not tasks:
        write_empty(df, a.output, "skipped: no prepared receptors", a.long)
        return
    print(f"resistance docks {len(tasks)} ligands x receptors, exhaustiveness {exh}, workers {workers}")
    if workers == 1:
        rows = [_dock_one(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            rows = list(pool.map(_dock_one, tasks, chunksize=4))
    long_df = pd.DataFrame(rows)
    if a.long:
        long_df.to_csv(a.long, index=False)
    ok = int((long_df["dock_status"] == "ok").sum()) if "dock_status" in long_df else 0
    print(f"resistance docking ok {ok}/{len(long_df)}")
    ranked = rank_panel(long_df, manifest, rc)
    ranked.to_csv(a.output, index=False)
    if not ranked.empty:
        both = int((ranked["pass_selectivity"] & ranked["pass_resistance"]).sum())
        print(f"shortlist {len(ranked)} pass selectivity and resistance {both}")


if __name__ == "__main__":
    main()
