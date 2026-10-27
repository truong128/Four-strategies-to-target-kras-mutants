#!/usr/bin/env python
"""RDKit ease score plus RetroChimera / Retro* route search.

The external search uses the local USPTO-50K RetroChimera checkpoint and the
Retro* purchasable-building-block inventory. A solved route raises the score.
Unsolved molecules keep the RDKit ease score and a route score of zero.
"""
import argparse, math, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")


def project_file(path, config_path):
    p = Path(path)
    if p.is_file():
        return p
    for base in (Path(config_path).resolve().parent, Path(__file__).resolve().parents[1]):
        cand = base / path
        if cand.is_file():
            return cand
    return p


def ease_score(mol):
    comp = (
        0.015 * Descriptors.MolWt(mol)
        + 0.35 * Lipinski.RingCount(mol)
        + 0.2 * Lipinski.NumRotatableBonds(mol)
        + 0.45 * len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        + 0.6 * rdMolDescriptors.CalcNumSpiroAtoms(mol)
        + 0.6 * rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    )
    return float(np.exp(-comp / 10.0)), float(comp)


def shortlist(df, n):
    if len(df) <= n:
        return df
    if "docking_score" in df.columns:
        return df.assign(_s=pd.to_numeric(df["docking_score"], errors="coerce")).sort_values("_s").head(n)
    if "ml_probability" in df.columns:
        return df.nlargest(n, "ml_probability")
    return df.head(n)


def load_inventory(path):
    from syntheseus.search.mol_inventory import SmilesListInventory
    smiles = [ln.strip() for ln in open(path) if ln.strip()]
    return SmilesListInventory(smiles, canonicalize=False)


def route_steps(graph):
    from syntheseus.search.analysis.route_extraction import iter_routes_time_order
    from syntheseus.search.graph.and_or import AndNode
    route = next(iter_routes_time_order(graph, max_routes=1), None)
    if route is None:
        return float("nan"), ""
    reactions = []
    for node in route:
        if isinstance(node, AndNode):
            reactions.append(str(node.reaction))
    text = " >> ".join(reactions)
    return float(len(reactions)), text[:2000]


def search_batch(df, cfg, config_path):
    import torch
    from retrochimera import RetroChimeraModel
    from syntheseus import Molecule
    from syntheseus.search.algorithms.best_first.retro_star import RetroStarSearch
    from syntheseus.search.analysis.solution_time import get_first_solution_time
    from syntheseus.search.node_evaluation.common import ConstantNodeEvaluator, ReactionModelLogProbCost

    torch.set_float32_matmul_precision("high")
    model_dir = project_file(cfg["model_dir"], config_path)
    inventory_path = project_file(cfg["inventory"], config_path)
    if not model_dir.is_dir():
        raise SystemExit(f"RetroChimera checkpoint directory is missing: {model_dir}")
    if not inventory_path.is_file():
        raise SystemExit(f"Purchasable inventory is missing: {inventory_path}")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"loading RetroChimera from {model_dir} on {device}")
    model = RetroChimeraModel(
        model_dir=str(model_dir),
        device=device,
        default_num_results=int(cfg.get("num_results", 10)),
        use_cache=True,
        smiles_transformer={"canonicalization_processes": 1},
    )
    print(f"loading inventory {inventory_path}")
    inventory = load_inventory(inventory_path)
    print(f"inventory size {len(inventory)}")
    RDLogger.DisableLog("rdApp.*")
    alg = RetroStarSearch(
        reaction_model=model,
        mol_inventory=inventory,
        and_node_cost_fn=ReactionModelLogProbCost(),
        value_function=ConstantNodeEvaluator(0.0),
        time_limit_s=float(cfg.get("time_limit_s", 60)),
        max_expansion_depth=int(cfg.get("max_depth", 6)),
        limit_iterations=int(cfg.get("limit_iterations", 100000)),
        limit_reaction_model_calls=int(cfg.get("limit_model_calls", 100000)),
        stop_on_first_solution=bool(cfg.get("stop_on_first_solution", True)),
        prevent_repeat_mol_in_trees=True,
    )
    rows = []
    for i, rec in enumerate(df.itertuples(index=False), 1):
        mol = Chem.MolFromSmiles(str(rec.smiles))
        ease, comp = ease_score(mol) if mol is not None else (float("nan"), float("nan"))
        solved, steps, route, status, calls = False, float("nan"), "", "fail", 0
        if mol is None:
            status = "invalid_smiles"
        else:
            try:
                alg.reset()
                graph, _ = alg.run_from_mol(Molecule(str(rec.smiles)))
                for node in graph.nodes():
                    node.data["analysis_time"] = node.data.get("num_calls_rxn_model", 0)
                solved = math.isfinite(get_first_solution_time(graph))
                calls = int(alg.reaction_model.num_calls())
                if solved:
                    steps, route = route_steps(graph)
                    status = "solved"
                else:
                    status = "unsolved"
            except Exception as exc:
                status = f"fail:{type(exc).__name__}"
                print(f"retro fail {rec.id}: {exc}", file=sys.stderr)
        score = (1.0 / (1.0 + steps)) if solved and np.isfinite(steps) else (1.0 if solved else 0.0)
        rows.append({
            "id": rec.id,
            "smiles": rec.smiles,
            "synthetic_accessibility": ease,
            "synthetic_complexity_raw": comp,
            "retrosynthesis_solved": bool(solved),
            "retrosynthesis_steps": steps,
            "retrosynthesis_score": score,
            "retrosynthesis_model_calls": calls,
            "retrosynthesis_route": route,
            "retrosynthesis_status": status,
        })
        print(f"[{i}/{len(df)}] {rec.id} {status} steps={steps} score={score:.3f}")
    return pd.DataFrame(rows)


def rdkit_only(df, status):
    rows = []
    for rec in df.itertuples(index=False):
        mol = Chem.MolFromSmiles(str(rec.smiles))
        if mol is None:
            continue
        ease, comp = ease_score(mol)
        rows.append({
            "id": rec.id,
            "smiles": rec.smiles,
            "synthetic_accessibility": ease,
            "synthetic_complexity_raw": comp,
            "retrosynthesis_solved": False,
            "retrosynthesis_steps": float("nan"),
            "retrosynthesis_score": 0.0,
            "retrosynthesis_model_calls": 0,
            "retrosynthesis_route": "",
            "retrosynthesis_status": status,
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    import yaml
    cfg = yaml.safe_load(open(a.config))["retrosynthesis"]
    py = cfg.get("python", "")
    if cfg.get("enabled") and py and os.path.isfile(py) and os.path.realpath(py) != os.path.realpath(sys.executable):
        os.execv(py, [py, os.path.abspath(__file__), *sys.argv[1:]])
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from common import read_any_table
    df = read_any_table(a.input)
    if "id" not in df.columns or "smiles" not in df.columns:
        raise SystemExit("Synthesis scoring needs id and smiles columns")
    df = shortlist(df, int(cfg.get("max_input", 50)))
    if not cfg.get("enabled", False):
        out = rdkit_only(df, "skipped: retrosynthesis disabled")
    else:
        out = search_batch(df, cfg, a.config)
    out.to_csv(a.output, index=False)
    solved = int(out["retrosynthesis_solved"].sum()) if len(out) else 0
    print(f"synthesis {solved}/{len(out)} routes solved")


if __name__ == "__main__":
    main()
