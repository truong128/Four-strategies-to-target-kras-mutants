#!/usr/bin/env python
"""Absolute binding free energy (ABFE) with OpenFE, as an independent check on MM/GBSA.

MM/GBSA is a single end-point estimate. OpenFE's AbsoluteBindingProtocol is an
alchemical calculation (decoupling cycle, Boresch restraints, HREX, MBAR) that
needs far more sampling and gives a free energy with an uncertainty. It is meant
for the short list only, and it is reported next to MM/GBSA, not averaged with it.

Stages (all per compound, resumable):
  plan     build an OpenFE Transformation JSON from the receptor and the docked pose
  run      `openfe quickrun` on each JSON (GPU)
  collect  read the result JSONs into one table, optionally compare with MM/GBSA

Status values: DONE, SKIPPED, NOT_INSTALLED, NOT_READY, FAILED. Only DONE carries a number.

Requirements for a compound to be planned:
  * a 3D pose with hydrogens in the `pose_sdf` column (a docked or Boltz pose);
  * a receptor with hydrogens. If the PDB has none, PDBFixer adds them (pH 7.4)
    when it is installed in the OpenFE environment, otherwise the compound is NOT_READY;
  * every cofactor that sits in the pocket (for KRAS: GDP, Mg2+) must be supplied
    as an SDF in `openfe.cofactors_sdf`, or the calculation describes a different system.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
KCAL_PER_KJ = 1.0 / 4.184


def _read_table(path):
    import pandas as pd
    sep = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
    return pd.read_csv(path, sep=sep)


def _load_yaml(path):
    import yaml
    with open(path) as fh:
        return yaml.safe_load(fh)


def find_openfe_env_exe(name: str, configured: str | None = None) -> str | None:
    """Return the path to `name` in an OpenFE conda env, an explicit path, or PATH."""
    if configured:
        c = os.path.expanduser(configured)
        if os.path.isabs(c) and os.path.isfile(c):
            return c
        w = shutil.which(c)
        if w and "openfe" in w:
            return w
    for env in ("openfe", "openfe_env"):
        for root in ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "micromamba"):
            cand = Path.home() / root / "envs" / env / "bin" / name
            if cand.is_file():
                return str(cand)
    return shutil.which(name)


def _slug(x) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in str(x))


# --------------------------------------------------------------------------- plan
def plan_one(args) -> int:
    """Runs inside the OpenFE environment. Writes abfe.json for one ligand."""
    import openfe
    from openfe.protocols.openmm_afe import AbsoluteBindingProtocol
    from openff.units import unit as offunit
    from rdkit import Chem

    work = Path(args.dir)
    work.mkdir(parents=True, exist_ok=True)
    prot_path = Path(args.protein)
    lines = prot_path.read_text().splitlines()
    n_h = sum(1 for ln in lines if ln.startswith(("ATOM", "HETATM")) and ln[76:78].strip() in {"H", "D"})
    if n_h == 0:
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
        except ImportError:
            print("NOT_READY: receptor has no hydrogens and pdbfixer is not installed")
            return 3
        fx = PDBFixer(filename=str(prot_path))
        fx.findMissingResidues()
        fx.missingResidues = {}
        fx.findNonstandardResidues()
        fx.replaceNonstandardResidues()
        fx.removeHeterogens(False)
        fx.findMissingAtoms()
        fx.addMissingAtoms()
        fx.addMissingHydrogens(7.4)
        prot_path = work / "receptor_H.pdb"
        with open(prot_path, "w") as fh:
            PDBFile.writeFile(fx.topology, fx.positions, fh, keepIds=True)

    protein = openfe.ProteinComponent.from_pdb_file(str(prot_path))

    suppl = Chem.SDMolSupplier(args.ligand, removeHs=False)
    mol = suppl[0] if len(suppl) else None
    if mol is None:
        print("NOT_READY: ligand SDF could not be read")
        return 3
    if mol.GetNumConformers() == 0 or not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        print("NOT_READY: ligand pose needs 3D coordinates and explicit hydrogens")
        return 3
    mol.SetProp("_Name", str(args.id))
    ligand = openfe.SmallMoleculeComponent.from_rdkit(mol, name=str(args.id))

    cofactors = []
    for p in args.cofactor or []:
        for m in Chem.SDMolSupplier(p, removeHs=False):
            if m is not None:
                cofactors.append(openfe.SmallMoleculeComponent.from_rdkit(m))

    solvent = openfe.SolventComponent()
    comps_a = {"ligand": ligand, "protein": protein, "solvent": solvent}
    comps_b = {"protein": protein, "solvent": solvent}
    for i, c in enumerate(cofactors):
        comps_a[f"cofactor{i}"] = c
        comps_b[f"cofactor{i}"] = c
    state_a = openfe.ChemicalSystem(comps_a, name=f"{args.id}_bound")
    state_b = openfe.ChemicalSystem(comps_b, name=f"{args.id}_apo")

    s = AbsoluteBindingProtocol.default_settings()
    s.protocol_repeats = int(args.repeats)
    s.complex_simulation_settings.production_length = args.complex_ns * offunit.nanosecond
    s.solvent_simulation_settings.production_length = args.solvent_ns * offunit.nanosecond
    s.complex_simulation_settings.equilibration_length = min(1.0, args.complex_ns) * offunit.nanosecond
    s.solvent_simulation_settings.equilibration_length = min(1.0, args.solvent_ns) * offunit.nanosecond
    if args.preequil_ns > 0:  # shorten the plain-MD pre-equilibration of both legs (smoke tests)
        for leg in (s.complex_equil_simulation_settings, s.solvent_equil_simulation_settings):
            leg.equilibration_length_nvt = args.preequil_ns * offunit.nanosecond
            leg.equilibration_length = args.preequil_ns * offunit.nanosecond
            leg.production_length = args.preequil_ns * offunit.nanosecond
    if args.charge_method == "nagl":  # fast AM1-BCC-trained GNN charges, same model the MM/GBSA ligands use
        s.partial_charge_settings.partial_charge_method = "nagl"
        s.partial_charge_settings.nagl_model = "openff-gnn-am1bcc-1.0.0.pt"
    if args.platform:
        s.engine_settings.compute_platform = args.platform.lower()
    protocol = AbsoluteBindingProtocol(settings=s)
    tf = openfe.Transformation(state_a, state_b, protocol, mapping=None, name=f"abfe_{args.id}")
    tf.to_json(work / "abfe.json")
    print("DONE: planned")
    return 0


# --------------------------------------------------------------------------- collect
def _quantity_kcal(q):
    """Convert an OpenFE serialized pint quantity (or a float) to kcal/mol."""
    if q is None:
        return float("nan")
    if isinstance(q, (int, float)):
        return float(q)
    if isinstance(q, dict):
        mag = q.get("magnitude")
        unit = str(q.get("unit", "kilocalorie_per_mole"))
    else:
        mag = getattr(q, "m", None)
        unit = str(getattr(q, "u", "kilocalorie_per_mole"))
    if mag is None:
        return float("nan")
    mag = float(mag)
    if "kilojoule" in unit:
        return mag * KCAL_PER_KJ
    return mag


def read_result_json(path: Path):
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return float("nan"), float("nan"), "FAILED", "result file unreadable"
    est, unc = d.get("estimate"), d.get("uncertainty")
    if est is None:
        return float("nan"), float("nan"), "FAILED", "no estimate in result (simulation failed)"
    return _quantity_kcal(est), _quantity_kcal(unc), "DONE", ""


def compare_with_mmgbsa(df, mm_path, out_json):
    import pandas as pd
    from scipy.stats import kendalltau, spearmanr
    mm = _read_table(mm_path)
    if "mmgbsa" not in mm.columns:
        return None
    m = df.merge(mm[["id", "mmgbsa"]], on="id", how="inner")
    m = m[(m["openfe_status"] == "DONE") & pd.to_numeric(m["mmgbsa"], errors="coerce").notna()]
    out = {"n_paired": int(len(m)), "note": "both scores are more negative for stronger binding"}
    if len(m) >= 3:
        rho, p = spearmanr(m["openfe_abfe_dg_kcal"], m["mmgbsa"])
        tau, pt = kendalltau(m["openfe_abfe_dg_kcal"], m["mmgbsa"])
        out.update(spearman_rho=float(rho), spearman_p=float(p), kendall_tau=float(tau), kendall_p=float(pt))
    else:
        out["status"] = "NOT_READY: fewer than 3 compounds have both OpenFE and MM/GBSA values"
    Path(out_json).write_text(json.dumps(out, indent=2))
    return out


# --------------------------------------------------------------------------- driver
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input")
    ap.add_argument("--protein")
    ap.add_argument("--config")
    ap.add_argument("--output")
    ap.add_argument("--workdir")
    ap.add_argument("--mmgbsa", help="mmpbsa_scores.csv to compare against")
    ap.add_argument("--stage", choices=["all", "plan", "run", "collect"], default="all")
    # internal: executed inside the OpenFE env
    ap.add_argument("--plan-one", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--id", help=argparse.SUPPRESS)
    ap.add_argument("--ligand", help=argparse.SUPPRESS)
    ap.add_argument("--dir", help=argparse.SUPPRESS)
    ap.add_argument("--cofactor", action="append", help=argparse.SUPPRESS)
    ap.add_argument("--repeats", type=int, default=3, help=argparse.SUPPRESS)
    ap.add_argument("--complex-ns", type=float, default=5.0, help=argparse.SUPPRESS)
    ap.add_argument("--solvent-ns", type=float, default=10.0, help=argparse.SUPPRESS)
    ap.add_argument("--platform", default="", help=argparse.SUPPRESS)
    ap.add_argument("--charge-method", default="am1bcc", help=argparse.SUPPRESS)
    ap.add_argument("--preequil-ns", type=float, default=0.0, help=argparse.SUPPRESS)
    a = ap.parse_args()

    if a.plan_one:
        return plan_one(a)

    import pandas as pd
    for need in ("input", "protein", "config", "output", "workdir"):
        if getattr(a, need) is None:
            ap.error(f"--{need} is required")
    cfg = _load_yaml(a.config).get("openfe", {}) or {}
    df = _read_table(a.input)
    cols = ["id", "smiles"]

    def write(rows):
        out = pd.DataFrame(rows)
        for c in ("openfe_abfe_dg_kcal", "openfe_abfe_dg_sem"):
            if c not in out:
                out[c] = float("nan")
        out.to_csv(a.output, index=False)
        return out

    def all_status(status, note):
        rows = [{"id": r.id, "smiles": r.smiles, "openfe_abfe_dg_kcal": float("nan"),
                 "openfe_abfe_dg_sem": float("nan"), "openfe_status": status, "openfe_note": note}
                for r in df[cols].itertuples(index=False)]
        write(rows)
        print(f"{status}: {note}")
        return 0

    if not cfg.get("enabled", False):
        return all_status("SKIPPED", "openfe.enabled is false")

    of_exe = find_openfe_env_exe("openfe", cfg.get("exe"))
    of_py = find_openfe_env_exe("python", cfg.get("python"))
    if not of_exe or not of_py or "openfe" not in of_exe:
        return all_status("NOT_INSTALLED", "OpenFE environment not found (see docs/INSTALL.md)")

    # AmberTools (antechamber/sqm) must be on PATH for AM1-BCC charges, as if the env were activated.
    env = os.environ.copy()
    env["PATH"] = str(Path(of_exe).parent) + os.pathsep + env.get("PATH", "")
    env["CONDA_PREFIX"] = str(Path(of_exe).parent.parent)
    env.setdefault("OMP_NUM_THREADS", str(cfg.get("threads", 8)))  # sqm and OpenMM CPU parts oversubscribe otherwise

    n = int(cfg.get("max_input", 5))
    if "docking_score" in df.columns:
        df = df.assign(_s=pd.to_numeric(df["docking_score"], errors="coerce")).sort_values("_s").head(n)
    else:
        df = df.head(n)
    work = Path(a.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    cof = [os.path.abspath(os.path.expanduser(p)) for p in (cfg.get("cofactors_sdf") or [])]

    rows = []
    for r in df.itertuples(index=False):
        tag = _slug(r.id)
        d = work / tag
        d.mkdir(parents=True, exist_ok=True)
        row = {"id": r.id, "smiles": r.smiles, "openfe_abfe_dg_kcal": float("nan"),
               "openfe_abfe_dg_sem": float("nan"), "openfe_status": "NOT_READY", "openfe_note": ""}
        pose = getattr(r, "pose_sdf", None)
        if not (isinstance(pose, str) and os.path.isfile(pose)):
            row["openfe_note"] = "no bound pose (pose_sdf); ABFE needs a docked or predicted pose"
            rows.append(row)
            continue
        res = d / "result.json"
        if a.stage in ("all", "plan") and not (d / "abfe.json").is_file() and not res.is_file():
            cmd = [of_py, str(Path(__file__).resolve()), "--plan-one", "--id", str(r.id), "--ligand", pose,
                   "--protein", os.path.abspath(a.protein), "--dir", str(d),
                   "--repeats", str(cfg.get("protocol_repeats", 3)),
                   "--complex-ns", str(cfg.get("complex_production_ns", 5.0)),
                   "--solvent-ns", str(cfg.get("solvent_production_ns", 10.0)),
                   "--platform", str(cfg.get("platform", "") or ""),
                   "--preequil-ns", str(cfg.get("preequil_ns", 0) or 0),
                   "--charge-method", str(cfg.get("charge_method", "am1bcc"))]
            for c in cof:
                cmd += ["--cofactor", c]
            q = subprocess.run(cmd, capture_output=True, text=True, env=env)
            if q.returncode != 0:
                msg = (q.stdout.strip().splitlines() or q.stderr.strip().splitlines() or ["planning failed"])[-1]
                row["openfe_status"] = "NOT_READY" if q.returncode == 3 else "FAILED"
                row["openfe_note"] = msg[:300]
                rows.append(row)
                continue
        if a.stage in ("all", "run") and (d / "abfe.json").is_file() and not res.is_file():
            q = subprocess.run([of_exe, "quickrun", str(d / "abfe.json"), "-o", str(res), "-d", str(d)],
                               capture_output=True, text=True, env=env)
            (d / "quickrun.log").write_text((q.stdout or "") + "\n" + (q.stderr or ""))
        if res.is_file():
            dg, sem, st, note = read_result_json(res)
            row.update(openfe_abfe_dg_kcal=dg, openfe_abfe_dg_sem=sem, openfe_status=st, openfe_note=note)
        elif (d / "abfe.json").is_file():
            row["openfe_status"] = "NOT_RUN" if a.stage == "plan" else "FAILED"
            row["openfe_note"] = "planned; run `--stage run`" if a.stage == "plan" else "see quickrun.log"
        rows.append(row)
    out = write(rows)
    ok = int((out["openfe_status"] == "DONE").sum())
    print(f"openfe_abfe DONE {ok}/{len(out)}")
    if a.mmgbsa and os.path.isfile(a.mmgbsa):
        cmp_ = compare_with_mmgbsa(out, a.mmgbsa, Path(a.output).with_name("openfe_vs_mmgbsa.json"))
        if cmp_:
            print("comparison:", json.dumps(cmp_))
    return 0


if __name__ == "__main__":
    sys.exit(main())
