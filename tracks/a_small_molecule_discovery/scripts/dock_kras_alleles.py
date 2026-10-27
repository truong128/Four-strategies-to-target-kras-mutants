#!/usr/bin/env python
"""Dock the final KRAS G12D leads into one GDP monomer per common KRAS allele.

Experimental receptors are aligned onto 7RPZ and docked in the same switch-II
box used for the G12D screen. Alleles that lack a deposited GDP structure with
an intact switch II are built by one side-chain substitution on wild-type 6MBU
and labeled modeled_from_6MBU. Every receptor keeps chain A only. 6BOF is a
crystal homodimer, and its second chain stays outside this switch-II box.

The large library screen used AutoDock Vina. This shortlist uses GNINA 1.3
with --cnn_scoring all, so the convolutional network scores the whole search.
CNNaffinity is a predicted pK (higher means tighter). On G12C the pose is
still noncovalent.

Prepare receptors, then dock:

    python scripts/dock_kras_alleles.py --all

Or separately:

    python scripts/dock_kras_alleles.py --prepare
    python scripts/dock_kras_alleles.py --dock
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
from envtools import env_tool

MD_PYTHON = os.environ.get("PDAC_MD_PYTHON") or env_tool("pdac-md", "python")
DD_PYTHON = os.environ.get("PDAC_DD_PYTHON") or env_tool("pdac-dd", "python")
REFERENCE = "7RPZ"
PARENT = "6MBU"
WT_SITES = {12: "GLY", 13: "GLY", 61: "GLN", 117: "LYS", 146: "ALA"}
ALLELE_SITE = {
    "G12C": (12, "CYS"), "G12D": (12, "ASP"), "G12V": (12, "VAL"), "G12R": (12, "ARG"),
    "G12A": (12, "ALA"), "G12S": (12, "SER"), "G12F": (12, "PHE"),
    "G13D": (13, "ASP"), "G13C": (13, "CYS"), "G13R": (13, "ARG"),
    "Q61H": (61, "HIS"), "Q61L": (61, "LEU"), "Q61R": (61, "ARG"), "Q61K": (61, "LYS"),
    "A146T": (146, "THR"), "A146V": (146, "VAL"), "K117N": (117, "ASN"),
}


def python_for(kind):
    path = MD_PYTHON if kind == "md" else DD_PYTHON
    if not os.path.isfile(path):
        raise SystemExit(f"Missing Python for {kind}: {path}")
    return path


def has_modules(names):
    for name in names:
        try:
            __import__(name)
        except ImportError:
            return False
    return True


def residue_names(pdb_path):
    out = {}
    for ln in Path(pdb_path).read_text().splitlines():
        if not ln.startswith("ATOM") or ln[21:22] != "A":
            continue
        resi = int(ln[22:26])
        out.setdefault(resi, ln[17:20].strip())
    return out


def expected_sites(allele):
    sites = dict(WT_SITES)
    change = ALLELE_SITE.get(allele)
    if change:
        sites[change[0]] = change[1]
    return sites


def assert_allele(pdb_path, allele):
    names = residue_names(pdb_path)
    missing = [i for i in list(range(30, 41)) + list(range(60, 77)) if i not in names]
    if missing:
        raise SystemExit(f"{pdb_path.name} is missing switch residues {missing}")
    for resi, aa in expected_sites(allele).items():
        got = names.get(resi)
        if got != aa:
            raise SystemExit(f"{pdb_path.name} allele {allele}: residue {resi} is {got}, expected {aa}")
    return names


def prepare(panel_path, config_path, out):
    sys.path.insert(0, str(ROOT / "scripts"))
    from build_resistance_panel import (
        align_pdb, download, fix_protein, het_lines, meeko_prepare,
        protein_lines, relax_sidechains, write_receptor_pdb,
    )
    import yaml
    from Bio.PDB import PDBParser

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    box = cfg["target"]["box"]
    center = [float(x) for x in box["center"]]
    size = [float(x) for x in box["size"]]
    meeko = cfg.get("resistance", {}).get("meeko_receptor_exe", "mk_prepare_receptor.py")
    if not os.path.isfile(meeko) and not shutil_which(meeko):
        meeko = env_tool("pdac-dd", meeko)
    rows = list(csv.DictReader(open(panel_path), delimiter="\t"))
    if not any(r["pdb_id"] == REFERENCE for r in rows):
        raise SystemExit(f"Panel must include reference {REFERENCE}")
    experimental = [r for r in rows if r["source_type"] == "experimental"]
    modeled = [r for r in rows if r["source_type"] != "experimental"]

    raw, aligned, receptors = out / "raw", out / "aligned", out / "receptors"
    for p in (raw, aligned, receptors):
        p.mkdir(parents=True, exist_ok=True)
    for r in experimental:
        download(f"https://files.rcsb.org/download/{r['pdb_id']}.pdb", raw / f"{r['pdb_id']}.pdb")

    parser = PDBParser(QUIET=True)
    ref = parser.get_structure("ref", str(raw / f"{REFERENCE}.pdb"))
    from build_resistance_panel import ca_map
    refmap = ca_map(ref, "A")
    metrics = []
    for r in experimental:
        pid = r["pdb_id"]
        rms, n = align_pdb(raw / f"{pid}.pdb", aligned / f"{pid}_aligned.pdb", refmap, parser)
        metrics.append({"pdb_id": pid, "allele": r["allele"], "rmsd_to_7RPZ_CA": f"{rms:.4f}", "n_ca": n})
        print(f"aligned {pid} {r['allele']} rmsd {rms:.3f} n_ca {n}")

    manifest = []
    for r in experimental:
        built = build_from_aligned(r, aligned / f"{r['pdb_id']}_aligned.pdb", None, receptors, center, size, meeko,
                                   fix_protein, het_lines, protein_lines, relax_sidechains, write_receptor_pdb, meeko_prepare)
        manifest.append(built)
    parent = aligned / f"{PARENT}_aligned.pdb"
    if not parent.exists():
        raise SystemExit(f"Wild-type parent {PARENT} was not aligned")
    for r in modeled:
        if not r["mutation"]:
            raise SystemExit(f"{r['pdb_id']} needs a mutation")
        built = build_from_aligned(r, parent, r["mutation"], receptors, center, size, meeko,
                                   fix_protein, het_lines, protein_lines, relax_sidechains, write_receptor_pdb, meeko_prepare)
        manifest.append(built)

    (out / "panel_manifest.json").write_text(json.dumps({
        "center": center, "size": size, "reference_pdb": REFERENCE,
        "note": "Modeled alleles are one substitution on wild-type 6MBU. They are not crystal structures.",
        "receptors": manifest,
    }, indent=2))
    with open(out / "alignment_metrics.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["pdb_id", "allele", "rmsd_to_7RPZ_CA", "n_ca"])
        writer.writeheader()
        writer.writerows(metrics)
    print("allele receptors ready", len(manifest))


def shutil_which(name):
    from shutil import which
    return which(name)


def build_from_aligned(row, src, mutation, receptors, center, size, meeko,
                       fix_protein, het_lines, protein_lines, relax_sidechains, write_receptor_pdb, meeko_prepare):
    pid = row["pdb_id"]
    text = src.read_text()
    cofactors, dropped = het_lines(text)
    if mutation:
        fixed = receptors / f"{pid}_fixed.pdb"
        fix_protein(src, fixed, mutation=mutation)
        focus = int(mutation.split("-")[1])
        relax_sidechains(fixed, cofactors, focus)
        body = fixed.read_text()
    else:
        prot = receptors / f"{pid}_protein.pdb"
        prot.write_text("\n".join(protein_lines(text) + ["END"]) + "\n")
        fixed = receptors / f"{pid}_fixed.pdb"
        fix_protein(prot, fixed)
        # PDBFixer builds missing side chains in place. A short restrained
        # minimization moves those atoms off the rest of the protein. CA stays fixed.
        relax_sidechains(fixed, cofactors, None)
        body = fixed.read_text()
    receptor_pdb = receptors / f"{pid}.pdb"
    write_receptor_pdb(body, cofactors, receptor_pdb)
    assert_allele(receptor_pdb, row["allele"])
    meeko_prepare(receptor_pdb, receptors / pid, center, size, meeko)
    kept = sorted({ln[17:20].strip() for ln in cofactors})
    print("prepared", pid, row["allele"], row["source_type"], "kept", kept)
    return {
        "pdb_id": pid,
        "allele": row["allele"],
        "state": row["state"],
        "source_type": row["source_type"],
        "mutation": mutation or "",
        "pdbqt": f"receptors/{pid}.pdbqt",
        "notes": row["notes"],
        "dropped_hets": ",".join(sorted(dropped)),
    }


def _gnina_one(task):
    import re
    import shutil
    import subprocess
    import tempfile
    smiles, mid, extra, receptor, center, size, exh, modes, gnina, sdf_script, no_gpu = task
    score_re = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", re.M)
    cnn_re = re.compile(r"CNNscore[:\s=]+(-?\d+(?:\.\d+)?)", re.I)
    aff_re = re.compile(r"CNNaffinity[:\s=]+(-?\d+(?:\.\d+)?)", re.I)
    td = tempfile.mkdtemp(prefix="gnina_")
    try:
        sdf = os.path.join(td, "lig.sdf")
        outp = os.path.join(td, "out.sdf")
        embed = subprocess.run(
            [sys.executable, sdf_script, "--smiles", str(smiles), "--id", str(mid), "--output", sdf],
            text=True, capture_output=True,
        )
        if embed.returncode:
            status, vina_aff, cnn, pk = "embed_fail", float("nan"), float("nan"), float("nan")
        else:
            cmd = [
                gnina, "-r", receptor, "-l", sdf, "-o", outp,
                "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
                "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
                "--cnn_scoring", "all",
                "--exhaustiveness", str(exh), "--num_modes", str(modes),
                "--seed", "42", "--cpu", "4",
            ]
            if no_gpu:
                cmd.append("--no_gpu")
            q = subprocess.run(cmd, text=True, capture_output=True)
            text = (q.stdout or "") + "\n" + (q.stderr or "")
            ms, mc, ma = score_re.search(text), cnn_re.search(text), aff_re.search(text)
            vina_aff = float(ms.group(1)) if ms else float("nan")
            cnn = float(mc.group(1)) if mc else float("nan")
            pk = float(ma.group(1)) if ma else float("nan")
            status = "ok" if ma else "parse_fail"
            if status != "ok":
                tail = " | ".join(text.strip().splitlines()[-6:])
                print(
                    f"gnina_fail {extra.get('allele')} {mid} rc={q.returncode} {tail}",
                    flush=True,
                )
        row = {
            "smiles": smiles, "id": mid, "dock_status": status,
            "gnina_affinity_kcal": vina_aff, "gnina_cnn_score": cnn, "gnina_cnn_affinity": pk,
        }
        row.update(extra)
        return row
    finally:
        shutil.rmtree(td, ignore_errors=True)


def dock(panel_path, config_path, leads_path, out, workers, exhaustiveness, no_gpu):
    sys.path.insert(0, str(ROOT / "scripts"))
    import pandas as pd
    from common import load_yaml, read_any_table
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    from shutil import which

    cfg = load_yaml(config_path)
    box = cfg["target"]["box"]
    center = tuple(float(x) for x in box["center"])
    size = tuple(float(x) for x in box["size"])
    gnina = cfg.get("gnina", {}).get("exe", "gnina")
    if which(gnina) is None and os.environ.get("GNINA_EXE"):
        gnina = os.environ["GNINA_EXE"]
    if which(gnina) is None and not os.path.isfile(gnina):
        raise SystemExit(f"GNINA not found: {gnina}")
    exh = int(exhaustiveness or 16)
    leads = read_any_table(leads_path)
    if not {"id", "smiles"} <= set(leads.columns):
        raise SystemExit("Leads table needs id and smiles")
    manifest = json.loads((out / "panel_manifest.json").read_text())
    sdf_script = str(ROOT / "scripts" / "smiles_to_sdf.py")
    tasks = []
    for rec in manifest["receptors"]:
        receptor = str((out / rec["pdbqt"]).resolve())
        if not os.path.isfile(receptor):
            raise SystemExit(f"Missing receptor {receptor}. Run --prepare first.")
        extra_base = {"allele": rec["allele"], "pdb_id": rec["pdb_id"], "source_type": rec["source_type"]}
        for lead in leads.itertuples(index=False):
            tasks.append((lead.smiles, lead.id, dict(extra_base), receptor, center, size, exh, 9, gnina, sdf_script, no_gpu))
    n_workers = 1 if not no_gpu else max(1, min(int(workers), len(tasks)))
    print(
        f"GNINA cnn_scoring=all, {len(leads)} leads x {len(manifest['receptors'])} alleles, "
        f"workers {n_workers}, exhaustiveness {exh}, gpu {not no_gpu}",
        flush=True,
    )
    progress_path = out / "allele_dock_progress.csv"
    if n_workers == 1:
        rows = []
        for i, t in enumerate(tasks, 1):
            row = _gnina_one(t)
            rows.append(row)
            pd.DataFrame(rows).to_csv(progress_path, index=False)
            print(
                f"[{i}/{len(tasks)}] {row['dock_status']} {row['allele']} {row['pdb_id']} {row['id']} "
                f"vina={row['gnina_affinity_kcal']} cnn_affinity={row['gnina_cnn_affinity']}",
                flush=True,
            )
    else:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=mp.get_context("spawn")) as pool:
            rows = list(pool.map(_gnina_one, tasks, chunksize=1))
    long_df = pd.DataFrame(rows)
    g12d = long_df.loc[long_df["allele"] == "G12D", ["id", "gnina_cnn_affinity"]].rename(
        columns={"gnina_cnn_affinity": "G12D_cnn_affinity"}
    )
    long_df = long_df.merge(g12d, on="id", how="left")
    long_df["delta_cnn_affinity_vs_G12D"] = long_df["gnina_cnn_affinity"] - long_df["G12D_cnn_affinity"]
    long_path = out / "allele_dock_long.csv"
    long_df.to_csv(long_path, index=False)
    wide = long_df.pivot(index="id", columns="allele", values="gnina_cnn_affinity")
    wide.to_csv(out / "allele_dock_wide.csv")
    ok = int((long_df["dock_status"] == "ok").sum())
    print(f"docked {ok}/{len(long_df)}")
    print("wrote", long_path)
    for lig, part in long_df.groupby("id"):
        held = part.loc[part["delta_cnn_affinity_vs_G12D"] >= -1.0, "allele"].tolist()
        print(f"{lig} within 1 pK unit of its G12D CNNaffinity: {', '.join(held)}")


def rerun(kind, extra):
    py = python_for(kind)
    cmd = [py, str(Path(__file__).resolve()), *extra]
    print("exec", " ".join(cmd))
    subprocess.check_call(cmd)


def main():
    ap = argparse.ArgumentParser(description="Dock the four leads against common KRAS alleles")
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--dock", action="store_true")
    mode.add_argument("--all", action="store_true")
    ap.add_argument("--panel", default=str(ROOT / "config" / "kras_allele_panel.tsv"))
    ap.add_argument("--config", default=str(ROOT / "config.yaml"))
    ap.add_argument("--leads", default=str(ROOT / "results" / "cluster" / "diverse_hits.csv"))
    ap.add_argument("--output-dir", default=str(ROOT / "results" / "allele_dock"))
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--no-gpu", action="store_true")
    a = ap.parse_args()
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    shared = ["--panel", a.panel, "--config", a.config, "--output-dir", a.output_dir]
    if a.all:
        rerun("md", ["--prepare", *shared])
        dock_args = ["--dock", *shared, "--leads", a.leads, "--workers", str(a.workers),
                     "--exhaustiveness", str(a.exhaustiveness)]
        if a.no_gpu:
            dock_args.append("--no-gpu")
        rerun("dd", dock_args)
        return
    if a.prepare:
        if not has_modules(["pdbfixer", "Bio"]):
            os.execv(python_for("md"), [python_for("md"), *sys.argv])
        prepare(Path(a.panel), Path(a.config), out)
        return
    if not has_modules(["rdkit"]):
        os.execv(python_for("dd"), [python_for("dd"), *sys.argv])
    dock(Path(a.panel), Path(a.config), Path(a.leads), out, a.workers, a.exhaustiveness, a.no_gpu)


if __name__ == "__main__":
    main()
