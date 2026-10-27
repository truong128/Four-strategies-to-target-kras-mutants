#!/usr/bin/env python
"""Download KRAS G12D PDB 7RPZ from RCSB and prepare the docking receptor.

Writes into structures/ (git-ignored, regenerated on demand):

    7RPZ.pdb                          raw RCSB entry
    7RPZ_KRAS_G12D_protein.pdb        chain A protein only (Boltz, OpenMM and MM/GBSA input)
    7RPZ_KRAS_G12D_with_GDP.pdb       chain A protein + GDP + Mg2+ (docking receptor)
    7RPZ_MRTX1133.pdb / .sdf          crystallographic ligand (6IC), used as reference and positive control
    KRAS_G12D.pdbqt                   Meeko receptor with the box from config.yaml
    KRAS_G12D.box.txt                 docking box (center and size)

The inhibitor is removed from the receptor. GDP and Mg are kept so that docking searches
the switch-II pocket of the inactive, nucleotide-bound state.

Run inside the pdac-dd environment:  python scripts/fetch_structures.py
(the script re-launches itself in the pdac-md environment for PDBFixer, then uses Meeko from pdac-dd)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from build_resistance_panel import (  # noqa: E402
    download, ensure_structure_python, fix_protein, het_lines, meeko_prepare, protein_lines, write_receptor_pdb,
)
from common import load_yaml  # noqa: E402
from envtools import env_tool  # noqa: E402


def ligand_pdb(text: str, resname: str = "6IC") -> list[str]:
    keep = []
    for ln in text.splitlines():
        if ln.startswith("HETATM") and ln[17:20].strip() == resname and ln[21:22] == "A" and ln[16] in " A":
            keep.append(ln[:16] + " " + ln[17:])
    return keep


def ligand_sdf(pdb_lines: list[str], out_sdf: Path, ccd: str = "6IC") -> bool:
    """Bond orders come from the RCSB ideal CCD SDF so that the pose stays crystallographic."""
    try:
        import urllib.request
        from rdkit import Chem
        from rdkit.Chem import AllChem

        tmp = out_sdf.with_suffix(".pdb")
        tmp.write_text("\n".join(pdb_lines + ["END"]) + "\n")
        url = f"https://files.rcsb.org/ligands/download/{ccd}_ideal.sdf"
        req = urllib.request.Request(url, headers={"User-Agent": "pdac-dd-e2e"})
        template = Chem.MolFromMolBlock(urllib.request.urlopen(req, timeout=60).read().decode(), removeHs=True)
        mol = Chem.MolFromPDBFile(str(tmp), removeHs=True, sanitize=False)
        mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        Chem.MolToMolFile(mol, str(out_sdf))
        return True
    except Exception as exc:  # reference ligand is optional (redocking QC only)
        print(f"note: reference ligand SDF not written ({exc})")
        return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=str(ROOT / "config.yaml"))
    ap.add_argument("--pdb-id", default="7RPZ")
    ap.add_argument("--output-dir", default=str(ROOT / "structures"))
    ap.add_argument("--skip-pdbqt", action="store_true")
    a = ap.parse_args()

    cfg = load_yaml(a.config)
    box = cfg["target"]["box"]
    center = [float(x) for x in box["center"]]
    size = [float(x) for x in box["size"]]
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    raw = out / f"{a.pdb_id}.pdb"
    download(f"https://files.rcsb.org/download/{a.pdb_id}.pdb", raw)
    text = raw.read_text()

    ensure_structure_python()  # PDBFixer lives in pdac-md; re-executes this script there if needed
    raw_prot = out / "7RPZ_KRAS_G12D_protein.raw.pdb"
    raw_prot.write_text("\n".join(protein_lines(text) + ["END"]) + "\n")
    protein = out / "7RPZ_KRAS_G12D_protein.pdb"
    fix_protein(raw_prot, protein)  # adds missing heavy atoms, keeps residue ids; no loop rebuilding
    raw_prot.unlink()
    cofactors, dropped = het_lines(text)
    write_receptor_pdb(protein.read_text(), cofactors, out / "7RPZ_KRAS_G12D_with_GDP.pdb")
    print(f"cofactor atoms {len(cofactors)}; ignored hets {sorted(dropped)}")

    lig = ligand_pdb(text)
    if lig:
        (out / "7RPZ_MRTX1133.pdb").write_text("\n".join(lig + ["END"]) + "\n")
        ligand_sdf(lig, out / "7RPZ_MRTX1133.sdf")

    (out / "KRAS_G12D.box.txt").write_text(
        "".join(f"center_{k} = {v:.3f}\n" for k, v in zip("xyz", center))
        + "".join(f"size_{k} = {v:.3f}\n" for k, v in zip("xyz", size))
    )

    if not a.skip_pdbqt:
        meeko = cfg.get("resistance", {}).get("meeko_receptor_exe", "mk_prepare_receptor.py")
        meeko = env_tool("pdac-dd", meeko)
        meeko_prepare(out / "7RPZ_KRAS_G12D_with_GDP.pdb", out / "KRAS_G12D", center, size, meeko)
        print("receptor ready:", out / "KRAS_G12D.pdbqt")


if __name__ == "__main__":
    main()
