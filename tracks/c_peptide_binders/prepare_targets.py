#!/usr/bin/env python3
"""Build all 18 KRAS allele targets and an RFD4 condition spec for short blocking peptides.

The scaffold is the KRAS G12D G domain from PDB 7RPZ (chain A, residues 1-169),
with GDP, magnesium, and MRTX1133 already removed. G12D keeps the crystal side
chain. Every other allele is the same backbone with an ideal side chain grafted
onto the mutated residue. Alleles away from codon 12 also put Gly back at 12,
so each target is a single mutant. Wild type is Gly12 on that same backbone.

The designed chain is a new 8-16 residue peptide. The target backbone is one
rigid group and its sequence is fixed. Residue 2 has no coordinates in 7RPZ,
which leaves residue 1 disconnected, so the target contig is A3-169. Hotspots
are the mutant side chain plus the two nearest side chains at that site.
Codons 12, 13, and 61 use the P-loop / switch-II pocket. A146 and K117 use the
local nucleotide-site residues instead of that pocket.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import biotite.structure.info as info
import gemmi
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent
# Produced by tracks/a_small_molecule_discovery/scripts/fetch_structures.py (scripts/fetch_public_inputs.sh)
DEFAULT_TEMPLATE = (
    ROOT.parent
    / "a_small_molecule_discovery"
    / "structures"
    / "7RPZ_KRAS_G12D_protein.pdb"
)
BACKBONE = {"N", "CA", "C", "O"}
# P-loop and switch II, the pocket that contains codon 12/13/61.
SWITCH_POCKET = {10, 11, 12, 13, 14, 16, 34, 58, 59, 60, 61, 62, 63, 64, 68, 72, 95, 96, 99}
SWITCH_AIM = {16, 59, 61, 95, 96, 99}
# G4 motif around Lys117, on the nucleotide-base side of the G domain.
G4_POCKET = {14, 15, 83, 84, 85, 116, 118, 119, 144, 145, 146, 147}
G4_AIM = {14, 15, 83, 85, 145, 146}
# G5 motif around Ala146. The aim set points into the nucleotide cleft
# rather than along the 145-147 backbone.
G5_POCKET = {15, 16, 18, 19, 22, 28, 116, 117, 144, 145, 147, 148, 149}
G5_AIM = {15, 18, 19, 28, 116, 117}
SITE_CONTEXT = {
    12: {"pocket": SWITCH_POCKET, "aim": SWITCH_AIM},
    13: {"pocket": SWITCH_POCKET, "aim": SWITCH_AIM},
    61: {"pocket": SWITCH_POCKET, "aim": SWITCH_AIM},
    117: {"pocket": G4_POCKET, "aim": G4_AIM},
    146: {"pocket": G5_POCKET, "aim": G5_AIM},
}
TERMINAL = {
    "GLY": ["CA"],
    "ALA": ["CB"],
    "VAL": ["CG1", "CG2"],
    "LEU": ["CD1", "CD2"],
    "PHE": ["CZ"],
    "CYS": ["SG"],
    "SER": ["OG"],
    "THR": ["OG1", "CG2"],
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
    "ASN": ["OD1", "ND2"],
    "GLN": ["OE1", "NE2"],
    "HIS": ["ND1", "NE2"],
    "LYS": ["NZ"],
    "ARG": ["NH1", "NH2"],
}
# Role of the mutant side chain in a hydrogen bond the peptide should satisfy.
HBOND = {
    "ASP": ("ACCEPTOR", ["OD1", "OD2"]),
    "GLU": ("ACCEPTOR", ["OE1", "OE2"]),
    "SER": ("BOTH", ["OG"]),
    "THR": ("BOTH", ["OG1"]),
    "CYS": ("BOTH", ["SG"]),
    "ASN": ("BOTH", ["OD1", "ND2"]),
    "GLN": ("BOTH", ["OE1", "NE2"]),
    "HIS": ("BOTH", ["ND1", "NE2"]),
    "LYS": ("DONOR", ["NZ"]),
    "ARG": ("DONOR", ["NH1", "NH2"]),
}
# Same 18 names as the FDA allele screen, in that screen's column order.
# name, residue number, amino acid after the edit, crystal (no graft) or not
ALLELES = [
    ("G12D", 12, "ASP", True),
    ("WT", 12, "GLY", False),
    ("G12C", 12, "CYS", False),
    ("G12V", 12, "VAL", False),
    ("G12R", 12, "ARG", False),
    ("G12A", 12, "ALA", False),
    ("G13D", 13, "ASP", False),
    ("A146T", 146, "THR", False),
    ("G12S", 12, "SER", False),
    ("G12F", 12, "PHE", False),
    ("G13C", 13, "CYS", False),
    ("G13R", 13, "ARG", False),
    ("Q61H", 61, "HIS", False),
    ("Q61L", 61, "LEU", False),
    ("Q61R", 61, "ARG", False),
    ("Q61K", 61, "LYS", False),
    ("A146V", 146, "VAL", False),
    ("K117N", 117, "ASN", False),
]
DEFAULT_ALLELES = [row[0] for row in ALLELES]


def _clean(name: str) -> str:
    return name.strip()


def _element(atom_name: str) -> gemmi.Element:
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if letters[:2] in {"CL", "BR", "FE", "MG", "ZN"}:
        return gemmi.Element(letters[:2].title())
    return gemmi.Element(letters[0])


def _kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (target - target_center)
    u_mat, _, vt_mat = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(vt_mat.T @ u_mat.T))
    rotation = vt_mat.T @ np.diag([1.0, 1.0, sign]) @ u_mat.T
    translation = target_center - rotation @ mobile_center
    return rotation, translation


def _ideal_heavy_atoms(resname: str) -> dict[str, np.ndarray]:
    residue = info.residue(resname)
    coords = {}
    for atom_name, coord in zip(residue.atom_name, residue.coord, strict=True):
        if atom_name[0] == "H" or atom_name in {"OXT", "HXT"}:
            continue
        coords[atom_name] = np.asarray(coord, dtype=float)
    return coords


def _is_polymer(res: gemmi.Residue) -> bool:
    # gemmi marks ATOM records as het_flag "A" and HETATM as "H".
    return res.het_flag != "H"


def _residue(chain: gemmi.Chain, number: int) -> gemmi.Residue:
    matches = [res for res in chain if res.seqid.num == number and _is_polymer(res)]
    if len(matches) != 1:
        raise SystemExit(f"chain {chain.name} residue {number}: expected 1 match, found {len(matches)}")
    return matches[0]


def target_contig(chain: gemmi.Chain) -> str:
    """Preserved KRAS segments. 7RPZ has no coordinates for residue 2, so a single A1-169 range is invalid."""
    numbers = sorted(res.seqid.num for res in chain if _is_polymer(res))
    spans: list[tuple[int, int]] = []
    start = previous = numbers[0]
    for number in numbers[1:]:
        if number == previous + 1:
            previous = number
            continue
        spans.append((start, previous))
        start = previous = number
    spans.append((start, previous))
    parts = [f"A{begin}" if begin == end else f"A{begin}-{end}" for begin, end in spans]
    return ",".join(parts)


def _coord(res: gemmi.Residue, atom_name: str) -> np.ndarray:
    for atom in res:
        if _clean(atom.name) == atom_name:
            return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)
    raise SystemExit(f"{res.name} {res.seqid.num} is missing {atom_name}")


def _drop_side_chain(res: gemmi.Residue) -> None:
    for index in range(len(res) - 1, -1, -1):
        if _clean(res[index].name) not in BACKBONE:
            del res[index]


def graft_side_chain(res: gemmi.Residue, resname: str) -> None:
    """Replace the side chain. Backbone N, CA, C, O stay on the crystal coordinates."""
    if resname == "GLY":
        _drop_side_chain(res)
        res.name = "GLY"
        return

    ideal = _ideal_heavy_atoms(resname)
    mobile = np.vstack([ideal[name] for name in ("N", "CA", "C")])
    target = np.vstack([_coord(res, name) for name in ("N", "CA", "C")])
    rotation, translation = _kabsch(mobile, target)
    _drop_side_chain(res)
    for atom_name, coord in ideal.items():
        if atom_name in BACKBONE:
            continue
        placed = rotation @ coord + translation
        atom = gemmi.Atom()
        atom.name = atom_name
        atom.element = _element(atom_name)
        atom.pos = gemmi.Position(float(placed[0]), float(placed[1]), float(placed[2]))
        atom.occ = 1.0
        atom.b_iso = 20.0
        res.add_atom(atom)
    res.name = resname


def _aim_side_chain_at_pocket(chain: gemmi.Chain, res: gemmi.Residue, aim: set[int]) -> None:
    """Rotate chi1 so a grafted side chain points at the local mutation site."""
    terminal_names = [name for name in TERMINAL.get(res.name, []) if name not in {"CA", "CB"}]
    if not terminal_names or not any(_clean(atom.name) == "CB" for atom in res):
        return
    origin = _coord(res, "CA")
    axis = _coord(res, "CB") - origin
    axis_norm = np.linalg.norm(axis)
    if axis_norm < 1e-6:
        return
    axis = axis / axis_norm
    points = []
    for other in chain:
        if not _is_polymer(other) or other.seqid.num not in aim:
            continue
        if other.seqid.num == res.seqid.num:
            continue
        try:
            points.append(_coord(other, "CA"))
        except SystemExit:
            continue
    if len(points) < 3:
        return
    center = np.mean(np.vstack(points), axis=0)

    def rotate(point: np.ndarray, angle: float) -> np.ndarray:
        vector = point - origin
        cosine, sine = np.cos(angle), np.sin(angle)
        return (
            origin
            + vector * cosine
            + np.cross(axis, vector) * sine
            + axis * np.dot(axis, vector) * (1.0 - cosine)
        )

    current = {_clean(atom.name): _xyz(atom) for atom in res}
    best_angle = 0.0
    best_score = None
    for step in range(36):
        angle = np.deg2rad(step * 10.0)
        score = float(np.mean([
            np.linalg.norm(rotate(current[name], angle) - center)
            for name in terminal_names
            if name in current
        ]))
        if best_score is None or score < best_score:
            best_score = score
            best_angle = angle
    for atom in res:
        name = _clean(atom.name)
        if name in BACKBONE or name == "CB":
            continue
        placed = rotate(current[name], best_angle)
        atom.pos = gemmi.Position(float(placed[0]), float(placed[1]), float(placed[2]))


def _side_chain_atoms(res: gemmi.Residue) -> list[gemmi.Atom]:
    return [atom for atom in res if _clean(atom.name) not in BACKBONE | {"OXT"}]


def _xyz(atom: gemmi.Atom) -> np.ndarray:
    return np.array([atom.pos.x, atom.pos.y, atom.pos.z], dtype=float)


def hotspot_atoms(chain: gemmi.Chain, site: int, pocket: set[int]) -> list[tuple[int, str, str]]:
    """Mutant terminal atoms, plus one side-chain atom from each of the two nearest pocket residues."""
    mutant = _residue(chain, site)
    terminal_names = [name for name in TERMINAL[mutant.name] if any(_clean(a.name) == name for a in mutant)]
    if not terminal_names:
        terminal_names = ["CA"]
    terminals = [_coord(mutant, name) for name in terminal_names]
    chosen = [(site, mutant.name, name) for name in terminal_names]

    nearest: list[tuple[float, int, str, str]] = []
    for res in chain:
        if not _is_polymer(res) or res.seqid.num not in pocket or res.seqid.num == site:
            continue
        for atom in _side_chain_atoms(res):
            distance = min(float(np.linalg.norm(_xyz(atom) - point)) for point in terminals)
            nearest.append((distance, res.seqid.num, res.name, _clean(atom.name)))
    nearest.sort()
    seen_residues: set[int] = set()
    for distance, number, resname, atom_name in nearest:
        # Arg side chains sit just past 8 A from the second neighbor (G13R 8.6, Q61R 8.8).
        if distance > 9.0 or number in seen_residues:
            continue
        seen_residues.add(number)
        chosen.append((number, resname, atom_name))
        if len(seen_residues) == 2:
            break
    return chosen


def select_for(atoms: list[tuple[int, str, str]]) -> list[str]:
    grouped: dict[int, list[str]] = {}
    for number, _resname, atom_name in atoms:
        grouped.setdefault(number, []).append(atom_name)
    return [f"sel('A/*/{number}/[{','.join(names)}]')" for number, names in grouped.items()]


def build_structure(template: Path, allele: str, site: int, resname: str, crystal: bool) -> gemmi.Structure:
    structure = gemmi.read_structure(str(template))
    model = structure[0]
    chain = model["A"] if "A" in [ch.name for ch in model] else model[0]
    chain.name = "A"
    # Residue 2 is absent, so residue 1 is not bonded to the rest of the chain.
    # Drop it. The pocket numbering (12, 13, 61) is unchanged.
    for index in range(len(chain) - 1, -1, -1):
        if chain[index].seqid.num < 3:
            del chain[index]
    # Single mutants only. The template is G12D, so other alleles put Gly back at 12.
    if not crystal and site != 12:
        graft_side_chain(_residue(chain, 12), "GLY")
    if not crystal:
        graft_side_chain(_residue(chain, site), resname)
        _aim_side_chain_at_pocket(chain, _residue(chain, site), SITE_CONTEXT[site]["aim"])
    elif _residue(chain, site).name != resname:
        raise SystemExit(f"{allele}: template residue {site} is {_residue(chain, site).name}, expected {resname}")
    structure.setup_entities()
    return structure


def write_cif(structure: gemmi.Structure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.setup_entities()
    for model in structure:
        for chain in model:
            for res in chain:
                # Atomworks reads label_asym_id and label_seq_id. Keep both on the crystal numbering.
                res.subchain = chain.name
                res.entity_id = "1"
                res.label_seq = res.seqid.num
    structure.make_mmcif_document().write_file(str(path))


def condition_block(hot_selects: list[str], hbond: tuple[str, list[str]] | None, site: int) -> dict:
    conditions: dict = {
        "C_CRDGRP": {"1": [{"select": "sel('A/*/*/[N,CA,C,O]') and not has_nan_coord()"}]},
        "C_SEQ": {"true": [{"select": "sel('A')"}]},
        "C_HOT": {"true": [{"select": item} for item in hot_selects]},
    }
    if hbond is not None:
        role, names = hbond
        present = [name for name in names if any(name in item for item in hot_selects)]
        if present:
            conditions["C_HB"] = {role: [{"select": f"sel('A/*/{site}/[{','.join(present)}]')"}]}
    return conditions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path, default=DEFAULT_TEMPLATE)
    parser.add_argument("--alleles", default=",".join(DEFAULT_ALLELES), help="Comma-separated allele names. Default is all 18.")
    parser.add_argument("--peptide", default="8-16", help="Designed peptide length range, inclusive.")
    args = parser.parse_args()
    requested = [item.strip() for item in args.alleles.split(",") if item.strip()]
    known = {row[0]: row for row in ALLELES}
    unknown = [name for name in requested if name not in known]
    if unknown:
        raise SystemExit(f"Unknown alleles {unknown}. Choices: {', '.join(known)}")

    targets = ROOT / "targets"
    specs = ROOT / "specs"
    specs.mkdir(parents=True, exist_ok=True)
    designs: dict[str, dict] = {}
    manifest_lines = ["allele\tsite\tamino_acid\tsource\thotspots\thbond\tcif"]

    for name in requested:
        _allele, site, resname, crystal = known[name]
        structure = build_structure(args.template, name, site, resname, crystal)
        cif_path = targets / f"KRAS_{name}.cif"
        write_cif(structure, cif_path)
        atoms = hotspot_atoms(structure[0]["A"], site, SITE_CONTEXT[site]["pocket"])
        selects = select_for(atoms)
        hbond = HBOND.get(resname) if site == known[name][1] else None
        # VAL has no polar side chain; do not invent a hydrogen-bond condition.
        if resname not in HBOND:
            hbond = None
        designs[f"KRAS_{name}_peptide"] = {
            "input": f"../targets/{cif_path.name}",
            "contig": f"{args.peptide}, /0, {target_contig(structure[0]['A'])}",
            "conditions": condition_block(selects, hbond, site),
        }
        if crystal:
            source = "7RPZ crystal side chain"
        elif resname == "GLY":
            source = "7RPZ backbone, Gly restored at residue 12"
        else:
            source = "7RPZ backbone, ideal side chain grafted at the mutation"
        hotspot_text = "; ".join(f"{res}{num}.{atom}" for num, res, atom in atoms)
        hbond_text = "" if hbond is None else f"{hbond[0]} {' '.join(hbond[1])}"
        manifest_lines.append(f"{name}\t{site}\t{resname}\t{source}\t{hotspot_text}\t{hbond_text}\t{cif_path.name}")
        print(f"{name}: residue {site} = {resname}; hotspots {hotspot_text}")

    spec_path = specs / "kras_mutant_peptides.yaml"
    spec_path.write_text(
        "# Short peptides aimed at each KRAS allele's mutation site.\n"
        "# Generated by prepare_targets.py. Paths resolve from this file's directory.\n"
        + yaml.safe_dump(designs, sort_keys=False)
    )
    (targets / "manifest.tsv").write_text("\n".join(manifest_lines) + "\n")
    print(f"wrote {spec_path}")

    from rfproteina.inference_utils.condition_spec import build_annotated_atom_array_from_condition_spec

    for example_id, spec in designs.items():
        built = build_annotated_atom_array_from_condition_spec(
            spec, spec_dir=specs, example_id=example_id, validate=True
        )
        chains = sorted(set(built["atom_array"].chain_id))
        print(f"validated {example_id}: chains {chains}, atoms {len(built['atom_array'])}")


if __name__ == "__main__":
    main()
