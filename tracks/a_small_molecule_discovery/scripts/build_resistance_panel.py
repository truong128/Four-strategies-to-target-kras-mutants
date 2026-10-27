#!/usr/bin/env python
"""Build the KRAS switch-II resistance panel and prepare each receptor with Meeko.

Experimental structures are aligned onto 7RPZ so they share the MRTX1133 box.
Modeled mutants are built from 7RPZ and are stress tests, not experimental structures.
Nucleotide and Mg stay in the receptor. The bound inhibitor is removed.
Missing side-chain atoms are filled. Missing loops are not rebuilt.
"""
import argparse, csv, json, os, re, shutil, subprocess, sys, urllib.request
from pathlib import Path

from envtools import env_tool

MD_PYTHON = os.environ.get("PDAC_MD_PYTHON") or env_tool("pdac-md", "python")
KEEP_HETS = {"GDP", "GCP", "GNP", "GTP", "MG", "GSP"}
MUTATIONS = {
    "MODEL_V9E": "VAL-9-GLU",
    "MODEL_V9W": "VAL-9-TRP",
    "MODEL_V9Q": "VAL-9-GLN",
    "MODEL_G13P": "GLY-13-PRO",
    "MODEL_T58Y": "THR-58-TYR",
    "MODEL_R68G": "ARG-68-GLY",
    "MODEL_Y96W": "TYR-96-TRP",
    "MODEL_H95Q": "HIS-95-GLN",
    "MODEL_H95L": "HIS-95-LEU",
}


def ensure_structure_python():
    try:
        import pdbfixer  # noqa: F401
        import Bio  # noqa: F401
    except ImportError:
        if not os.path.isfile(MD_PYTHON):
            raise SystemExit(f"pdbfixer is required and {MD_PYTHON} is missing")
        os.execv(MD_PYTHON, [MD_PYTHON, *sys.argv])


def download(url, dest: Path):
    if dest.exists() and dest.stat().st_size > 1000:
        print("exists", dest.name)
        return
    print("download", url)
    req = urllib.request.Request(url, headers={"User-Agent": "pdac-dd-e2e"})
    with urllib.request.urlopen(req, timeout=90) as resp:
        dest.write_bytes(resp.read())


def ca_map(struct, chain="A"):
    model = next(struct.get_models())
    if chain not in model:
        chain = list(model.child_dict)[0]
    out = {}
    for res in model[chain]:
        if res.id[0] == " " and "CA" in res:
            out[res.id[1]] = res["CA"]
    return out


def align_pdb(src, dst, refmap, parser):
    from Bio.PDB import PDBIO, Superimposer
    struct = parser.get_structure("x", str(src))
    mobile = ca_map(struct, "A")
    common = [i for i in sorted(set(refmap) & set(mobile)) if 5 <= i <= 165]
    if len(common) < 50:
        raise RuntimeError(f"{src.name}: only {len(common)} shared CA atoms with 7RPZ")
    sup = Superimposer()
    sup.set_atoms([refmap[i] for i in common], [mobile[i] for i in common])
    sup.apply(list(struct.get_atoms()))
    io = PDBIO()
    io.set_structure(struct)
    io.save(str(dst))
    return float(sup.rms), len(common)


def _xyz(ln):
    return (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))


def het_lines(pdb_text):
    """Keep one chain-A nucleotide and the Mg closest to its phosphates.

    Alternate conformers keep altloc blank or A, then the altloc column is cleared.
    7T47 stores GppCp (GCP) as altloc A and GDP as altloc B. The GDP copy is dropped
    so docking sees the GTP-like state.
    """
    dropped = set()
    candidates = []
    for ln in pdb_text.splitlines():
        if not ln.startswith("HETATM"):
            continue
        resn = ln[17:20].strip()
        chain = ln[21:22]
        alt = ln[16]
        if resn not in KEEP_HETS or chain not in {"A", " "}:
            dropped.add(resn)
            continue
        if alt not in {" ", "A"}:
            dropped.add(f"{resn}:alt{alt}")
            continue
        if alt != " ":
            ln = ln[:16] + " " + ln[17:]
        candidates.append(ln)
    nucleotides = [ln for ln in candidates if ln[17:20].strip() != "MG"]
    mgs = [ln for ln in candidates if ln[17:20].strip() == "MG"]
    if not nucleotides:
        return [], dropped
    names = {ln[17:20].strip() for ln in nucleotides}
    if "GDP" in names and len(names) > 1:
        nucleotides = [ln for ln in nucleotides if ln[17:20].strip() != "GDP"]
        dropped.add("GDP:secondary")
    res_ids = {(ln[21:22], ln[22:26]) for ln in nucleotides}
    if len(res_ids) > 1:
        biggest = max(res_ids, key=lambda key: sum(1 for ln in nucleotides if (ln[21:22], ln[22:26]) == key))
        nucleotides = [ln for ln in nucleotides if (ln[21:22], ln[22:26]) == biggest]
    phos = [_xyz(ln) for ln in nucleotides if ln[12:16].strip().startswith("P")]
    kept_mg = []
    if mgs and phos:
        def dist(ln):
            x, y, z = _xyz(ln)
            return min((x - px) ** 2 + (y - py) ** 2 + (z - pz) ** 2 for px, py, pz in phos)
        kept_mg = [min(mgs, key=dist)]
        if len(mgs) > 1:
            dropped.add(f"MG:extra{len(mgs) - 1}")
    elif mgs:
        kept_mg = [mgs[0]]
    # Crystal nucleotides sometimes carry extra protons (for example on the
    # phosphates). Meeko's GDP template is the heavy-atom residue, so those
    # hydrogens are dropped and Meeko supplies its own.
    def heavy(ln):
        elem = ln[76:78].strip().upper() if len(ln) >= 78 else ""
        name = ln[12:16].strip()
        return elem != "H" and not name.startswith("H")
    return [ln for ln in nucleotides + kept_mg if heavy(ln)], dropped


def protein_lines(pdb_text):
    """Chain A protein atoms. Prefer altloc blank, otherwise the first altloc."""
    chosen = {}
    order = []
    for ln in pdb_text.splitlines():
        if not ln.startswith("ATOM") or ln[21:22] != "A":
            continue
        alt = ln[16]
        key = (ln[22:26], ln[12:16])
        if key not in chosen:
            chosen[key] = ln
            order.append(key)
        elif alt == " ":
            chosen[key] = ln
    return [chosen[k] for k in order]


def write_receptor_pdb(protein_pdb_text, cofactors, dest: Path):
    body = [ln for ln in protein_pdb_text.splitlines() if ln.startswith(("ATOM", "HETATM", "TER"))]
    dest.write_text("\n".join(body + cofactors + ["END"]) + "\n")


def fix_protein(src: Path, dst: Path, mutation=None):
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=str(src))
    fixer.removeHeterogens(keepWater=False)
    if mutation:
        fixer.applyMutations([mutation], "A")
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    with open(dst, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)


def relax_sidechains(pdb_path: Path, cofactors=None, focus_residue=None):
    """Move a modeled side chain out of clashes with the protein and the nucleotide.

    PDBFixer builds the new side chain without a rotamer search. The backbone
    CA atoms stay fixed except at the mutated residue and its two neighbors,
    which proline in particular has to adjust. GDP/Mg are kept as steric walls
    so the new side chain does not fall into the nucleotide site.
    """
    from openmm.app import PDBFile, ForceField, Simulation, NoCutoff, Modeller
    from openmm import CustomExternalForce, Platform, VerletIntegrator
    pdb = PDBFile(str(pdb_path))
    model = Modeller(pdb.topology, pdb.positions)
    model.addHydrogens(pH=7.4)
    ff = ForceField("amber14-all.xml")
    system = ff.createSystem(model.topology, nonbondedMethod=NoCutoff, constraints=None)
    skip = set()
    if focus_residue:
        skip = {focus_residue - 1, focus_residue, focus_residue + 1}
    restraint = CustomExternalForce("krestr*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("krestr", 2000.0)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for atom in model.topology.atoms():
        if atom.name != "CA":
            continue
        if int(str(atom.residue.id)) in skip:
            continue
        pos = model.positions[atom.index]
        restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(restraint)
    walls = []
    for ln in cofactors or []:
        if not ln.startswith("HETATM"):
            continue
        elem = ln[76:78].strip().upper() if len(ln) >= 78 else ""
        if elem == "H":
            continue
        walls.append((float(ln[30:38]) * 0.1, float(ln[38:46]) * 0.1, float(ln[46:54]) * 0.1))
    if walls:
        terms, defs = [], []
        wall = CustomExternalForce("")
        wall.addGlobalParameter("kwall", 8000.0)
        wall.addGlobalParameter("r0", 0.30)
        for i, (x, y, z) in enumerate(walls):
            wall.addGlobalParameter(f"gx{i}", x)
            wall.addGlobalParameter(f"gy{i}", y)
            wall.addGlobalParameter(f"gz{i}", z)
            terms.append(f"step(r0-d{i})*(d{i}-r0)^2")
            defs.append(f"d{i}=sqrt((x-gx{i})^2+(y-gy{i})^2+(z-gz{i})^2)")
        wall.setEnergyFunction("kwall*(" + "+".join(terms) + "); " + "; ".join(defs))
        for atom in model.topology.atoms():
            if atom.element is None or atom.element.symbol == "H":
                continue
            if skip and int(str(atom.residue.id)) not in skip:
                continue
            wall.addParticle(atom.index)
        if wall.getNumParticles():
            system.addForce(wall)
    sim = Simulation(model.topology, system, VerletIntegrator(0.001), Platform.getPlatformByName("CPU"))
    sim.context.setPositions(model.positions)
    sim.minimizeEnergy(maxIterations=500)
    state = sim.context.getState(getPositions=True)
    tmp = pdb_path.with_suffix(".relaxed.pdb")
    with open(tmp, "w") as fh:
        PDBFile.writeFile(model.topology, state.getPositions(), fh, keepIds=True)
    heavy = []
    for ln in tmp.read_text().splitlines():
        if ln.startswith("ATOM"):
            elem = ln[76:78].strip().upper() if len(ln) >= 78 else ""
            if elem == "H" or ln[12:16].strip().startswith("H"):
                continue
        if ln.startswith(("ATOM", "TER", "END")):
            heavy.append(ln)
    pdb_path.write_text("\n".join(heavy) + "\n")
    tmp.unlink(missing_ok=True)


def meeko_prepare(pdb, prefix, center, size, meeko):
    templates = []
    last = ""
    for _ in range(8):
        cmd = [
            meeko, "-i", str(pdb), "-o", str(prefix), "-p", "-v",
            "--default_altloc", "A",
            "--box_center", str(center[0]), str(center[1]), str(center[2]),
            "--box_size", str(size[0]), str(size[1]), str(size[2]),
        ]
        for spec in templates:
            cmd.extend(["--set_template", spec])
        proc = subprocess.run(cmd, text=True, capture_output=True)
        pdbqt = Path(str(prefix) + ".pdbqt")
        last = (proc.stdout or "") + "\n" + (proc.stderr or "")
        if proc.returncode == 0 and pdbqt.exists() and pdbqt.stat().st_size >= 1000:
            if "GDP" not in pdbqt.read_text() and "MG" not in pdbqt.read_text():
                print("warning: no GDP/MG record found in", pdbqt.name)
            if templates:
                print("meeko templates:", ", ".join(templates))
            return
        added = False
        for key, names in re.findall(r"residue_key='([^']+)'.*?tied for fewest missing and excess H: ([A-Za-z ]+)", last):
            choice = "HIE" if "HIE" in names.split() else names.split()[0]
            spec = f"{key}={choice}"
            if spec not in templates:
                templates.append(spec)
                added = True
                print(f"histidine protonation tie at {key}; using {choice}")
        if not added:
            break
    log = Path(str(prefix) + ".meeko.log")
    log.write_text(last)
    raise RuntimeError(f"Meeko failed for {pdb.name}; see {log.name}")


def main():
    ensure_structure_python()
    import yaml
    from Bio.PDB import PDBParser

    def load_yaml(path):
        with open(path) as fh:
            return yaml.safe_load(fh)

    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()
    cfg = load_yaml(a.config)
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if not cfg.get("resistance", {}).get("enabled", True):
        (out / "panel_manifest.json").write_text(json.dumps({
            "disabled": True, "center": [], "size": [], "receptors": [],
        }))
        print("resistance panel disabled")
        return
    box = cfg["target"]["box"]
    center = [float(x) for x in box["center"]]
    size = [float(x) for x in box["size"]]
    meeko = cfg.get("resistance", {}).get("meeko_receptor_exe", "mk_prepare_receptor.py")
    if shutil.which(meeko) is None:
        meeko = env_tool("pdac-dd", meeko)
    if shutil.which(meeko) is None and not os.path.isfile(meeko):
        raise SystemExit(f"Missing Meeko receptor tool: {meeko}")

    raw = out / "raw"
    aligned = out / "aligned"
    receptors = out / "receptors"
    for p in (raw, aligned, receptors):
        p.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(a.panel), delimiter="\t"))
    experimental = [r for r in rows if not r["pdb_id"].startswith("MODEL_")]
    for r in experimental:
        download(f"https://files.rcsb.org/download/{r['pdb_id']}.pdb", raw / f"{r['pdb_id']}.pdb")

    parser = PDBParser(QUIET=True)
    ref = parser.get_structure("ref", str(raw / "7RPZ.pdb"))
    refmap = ca_map(ref, "A")
    metrics = []
    for r in experimental:
        pid = r["pdb_id"]
        rms, n = align_pdb(raw / f"{pid}.pdb", aligned / f"{pid}_aligned.pdb", refmap, parser)
        metrics.append({"id": pid, "rmsd_to_7RPZ_CA": rms, "n_ca": n, "source_type": "experimental"})
        print(f"aligned {pid} rmsd {rms:.3f} n_ca {n}")

    ref_text = (raw / "7RPZ.pdb").read_text()
    template_cofactors, dropped = het_lines(ref_text)
    print("7RPZ cofactors kept", len(template_cofactors), "dropped hets", sorted(dropped))
    for mid, mut in MUTATIONS.items():
        prot = aligned / f"{mid}_protein.pdb"
        fix_protein(raw / "7RPZ.pdb", prot, mutation=mut)
        write_receptor_pdb(prot.read_text(), template_cofactors, aligned / f"{mid}_aligned.pdb")
        prot.unlink(missing_ok=True)
        print("modeled", mid, mut)

    manifest_rows = []
    for r in rows:
        pid = r["pdb_id"]
        src = aligned / f"{pid}_aligned.pdb"
        text = src.read_text()
        cofactors, dropped = het_lines(text)
        prot = receptors / f"{pid}_protein.pdb"
        prot.write_text("\n".join(protein_lines(text) + ["END"]) + "\n")
        fix_protein(prot, receptors / f"{pid}_fixed.pdb")
        if pid.startswith("MODEL_"):
            focus = int(MUTATIONS[pid].split("-")[1])
            relax_sidechains(receptors / f"{pid}_fixed.pdb", cofactors, focus)
        receptor_pdb = receptors / f"{pid}.pdb"
        write_receptor_pdb((receptors / f"{pid}_fixed.pdb").read_text(), cofactors, receptor_pdb)
        prefix = receptors / pid
        meeko_prepare(receptor_pdb, prefix, center, size, meeko)
        manifest_rows.append({
            "pdb_id": pid,
            "pdbqt": f"receptors/{pid}.pdbqt",
            "label": r["label"],
            "state": r["state"],
            "role": r["role"],
            "mutation": r["mutation"],
            "source_type": r["source_type"],
            "notes": r["notes"],
            "dropped_hets": ",".join(sorted(dropped)),
        })
        kept_names = sorted({ln[17:20].strip() for ln in cofactors})
        print("prepared", pid, "kept", kept_names, "excluded", sorted(dropped))

    manifest = {
        "center": center,
        "size": size,
        "reference_pdb": "7RPZ",
        "note": "Modeled mutants are stress tests built from 7RPZ. They are not experimental structures.",
        "receptors": manifest_rows,
    }
    (out / "panel_manifest.json").write_text(json.dumps(manifest, indent=2))
    with open(out / "alignment_metrics.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "rmsd_to_7RPZ_CA", "n_ca", "source_type"])
        writer.writeheader()
        writer.writerows(metrics)
    print("panel ready", len(manifest_rows), "receptors")


if __name__ == "__main__":
    main()
