#!/usr/bin/env python
"""Average MM/GBSA over finished OpenMM trajectories.

The score is the mean DELTA TOTAL from AmberTools MMPBSA.py (igb=5, 0.15 M salt)
on frames taken from the production DCD. The first 20 ns are discarded. One frame
is kept per nanosecond from 21 through 99 ns (79 frames). Topologies use amber14
for the protein and OpenFF 2.2.1 with NAGL AM1-BCC charges for the ligand.

This is generalized Born (MM/GBSA), not Poisson-Boltzmann MM/PBSA.
"""
import argparse, os, re, subprocess, sys
import pandas as pd

AMINO = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
from envtools import env_prefix, env_tool

AMBER_ENVS = ("binding_energy", "amber", "pdac-md")
CPP = os.environ.get("CPPTRAJ_EXE") or env_tool(AMBER_ENVS, "cpptraj")
MMPBSA = os.environ.get("MMPBSA_EXE") or env_tool(AMBER_ENVS, "MMPBSA.py")
GB_IN = """&general
  startframe=1, interval=1, verbose=1,
/
&gb
  igb=5, saltcon=0.150, radiopt=0,
/
"""


def build_prmtop(d, lig_sdf, final_pdb):
    from openmm.app import PDBFile, ForceField, PME
    from openmm import unit
    from openff.toolkit import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    import parmed as pmd

    pdb = PDBFile(final_pdb)
    lig = Molecule.from_file(lig_sdf)
    # OpenFF 2.2.1 expects AM1-BCC. AmberTools sqm is not in this environment.
    # NAGL supplies the AM1-BCC model that OpenFF ships with the force field.
    lig.assign_partial_charges("openff-gnn-am1bcc-1.0.0.pt")
    ff = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    gen = SMIRNOFFTemplateGenerator(molecules=[lig], forcefield="openff-2.2.1")
    ff.registerTemplateGenerator(gen.generator)
    # No constraints: constrained water bonds have no Amber type, and MM/GBSA
    # evaluates energies on fixed coordinates rather than integrating them.
    system = ff.createSystem(
        pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * unit.nanometer, constraints=None
    )
    n_top = pdb.topology.getNumAtoms()
    n_sys = system.getNumParticles()
    if n_top != n_sys:
        raise RuntimeError(f"particle count {n_sys} does not match topology {n_top}")
    parm = pmd.openmm.load_topology(pdb.topology, system, xyz=pdb.positions)
    # Solvent is a contiguous block after the protein and ligand in this trajectory.
    n_complex = sum(1 for a in parm.atoms if a.residue.name.strip() not in {"HOH", "NA", "CL"})
    parm.strip(":HOH,NA,CL")
    if len(parm.atoms) != n_complex:
        raise RuntimeError(f"dry topology has {len(parm.atoms)} atoms, expected {n_complex}")
    parm.box = None
    from parmed.tools import changeRadii
    changeRadii(parm, "mbondi2").execute()
    dest = os.path.join(d, "dry.prmtop")
    parm.save(dest, overwrite=True)
    open(os.path.join(d, "n_complex.txt"), "w").write(str(n_complex) + "\n")
    return dest, n_complex


def write_frames(d, dcd, final_pdb, n_complex, smoke):
    import mdtraj as md
    # Reporter interval is 10 ps. stride 100 keeps one frame per nanosecond.
    # Frame i in this strided trajectory is production time i ns (21 through 99 ns).
    traj = md.load(dcd, top=final_pdb, stride=100, atom_indices=list(range(n_complex)))
    molecules = traj.topology.find_molecules()
    if len(molecules) < 2:
        raise RuntimeError("expected a protein molecule and a ligand molecule")
    # The DCD stores the ligand in a different periodic image from the protein.
    traj = traj.image_molecules(anchor_molecules=[molecules[0]], other_molecules=molecules[1:])
    keep = []
    for i in range(traj.n_frames):
        if smoke and i in (21, 22):
            keep.append(i)
        elif (not smoke) and 21 <= i <= 99:
            keep.append(i)
    if not keep:
        raise RuntimeError(f"no frames selected from {traj.n_frames}")
    traj[keep].save(os.path.join(d, "dry.nc"))
    return len(keep)


def run(cmd, cwd, log_name):
    env = os.environ.copy()
    amberhome = os.environ.get("AMBERHOME") or env_prefix(AMBER_ENVS)
    if amberhome:
        env["AMBERHOME"] = amberhome
    q = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env)
    with open(os.path.join(cwd, log_name), "w") as fh:
        fh.write(q.stdout or "")
        fh.write("\n---- stderr ----\n")
        fh.write(q.stderr or "")
    if q.returncode:
        tail = (q.stderr or q.stdout or "")[-1500:]
        raise RuntimeError(f"{cmd[0]} failed ({q.returncode}): {tail}")
    return q


def one(work, md_dir, smoke):
    d = os.path.abspath(work)
    os.makedirs(d, exist_ok=True)
    final = os.path.join(md_dir, "md_final.pdb")
    dcd = os.path.join(md_dir, "md.dcd")
    lig = os.path.join(md_dir, "lig.sdf")
    if not (os.path.isfile(final) and os.path.isfile(dcd) and os.path.isfile(lig)):
        return float("nan"), float("nan"), 0, "missing_md"
    prmtop = os.path.join(d, "dry.prmtop")
    if not os.path.isfile(prmtop):
        prmtop, n_complex = build_prmtop(d, lig, final)
    else:
        n_complex = int(open(os.path.join(d, "n_complex.txt")).read())
    n_written = write_frames(d, dcd, final, n_complex, smoke)
    open(os.path.join(d, "receptor.cpptraj"), "w").write(
        "parm dry.prmtop\nparmstrip :UNK\nparmwrite out receptor.prmtop\nrun\nquit\n"
    )
    open(os.path.join(d, "ligand.cpptraj"), "w").write(
        "parm dry.prmtop\nparmstrip !(:UNK)\nparmwrite out ligand.prmtop\nrun\nquit\n"
    )
    run([CPP, "-i", "receptor.cpptraj"], d, "receptor.log")
    run([CPP, "-i", "ligand.cpptraj"], d, "ligand.log")
    open(os.path.join(d, "gb.in"), "w").write(GB_IN)
    run([
        MMPBSA, "-O", "-i", "gb.in",
        "-cp", "dry.prmtop", "-rp", "receptor.prmtop", "-lp", "ligand.prmtop",
        "-y", "dry.nc", "-o", "FINAL_MMPBSA.dat", "-eo", "per_frame.csv",
    ], d, "mmpbsa.log")
    text = open(os.path.join(d, "FINAL_MMPBSA.dat")).read()
    hits = re.findall(r"DELTA TOTAL\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)", text)
    if not hits:
        return float("nan"), float("nan"), 0, "parse_fail"
    mean, sd = float(hits[-1][0]), float(hits[-1][1])
    nmatch = re.search(r"Calculations performed using\s+([0-9.]+)\s+complex frames", text)
    n = int(float(nmatch.group(1))) if nmatch else 0
    return mean, sd, n, "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md-root", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    scores = pd.read_csv(a.scores)
    if a.only:
        scores = scores[scores["id"] == a.only]
        if scores.empty:
            raise SystemExit(f"id not in scores: {a.only}")
    rows = []
    for r in scores.itertuples(index=False):
        tag = str(r.id).replace("/", "_")
        md_dir = os.path.join(a.md_root, tag)
        print(f"start {tag} smoke={a.smoke}", flush=True)
        try:
            mean, sd, n, status = one(os.path.join(a.workdir, tag), md_dir, a.smoke)
        except Exception as exc:
            mean, sd, n, status = float("nan"), float("nan"), 0, f"fail: {exc}"
            print(status, flush=True)
        rows.append({
            "id": r.id, "smiles": r.smiles,
            "mmgbsa_mean": mean, "mmgbsa_sd": sd, "n_frames": n, "mmgbsa_status": status,
        })
        print(f"done {tag} {status} mean={mean} sd={sd} n={n}", flush=True)
        pd.DataFrame(rows).to_csv(a.output, index=False)
    ok = sum(1 for r in rows if r["mmgbsa_status"] == "ok")
    print(f"mmgbsa_ok {ok}/{len(rows)}", flush=True)


if __name__ == "__main__":
    main()
