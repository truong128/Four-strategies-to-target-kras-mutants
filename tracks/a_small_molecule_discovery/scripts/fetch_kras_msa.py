#!/usr/bin/env python
"""Fetch one ColabFold MSA per KRAS allele and save Boltz csv files."""
import csv
import re
from pathlib import Path

from boltz.data.msa.mmseqs2 import run_mmseqs2

ROOT = Path(__file__).resolve().parents[1]
AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

OUT = ROOT / "results" / "kras_msa"
WT_PDB = ROOT / "results" / "allele_dock" / "receptors" / "6MBU.pdb"
PANEL = ROOT / "config" / "kras_allele_panel.tsv"
MAX_SEQS = 16384


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


def apply_mutation(seq, allele, mutation_field):
    if allele == "WT":
        return seq
    if mutation_field:
        wt3, resi, mut3 = mutation_field.split("-")
        pos, wt, mut = int(resi), AA3[wt3], AA3[mut3]
    else:
        match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", allele)
        pos, wt, mut = int(match.group(2)), match.group(1), match.group(3)
    if seq[pos - 1] != wt:
        raise SystemExit(f"{allele}: residue {pos} is {seq[pos - 1]}, expected {wt}")
    return seq[: pos - 1] + mut + seq[pos:]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    wt = "".join(chain_a_map(WT_PDB)[i] for i in range(1, 169))
    names, seqs = [], []
    with PANEL.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mutation = row["mutation"] if row["mutation"] else ""
            seq = apply_mutation(wt, row["allele"], mutation)
            names.append(row["allele"])
            seqs.append(seq)
    print(f"requesting MSA for {len(seqs)} sequences", flush=True)
    a3ms = run_mmseqs2(
        seqs,
        str(OUT / "mmseqs_tmp"),
        use_env=True,
        use_pairing=False,
        host_url="https://api.colabfold.com",
    )
    if len(a3ms) != len(names):
        raise SystemExit(f"Server returned {len(a3ms)} MSAs for {len(names)} sequences")
    for name, a3m in zip(names, a3ms):
        lines = a3m.strip().splitlines()
        sequences = lines[1::2][:MAX_SEQS]
        if not sequences:
            raise SystemExit(f"Empty MSA for {name}")
        body = ["key,sequence"] + [f"-1,{seq}" for seq in sequences]
        dest = OUT / f"{name}.csv"
        dest.write_text("\n".join(body) + "\n")
        print(f"wrote {dest.name} sequences {len(sequences)}", flush=True)


if __name__ == "__main__":
    main()
