#!/usr/bin/env bash
set -euo pipefail
REC="$1"; LIG="$2"; CX="$3"; CY="$4"; CZ="$5"; SX="$6"; SY="$7"; SZ="$8"
mk_prepare_ligand.py -i "$LIG" -o reference_ligand.pdbqt
vina --receptor "$REC" --ligand reference_ligand.pdbqt \
  --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
  --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
  --exhaustiveness 32 --num_modes 20 --out reference_redocked.pdbqt
mk_export.py reference_redocked.pdbqt -o reference_redocked.sdf
echo "Now compute atom-mapped heavy-atom RMSD to the crystallographic pose; do not accept by score alone."
