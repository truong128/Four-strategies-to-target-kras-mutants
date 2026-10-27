#!/usr/bin/env bash
set -euo pipefail
PDB="$1"; PREFIX="$2"; CX="$3"; CY="$4"; CZ="$5"; SX="$6"; SY="$7"; SZ="$8"
mk_prepare_receptor.py -i "$PDB" -o "$PREFIX" -p -v \
  --box_center "$CX" "$CY" "$CZ" --box_size "$SX" "$SY" "$SZ"
