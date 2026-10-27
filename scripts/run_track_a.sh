#!/usr/bin/env bash
# Track A: ML-guided ultra-large docking campaign (Nextflow).
#
#   scripts/run_track_a.sh smoke                       # unit smoke test, no external tools
#   scripts/run_track_a.sh toy                         # 5-molecule toy library through library QC only
#   scripts/run_track_a.sh screen --library PATH       # full campaign; PATH = directory of Enamine
#                                                      #   *.cxsmiles.bz2 shards, or any .smi/.csv library
#   scripts/run_track_a.sh screen --library PATH --no-resume
#   scripts/run_track_a.sh alleles                     # 18-allele dock + Boltz-2 panel on the final leads
#
# Extra arguments after the mode are passed to `nextflow run` (for example -with-report).
set -euo pipefail
source "$(dirname "$0")/lib.sh"
A="$TRACKS/a_small_molecule_discovery"
MODE="${1:-}"; shift || true

DD_PY="$(env_exe pdac-dd python)"
DD_BIN="$(dirname "$DD_PY")"
export PDAC_DD_BIN="$DD_BIN"

case "$MODE" in
  smoke)
    log "Track A smoke test"
    ( cd "$A" && "$DD_PY" scripts/smoke_test.py )
    ;;
  toy)
    log "Track A toy library QC"
    ( cd "$A" && "$DD_PY" scripts/prepare_library.py --input example/toy_library.smi \
        --output /tmp/pdac_toy_standardized.smi --config config.yaml --qc /tmp/pdac_toy_qc.json && cat /tmp/pdac_toy_qc.json )
    ;;
  screen)
    LIB=""; RESUME="-resume"; PASS=()
    while [[ $# -gt 0 ]]; do
      case "$1" in
        --library) LIB="$2"; shift 2 ;;
        --no-resume) RESUME=""; shift ;;
        *) PASS+=("$1"); shift ;;
      esac
    done
    [[ -n "$LIB" ]] || die "--library PATH is required (see docs/DATA_ACCESS.md for Enamine REAL)"
    need nextflow "Install: curl -s https://get.nextflow.io | bash"
    "$DD_PY" "$REPO/scripts/configure_local.py"
    [[ -s "$A/structures/KRAS_G12D.pdbqt" ]] || "$REPO/scripts/fetch_public_inputs.sh" a
    read -r CX CY CZ SX SY SZ < <("$DD_PY" - <<PY
import yaml
b = yaml.safe_load(open("$A/config.yaml"))["target"]["box"]
print(*b["center"], *b["size"])
PY
)
    log "Launching Nextflow (box center $CX $CY $CZ, size $SX $SY $SZ)"
    ( cd "$A" && nextflow run main.nf -profile local $RESUME \
        --library "$LIB" \
        --config "$A/config.local.yaml" \
        --receptor_pdbqt "$A/structures/KRAS_G12D.pdbqt" \
        --receptor_pdb "$A/structures/7RPZ_KRAS_G12D_protein.pdb" \
        --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
        --size_x "$SX" --size_y "$SY" --size_z "$SZ" \
        "${PASS[@]}" )
    ;;
  alleles)
    log "Allele panel: dock leads on 18 alleles, then Boltz-2 affinity"
    "$DD_PY" "$REPO/scripts/configure_local.py"
    ( cd "$A"
      "$DD_PY" scripts/dock_kras_alleles.py --all --config config.local.yaml
      "$(env_exe boltz2 python)" scripts/fetch_kras_msa.py
      "$DD_PY" scripts/boltz_allele_panel.py --config config.local.yaml )
    ;;
  *) die "usage: run_track_a.sh smoke | toy | screen --library PATH | alleles" ;;
esac
