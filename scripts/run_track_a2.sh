#!/usr/bin/env bash
# Track A2: lightweight ensemble-docking triage (scripts 00-09). Use for libraries up to ~1e5 molecules,
# for example the final leads of Track A or any vendor subset.
#
#   scripts/run_track_a2.sh                  # uses tracks/a2_ensemble_docking_triage/inputs/library.csv
#   scripts/run_track_a2.sh --library my.csv # copies your CSV (columns: id,smiles) into place first
set -euo pipefail
source "$(dirname "$0")/lib.sh"
T="$TRACKS/a2_ensemble_docking_triage"
if [[ "${1:-}" == "--library" ]]; then cp "$2" "$T/inputs/library.csv"; shift 2; fi
PY_ENV="$(env_prefix krasg12d_e2e)" || die "env krasg12d_e2e missing. Run scripts/bootstrap.sh a2"
( cd "$T" && PATH="$PY_ENV/bin:$PATH" bash run_all.sh "$@" )
