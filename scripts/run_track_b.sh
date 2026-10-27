#!/usr/bin/env bash
# Track B: ADAR-recruiting guide RNA design for KRAS (G>A variants).
#
#   scripts/run_track_b.sh smoke                 # unit and regression tests (no data, no GPU)
#   scripts/run_track_b.sh setup                 # Ensembl reference, KRAS CDS, ADAR1 FASTAs
#   scripts/run_track_b.sh preflight             # checks every input before heavy steps
#   scripts/run_track_b.sh plan                  # design stages only (variant annotation ... ViennaRNA)
#   scripts/run_track_b.sh full                  # plus Bowtie2 off-target screen and Boltz-2 geometry (GPU)
#
# Every run writes to a NEW directory: edit `out:` in the config or pass --out. Existing outputs are refused.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
B="$TRACKS/b_adar_rna_editing"
CFG="${CFG:-config/real_run_pdac_kras.yaml}"
MODE="${1:-}"; shift || true
KENV="$(env_prefix kras-adar)" || die "env kras-adar missing. Run scripts/bootstrap.sh b"
export PATH="$KENV/bin:$PATH" KRAS_ADAR_ENV="$KENV"
cd "$B"
case "$MODE" in
  smoke)     ./run_all.sh smoke ;;
  setup)     bash scripts/fetch_ensembl_reference.sh && ./run_all.sh setup-real ;;
  preflight) ./run_all.sh preflight --config "$CFG" ;;
  plan)      ./run_all.sh full-e2e --config "$CFG" "$@" ;;
  full)      ./run_all.sh full-e2e --config "$CFG" --execute-heavy "$@" ;;
  *) die "usage: run_track_b.sh smoke | setup | preflight | plan | full" ;;
esac
