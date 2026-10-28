#!/usr/bin/env bash
# Create every conda environment and fetch every external code dependency.
#
#   scripts/bootstrap.sh            # everything
#   scripts/bootstrap.sh a|a2|b|c   # one track only
#   scripts/bootstrap.sh check      # report what is installed, change nothing
#
# Needs: conda or mamba (Miniforge recommended), git. Nextflow and CUDA are checked, not installed.
set -euo pipefail
source "$(dirname "$0")/lib.sh"

SOLVER="$(command -v mamba || command -v conda || true)"
WHAT="${1:-all}"

mk_env() {  # mk_env <env-name> <yml> [pip-extra...]
  local name="$1" yml="$2"
  if env_prefix "$name" >/dev/null; then log "env $name already present"; return; fi
  [[ -n "$SOLVER" ]] || die "conda/mamba not found. Install Miniforge first: https://github.com/conda-forge/miniforge"
  log "creating env $name from ${yml#"$REPO"/}"
  "$SOLVER" env create -y -n "$name" -f "$yml"
}

check() {
  log "environment report"
  for e in pdac-dd pdac-md admet-ai boltz2 amber openfe lbz retrochimera kras-adar; do
    if p="$(env_prefix "$e")"; then printf '  %-14s %s\n' "$e" "$p"; else printf '  %-14s MISSING\n' "$e"; fi
  done
  for t in nextflow git gh nvidia-smi; do
    if command -v "$t" >/dev/null; then printf '  %-14s %s\n' "$t" "$(command -v "$t")"; else printf '  %-14s MISSING\n' "$t"; fi
  done
  [[ -d "$REPO/external/rfd4-proteina" ]] && echo "  rfd4-proteina  $REPO/external/rfd4-proteina" || echo "  rfd4-proteina  MISSING"
}

track_a() {
  local A="$TRACKS/a_small_molecule_discovery"
  mk_env pdac-dd "$A/environment.yml"
  log "installing PyTorch (CUDA 12.8 wheels) and Meeko into pdac-dd"
  "$(env_prefix pdac-dd)/bin/pip" install --quiet torch --index-url https://download.pytorch.org/whl/cu128 || \
    warn "PyTorch CUDA wheel failed; install a build that matches your driver"
  mk_env pdac-md "$A/envs/md.yml"
  mk_env admet-ai "$A/envs/admet.yml"
  mk_env boltz2 "$A/envs/boltz2.yml"
  mk_env amber "$A/envs/amber.yml"
  mk_env openfe "$A/envs/openfe.yml"
  mk_env lbz "$A/envs/lbz.yml"
  mk_env retrochimera "$A/envs/retrochimera.yml"
  command -v nextflow >/dev/null || warn "Nextflow missing. Install: curl -s https://get.nextflow.io | bash && mv nextflow ~/bin/"
  warn "RetroChimera + Syntheseus + checkpoints are installed by hand: see docs/INSTALL.md (section RetroChimera)."
}

track_a2() { mk_env krasg12d_e2e "$TRACKS/a2_ensemble_docking_triage/environment.yml"; }

track_b() {
  mk_env kras-adar "$TRACKS/b_adar_rna_editing/environment-kras-adar.yml"
}

track_c() {
  local dest="$REPO/external/rfd4-proteina"
  local url="${RFD4_GIT_URL:-https://github.com/RosettaCommons/RFD4-Proteina.git}"
  mkdir -p "$REPO/external"
  if [[ -d "$dest/.git" ]]; then log "rfd4-proteina already cloned"; else
    log "cloning $url"
    git clone --depth 1 "$url" "$dest"
  fi
  log "RFD4-Proteina setup (creates .venv, downloads ~1.7 GB of released weights)"
  ( cd "$dest" && bash setup.sh --skip-search-envs )
  # prepare_targets.py needs gemmi, biotite and pyyaml inside that venv
  "$dest/.venv/bin/pip" install --quiet gemmi biotite pyyaml 2>/dev/null || \
    "$dest/.venv/bin/python" -m pip install --quiet gemmi biotite pyyaml || warn "install gemmi, biotite, pyyaml into the RFD4 venv by hand"
}

case "$WHAT" in
  check) check ;;
  a) track_a ;;
  a2) track_a2 ;;
  b) track_b ;;
  c) track_c ;;
  all) track_a; track_a2; track_b; track_c; check ;;
  *) die "usage: bootstrap.sh [all|a|a2|b|c|check]" ;;
esac
