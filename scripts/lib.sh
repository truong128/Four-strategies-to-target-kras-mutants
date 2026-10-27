#!/usr/bin/env bash
# Shared helpers for the wrapper scripts. Source this file; do not execute it.
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRACKS="$REPO/tracks"
export REPO TRACKS

log()  { printf '\033[1;34m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*"; }
warn() { printf '\033[1;33m[warn]\033[0m %s\n' "$*" >&2; }
die()  { printf '\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

# Print the prefix of conda env $1 (searches every common envs directory), or return 1.
env_prefix() {
  local name="$1" base roots=()
  [[ -n "${CONDA_ENVS_PATH:-}" ]] && IFS=: read -r -a roots <<< "$CONDA_ENVS_PATH"
  if command -v conda >/dev/null 2>&1; then roots+=("$(conda info --base 2>/dev/null)/envs"); fi
  roots+=("$HOME/miniforge3/envs" "$HOME/mambaforge/envs" "$HOME/miniconda3/envs" "$HOME/anaconda3/envs" "$HOME/micromamba/envs")
  for base in "${roots[@]}"; do
    if [[ -d "$base/$name/bin" ]]; then echo "$base/$name"; return 0; fi
  done
  return 1
}

# Print <env>/bin/<exe> or die with a bootstrap hint.
env_exe() {
  local env="$1" exe="$2" prefix
  prefix="$(env_prefix "$env")" || die "conda env '$env' not found. Run: scripts/bootstrap.sh"
  [[ -x "$prefix/bin/$exe" ]] || die "$exe not found in env '$env'"
  echo "$prefix/bin/$exe"
}

need() { command -v "$1" >/dev/null 2>&1 || die "'$1' is required but not on PATH. $2"; }
