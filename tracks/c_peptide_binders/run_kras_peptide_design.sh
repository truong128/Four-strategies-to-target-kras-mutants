#!/usr/bin/env bash
# Design 8-16 residue peptides against 18 KRAS allele targets (WT + 17 mutants) with RFD4-Proteina.
#
# Alleles: G12D, WT, G12C, G12V, G12R, G12A, G13D, A146T, G12S, G12F,
# G13C, G13R, Q61H, Q61L, Q61R, Q61K, A146V, K117N.
# Codons 12, 13, and 61 use the P-loop / switch-II pocket. A146 and K117
# use the local nucleotide-site residues. Metrics are off because
# AlphaFold3 / OpenFold3 was not installed.
#
# Usage:
#   ./run_kras_peptide_design.sh                  # prepare targets, then design
#   ./run_kras_peptide_design.sh --prepare-only   # structures and YAML only
#   N=4 BATCH=2 ./run_kras_peptide_design.sh      # fewer designs per allele
#   ALLELES=G12D ./run_kras_peptide_design.sh     # one allele
#   SMOKE=1 ./run_kras_peptide_design.sh          # one design, no compile warm-up
set -euo pipefail

CAMPAIGN="$(cd "$(dirname "$0")" && pwd)"
# rfd4-proteina is an external dependency: git clone https://github.com/RosettaCommons/RFD4-Proteina.git
# (scripts/bootstrap.sh does this into external/). Override with RFD4_REPO=/path/to/clone.
REPO="${RFD4_REPO:-$CAMPAIGN/../../external/rfd4-proteina}"
PREPARE_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --prepare-only) PREPARE_ONLY=1 ;;
    -h|--help)
      sed -n '2,16p' "$0"
      exit 0 ;;
    *) echo "Unknown argument: $arg" >&2; exit 2 ;;
  esac
done

# shellcheck disable=SC1091
source "$REPO/.venv/bin/activate"
export DATA_PATH="${DATA_PATH:-$REPO/data_path}"
export PATH="$HOME/.local/bin:$PATH"

ALLELES="${ALLELES:-G12D,WT,G12C,G12V,G12R,G12A,G13D,A146T,G12S,G12F,G13C,G13R,Q61H,Q61L,Q61R,Q61K,A146V,K117N}"
PEPTIDE="${PEPTIDE:-8-16}"
N="${N:-8}"
BATCH="${BATCH:-4}"
OUT="${OUT:-$CAMPAIGN/runs/kras_peptides}"

python "$CAMPAIGN/prepare_targets.py" --alleles "$ALLELES" --peptide "$PEPTIDE"

if [[ "$PREPARE_ONLY" == 1 ]]; then
  echo "Targets and spec are ready. Re-run without --prepare-only to generate peptides."
  exit 0
fi

cd "$REPO"
args=(
  inputs="$CAMPAIGN/specs/kras_mutant_peptides.yaml"
  num_replicates="$N"
  batch_size="$BATCH"
  out_dir="$OUT"
  '~metrics'
)
if [[ "${SMOKE:-0}" == 1 ]]; then
  args+=(inference_experiment=quickstart compile_lengths=null num_replicates=1 batch_size=1)
fi

echo "Generating ${N} peptide complexes per allele into $OUT"
rfd4-generate "${args[@]}"
echo "Designs: $OUT/inference_samples/"
