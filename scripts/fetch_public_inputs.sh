#!/usr/bin/env bash
# Download the PUBLIC inputs the pipelines need. Nothing is bundled in this repository.
#
#   scripts/fetch_public_inputs.sh [a|a2|b|c|all]
#
# a  : RCSB 7RPZ -> receptor (PDBFixer + Meeko), docking box, reference ligand
# a2 : nothing (01_download_structures.py fetches the PDB panel itself)
# b  : Ensembl release 113 cDNA + GTF, then KRAS CDS, gene map and ADAR1 FASTAs (UniProt P55265)
# c  : KRAS allele targets (7RPZ backbone, mutated residue grafted) via prepare_targets.py
#
# NOT fetched (license or size): Enamine REAL, ChEMBL/ZINC subsets, RESTORE 2.0 supplementary files,
# Boltz-2 weights (downloaded by Boltz on first use), RetroChimera checkpoints. See docs/DATA_ACCESS.md.
set -euo pipefail
source "$(dirname "$0")/lib.sh"
WHAT="${1:-all}"

a() {
  log "Track A: fetching 7RPZ and preparing the receptor"
  ( cd "$TRACKS/a_small_molecule_discovery" && "$(env_exe pdac-dd python)" scripts/fetch_structures.py )
}
b() {
  log "Track B: Ensembl reference + CDS + ADAR1 FASTAs"
  local B="$TRACKS/b_adar_rna_editing" py
  py="$(env_prefix kras-adar)/bin"
  ( cd "$B" && PATH="$py:$PATH" bash scripts/fetch_ensembl_reference.sh && PATH="$py:$PATH" ./run_all.sh setup-real )
}
c() {
  log "Track C: building the 18 KRAS allele targets"
  [[ -d "$REPO/external/rfd4-proteina/.venv" ]] || die "Run scripts/bootstrap.sh c first"
  a
  ( cd "$TRACKS/c_peptide_binders" && "$REPO/external/rfd4-proteina/.venv/bin/python" prepare_targets.py )
}
case "$WHAT" in
  a) a ;; a2) log "Track A2 downloads its own structures" ;; b) b ;; c) c ;;
  all) a; b; c ;;
  *) die "usage: fetch_public_inputs.sh [all|a|a2|b|c]" ;;
esac
