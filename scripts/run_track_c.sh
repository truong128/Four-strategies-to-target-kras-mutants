#!/usr/bin/env bash
# Track C: de novo peptide binders for 18 KRAS allele targets (RFdiffusion3 / Proteina via RFD4-Proteina).
#
#   scripts/run_track_c.sh prepare            # build targets + design spec only (CPU)
#   scripts/run_track_c.sh smoke              # one design for one allele (GPU)
#   scripts/run_track_c.sh design             # N=8 designs per allele, 18 alleles (GPU)
#   N=4 BATCH=2 ALLELES=G12D,G12C scripts/run_track_c.sh design
set -euo pipefail
source "$(dirname "$0")/lib.sh"
export RFD4_REPO="${RFD4_REPO:-$REPO/external/rfd4-proteina}"
[[ -d "$RFD4_REPO/.venv" ]] || die "RFD4-Proteina not installed. Run scripts/bootstrap.sh c"
[[ -s "$TRACKS/a_small_molecule_discovery/structures/7RPZ_KRAS_G12D_protein.pdb" ]] || "$REPO/scripts/fetch_public_inputs.sh" a
cd "$TRACKS/c_peptide_binders"
case "${1:-}" in
  prepare) ./run_kras_peptide_design.sh --prepare-only ;;
  smoke)   SMOKE=1 ALLELES="${ALLELES:-G12D}" ./run_kras_peptide_design.sh ;;
  design)  ./run_kras_peptide_design.sh ;;
  *) die "usage: run_track_c.sh prepare | smoke | design" ;;
esac
