#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
export PYTHONPATH="$ROOT/scripts:${PYTHONPATH:-}"
MODE="${1:-full}"
[[ -f inputs/library.csv ]] || { mkdir -p inputs; cp inputs/library.template.csv inputs/library.csv; echo "Created inputs/library.csv from the empty template. Add your id,smiles rows."; }
python scripts/00_check_environment.py
python scripts/01_download_structures.py
python scripts/02_build_structure_panel.py
python scripts/03_prepare_ligands.py
python scripts/04_prepare_receptors.py
python scripts/05_dock_panel.py
python scripts/06_rank_hits.py
python scripts/07_admet_ai.py || true
python scripts/08_make_summary.py
python scripts/09_export_top_poses.py
printf '\nDONE. Key outputs:\n  results/tables/ranked_hits.csv\n  results/tables/final_candidates.csv\n  results/figures/selectivity_vs_resistance.png\n  results/SUMMARY.md\n'
