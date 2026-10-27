#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
# Use the kras-adar conda environment when the caller has not activated it.
KRAS_ENV="${KRAS_ADAR_ENV:-$HOME/miniforge3/envs/kras-adar}"
if [[ "${CONDA_DEFAULT_ENV:-}" != "kras-adar" && -x "$KRAS_ENV/bin/python3" ]]; then export PATH="$KRAS_ENV/bin:$PATH"; fi
mode="${1:-}"; shift || true
case "$mode" in
 smoke) python3 -m unittest discover -s tests -p "test_*.py" -v ;;
 discover) python3 scripts/adar_pipeline.py "$@" ;;
 advanced) python3 scripts/advanced.py "$@" ;;
 restore2) python3 scripts/restore2_reference.py "$@" ;;
 publisher) python3 scripts/ingest_publisher_workbooks.py "$@" ;;
 transfer) python3 scripts/transfer_publisher_chemistry.py "$@" ;;
 adaredit-fetch) bash scripts/get_adaredit.sh ;;
 external) bash scripts/external_tools.sh "$@" ;;
 audit) python3 scripts/workflow_audit.py "$@" ;;
 structural-audit) python3 scripts/structural_ml_preflight.py "$@" ;;
 validate-inputs) python3 scripts/validate_inputs.py "$@" ;;
 adaredit-import) python3 scripts/external_adaredit.py "$@" ;;
 references) python3 -m full.references "$@" ;;
 doctor) python3 -m full.workflow doctor "$@" ;;
 preflight) python3 -m full.workflow preflight "$@" ;;
 setup-real) python3 scripts/setup_real_run.py "$@" ;;
 full-e2e) python3 -m full.workflow run "$@" ;;
 align) python3 -m full.alignment "$@" ;;
 somatic) python3 -m full.somatic "$@" ;;
 protein-rna) python3 -m full.structure "$@" ;;
 forcefield) python3 -m full.parameterize "$@" ;;
 evidence) python3 -m full.evidence "$@" ;;
 dynamics) python3 -m full.dynamics "$@" ;;
 ml) python3 -m full.learning "$@" ;;
 rna-offtarget) python3 -m full.offtarget "$@" ;;
 rna-assay) python3 -m full.experimental "$@" ;;
 ranking) python3 -m full.prioritize "$@" ;;
 e2e) python3 scripts/adar_pipeline.py --out results/run "$@"; python3 scripts/advanced.py --input results/run --out results/run; if [[ -f "reference/publisher/publisher_constructs.csv" ]]; then python3 scripts/transfer_publisher_chemistry.py --reference-dir reference/publisher --guides results/run/02_guides.csv --out results/run/restore2_published; fi; if [[ -n "${RESTORE2_REF_CSV:-}" ]]; then python3 scripts/restore2_reference.py --reference-csv "$RESTORE2_REF_CSV" --guides-csv results/run/02_guides.csv --out results/run/restore2; fi; python3 scripts/workflow_audit.py --results results/run --out results/audit; python3 scripts/structural_ml_preflight.py ;;
 *) echo 'Usage: ./run_all.sh smoke | doctor | preflight --config FILE | setup-real | discover [args] | advanced [args] | restore2 [args] | e2e [discovery args]' >&2; exit 2 ;;
esac
