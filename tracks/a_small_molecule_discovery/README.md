> **Part of [two-strategies-to-target-kras-mutants](../../README.md).** Use the wrappers in `scripts/` at the repository root (`bootstrap.sh`, `fetch_public_inputs.sh`, `run_track_a.sh`) instead of the manual steps below. Run `python ../../scripts/configure_local.py` once to write `config.local.yaml`, and pass that file as `--config`.

# PDAC-DD-E2E

Runnable, modular end-to-end small-molecule discovery workflow for a PDAC target.
The discovery screen is **KRAS G12D** (PDB 7RPZ, switch-II inhibitor MRTX1133). The receptor keeps GDP and magnesium so docking searches the switch-II pocket. After clustering, the shortlist is docked again across a resistance panel: G12D GDP, G12D GTP-like, wild-type KRAS, experimental resistance structures, and modeled resistance mutants.

## Workflow

Target/QC → library standardization → representative training sample → Vina docking →
top-docking virtual-active labels → 20-classifier benchmark → scaffold-separated
train/validation/calibration/locked-test → single-model vs top-3 soft vote →
isotonic calibration + class-conditional split conformal prediction →
chunked ultra-large inference → re-docking → ADMET-AI → diversity clustering →
resistance panel (selectivity + resistance robustness) →
Boltz-2 → RetroChimera retrosynthesis → OpenMM MD → MM/GBSA →
Pareto + rank aggregation.

## Classifiers

Training compares 20 models: the 16 classifiers from the BRFSS benchmark
(`logistic_regression`, `multinom_reg`, `svm_rbf`, `svm_linear`, `decision_tree`,
`random_forest`, `gradient_boosting`, `extra_trees`, `adaboost`, `bag_tree`,
`bagging`, `knn`, `lda`, `qda`, `naive_bayes`, `mlp`) and four neural nets
(`ft_transformer`, `ft_transformer_lite`, `realmlp`, `tabular_resnet`).
Names and the classical hyperparameters follow that repository's `get_classifiers`.
See `MODEL_LIST.md`.

The workflow keeps these selection rules:
- model selection on validation, not the locked test;
- best-single vs top-3 soft-vote comparison;
- probability calibration after model selection;
- one locked test evaluation.

Chemical-specific changes:
- scaffold-separated partitions replace temporal partitions;
- SMOTE is intentionally not used on ECFP bits;
- enrichment/recall metrics are added because the active class is ~1%;
- slow benchmark models are capped; only scalable models can become production models.

## Install

```bash
conda env create -f environment.yml
conda activate pdac-dd
pip install torch --index-url https://download.pytorch.org/whl/cu128
python scripts/smoke_test.py
```

Optional environments:

```bash
conda env create -f envs/admet.yml
conda env create -f envs/md.yml
```

Install Nextflow separately, then check:

```bash
nextflow -version
```

## Required target inputs for a real run

Prepare once, after structural review:

- `receptor.pdbqt` for Vina
- `receptor.pdb` for GNINA/MD
- docking box center `(x,y,z)` and size `(x,y,z)`
- preferably a co-crystal/reference ligand for redocking QC

A Meeko helper is included:

```bash
scripts/prepare_receptor_meeko.sh receptor.pdb RSK2 10 20 30 22 22 22
```

## Input library

Whitespace-delimited, no header:

```text
SMILES ID
CCO CMPD000001
CCN CMPD000002
```

For the Enamine REAL download, pass the directory of `.cxsmiles.bz2` files and leave the files separate. The shipped settings were tuned for a 16-core/32-thread workstation and are: sample the library in parallel with `lbzip2` (4 files at once, 7 threads each, conda env `lbz`; rewrite the thread counts in `config.yaml` for your hardware), dock 200,000 molecules (exhaustiveness 4, 28 workers), score the full library with 26 scoring processes plus 2 readers (rows that already fail MW, sLogP, TPSA or rotatable-bond limits in the Enamine columns are skipped before RDKit), keep the top 50,000, and redock those at exhaustiveness 8 with 28 workers. Do not decompress the library into one text file.

## Example Nextflow launch

If you already have a docking-labeled training table:

```bash
nextflow run main.nf -profile local \
  --library /path/library.smi \
  --training_docking_csv /path/training_docking.csv
```

To have the workflow generate the training docking labels:

```bash
nextflow run main.nf -profile local \
  --library /path/to/enamine_cxsmiles_dir \
  --receptor_pdbqt structures/KRAS_G12D.pdbqt \
  --receptor_pdb structures/7RPZ_KRAS_G12D_protein.pdb \
  --center_x 1.714 --center_y 4.927 --center_z -23.164 \
  --size_x 24 --size_y 24 --size_z 24
```

## Key outputs

```text
results/library/standardized.smi
results/sample/training_sample.smi
results/docking/training_docking.csv
results/ml/model_comparison.csv
results/ml/model_bundle.joblib
results/ml/locked_test_metrics.json
results/screen/ml_hits.csv.gz
results/redock/redocked_hits.csv
results/admet/admet_predictions.csv
results/admet/admet_filtered.csv
results/cluster/diverse_hits.csv
results/resistance/resistance_scores.csv
results/resistance/resistance_long.csv
results/synthesis/synthesis_scores.csv
results/final/ranked_candidates.csv
results/final/pareto_candidates.csv
```

## QC gates before a full ultra-large screen

1. Redock a known co-crystal ligand and verify a chemically sensible pose.
2. Benchmark known inhibitors/decoys if available.
3. Require useful enrichment on the **locked scaffold test**.
4. Confirm the production classifier is fast enough for the library scale.
5. Inspect false negatives among known actives before screening billions.
6. Treat docking/ML/ADMET/MD/MMGBSA as prioritization, not proof of efficacy.

## External-tool notes

- **Vina + Meeko**: main high-throughput docking adapter. The ultra-large screen uses only the G12D GDP receptor. The resistance stage then docks at most 200 clustered compounds into the aligned panel in `config/pdb_panel.tsv` (7RPZ, 7T47, 6MBU, 9U8J, 9U95, 9U99, 9U97, plus mutants modeled from 7RPZ). A larger wild-type margin is better. A smaller worst-mutant penalty is better. Modeled mutants are stress tests, not crystal structures. This panel is the noncovalent switch-II track. It does not model a RAS(ON) cyclophilin-A tri-complex, and it does not dock bypass pathways such as EGFR, PI3K, or YAP.
- **ADMET-AI**: called through `admet_predict`.
- **GNINA** is not in the automatic run. Vina already docks the short list; GNINA would only rescore those poses.
- **Boltz-2**: protein-ligand complex and binding affinity for at most 200 clustered compounds, using the `boltz2` environment. In the discovery pipeline the protein is passed as a single sequence (`msa: empty`), so that step does not wait on an MSA server. A higher affinity probability ranks higher. The predicted affinity value is a log10(IC50, µM), so a more negative value ranks higher. An optional drug-repurposing extension, not part of this Nextflow run, scores a user-supplied table of approved small molecules against KRAS G12D with a ColabFold MSA (`scripts/boltz_g12d_fda.py`; set `FDA_DRUG_PKL` and/or `FDA_DRUG_CSV`). Keep the rows you care about (for example affinity probability above 0.5) in `results/g12d_fda_boltz_msa/g12d_passed_gt_0.5.csv`, then score them on the other 17 alleles with `scripts/boltz_fda_passed_alleles.py` (run `scripts/fetch_kras_msa.py` first). See `docs/STEP_BY_STEP.md`.
- **OpenMM**: 100 ns for the top 20 docked compounds, using the `pdac-md` environment. The stability score is higher when the ligand moves less.
- **MM/GBSA**: AmberTools single-frame energy for those same 20 compounds. A more negative value ranks higher.
- **Retrosynthesis**: RetroChimera (local USPTO-50K checkpoint) plus Retro* search against the purchasable building-block inventory, for the 50 best-docked clustered compounds. A solved route gets a higher score when it is shorter. The RDKit ease score is still recorded. This checkpoint can hallucinate reactions, so a chemist has to check any route before it is treated as a synthesis plan. Boltz-2 finishes before this search, and the 100 ns MD starts after it, so the three GPU steps do not run together.
