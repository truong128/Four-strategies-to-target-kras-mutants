# Step by step

Every step of every track, in the order you run it. Commands assume you are in the repository root. "Gate" lines say what to check before moving on. Output paths are relative to the track directory.

Conventions: `A` = `tracks/a_small_molecule_discovery`, `A2` = `tracks/a2_ensemble_docking_triage`, `B` = `tracks/b_adar_rna_editing`, `C` = `tracks/c_peptide_binders`.

---

## Step 0. Common setup (once)

| # | Command | What it does | Gate |
|---|---|---|---|
| 0.1 | `scripts/bootstrap.sh check` | lists which environments and tools exist | read the report |
| 0.2 | `scripts/bootstrap.sh` | creates `pdac-dd`, `pdac-md`, `admet-ai`, `boltz2`, `amber`, `lbz`, `retrochimera`, `krasg12d_e2e`, `kras-adar`; clones RFD4-Proteina into `external/` | `bootstrap.sh check` shows no MISSING you need |
| 0.3 | install RetroChimera by hand | see `docs/INSTALL.md`, section RetroChimera | `python -c "import retrochimera"` in the env |
| 0.4 | `scripts/fetch_public_inputs.sh` | downloads public inputs for tracks A, B, C | no error; `tracks/*/structures` and `tracks/b_adar_rna_editing/reference` exist |
| 0.5 | `scripts/run_all.sh smoke` | Track A smoke test + 31 Track B regression tests | both report PASS / OK |

---

## Track A. Ultra-large small-molecule discovery (Nextflow)

Entry point: `scripts/run_track_a.sh screen --library PATH`. Behind it: `nextflow run tracks/a_small_molecule_discovery/main.nf`.

### A.1 Target preparation

`python scripts/fetch_structures.py` (run through `scripts/fetch_public_inputs.sh a`)

1. Download PDB 7RPZ (KRAS G12D with the switch-II inhibitor MRTX1133, ligand code 6IC).
2. Keep chain A, one altloc per atom; remove the inhibitor and waters; keep GDP and Mg2+.
3. PDBFixer adds missing heavy atoms (no loop rebuilding).
4. Meeko writes `structures/KRAS_G12D.pdbqt` for the box in `config.yaml` (center 1.714, 4.927, -23.164; size 24 A).
5. Writes the reference ligand (`7RPZ_MRTX1133.pdb/.sdf`) for redocking QC.

**Gate (recommended):** `scripts/redocking_qc.sh` redocks the crystal ligand; compute an atom-mapped heavy-atom RMSD to the crystal pose. Do not accept a receptor by score alone.

### A.2 Library access

Obtain Enamine REAL as `.cxsmiles.bz2` shards (see `docs/DATA_ACCESS.md`). Keep them compressed in one directory. The pipeline streams them with `lbzip2`; never decompress to one file.

Non-streaming libraries (`.smi`, `.csv`) are standardized first (`prepare_library.py`, process `STANDARDIZE_LIBRARY`). Try `scripts/run_track_a.sh toy` to see the QC JSON on a 20-line toy file.

### A.3 Configure for this machine

`python scripts/configure_local.py` writes `tracks/a_small_molecule_discovery/config.local.yaml` with absolute interpreter paths for each conda environment (git-ignored). Then edit the **hardware knobs** in `config.yaml` (workers, threads, chunk sizes) to match your CPU; the defaults target a 16-core / 32-thread workstation.

### A.4 Launch

```bash
scripts/run_track_a.sh screen --library /path/to/enamine_real          # resumable
scripts/run_track_a.sh screen --library /path/to/enamine_real --no-resume
```

The Nextflow processes run in this order (results are published under `results/<folder>/`):

| # | Process | Script | Folder | What happens | Key output |
|---|---|---|---|---|---|
| 1 | `STANDARDIZE_LIBRARY` | `prepare_library.py` | `library/` | only for non-streamed input: SMILES standardization plus charge and property filters (counts in the QC JSON) | `standardized.smi`, `library_qc.json` |
| 2 | `SAMPLE_LIBRARY` | `sample_library.py` | `sample/` | Bernoulli or reservoir sampling of 200,000 molecules; Enamine columns pre-filter MW, sLogP, TPSA, RotB before RDKit | `training_sample.smi` |
| 3 | `DOCK_TRAINING` | `dock_vina_batch.py` | `docking/` | Vina, exhaustiveness 4, 28 workers; skipped if you pass `--training_docking_csv` | `training_docking.csv` |
| 4 | `TRAIN_ML` | `benchmark_16_models.py` | `ml/` | label top 1% as active; scaffold split 0.65/0.15/0.10/0.10; fit 20 classifiers; select on validation; calibrate; conformal; one locked test | `model_comparison.csv`, `model_bundle.joblib`, `locked_test_metrics.json`, `split_manifest.csv` |
| 5 | `SCREEN_LIBRARY` | `predict_ultralarge.py` | `screen/` | 26 scoring processes + 2 readers over the streamed library; keep top 50,000 that pass the conformal gate | `ml_hits.csv.gz` |
| 6 | `REDOCK` | `dock_vina_batch.py --redock` | `redock/` | Vina, exhaustiveness 8 on the 50,000 | `redocked_hits.csv` |
| 7 | `ADMET` | `run_admet.py` | `admet/` | ADMET-AI; hard filters hERG, AMES, hepatotoxicity probability < 0.5 | `admet_predictions.csv`, `admet_filtered.csv` |
| 8 | `CLUSTER` | `cluster_hits.py` | `cluster/` | Butina clustering, Tanimoto 0.55, ECFP4; one representative per cluster | `diverse_hits.csv` |
| 9 | `PREPARE_RESISTANCE_PANEL` | `build_resistance_panel.py` | `resistance/panel/` | download 7 experimental receptors, align to 7RPZ, model 9 variants, Meeko each | panel manifest, PDBQTs |
| 10 | `DOCK_RESISTANCE` | `dock_resistance_panel.py` | `resistance/` | dock up to 200 compounds into all receptors; G12D affinity, WT margin, worst resistance penalty | `resistance_scores.csv`, `resistance_long.csv` |
| 11 | `BOLTZ` | `run_boltz_batch.py` | `boltz/` | Boltz-2 complex + affinity, up to 200 compounds, protein as single sequence | `boltz_scores.csv` |
| 12 | `SYNTHESIS` | `synthetic_feasibility.py` | `synthesis/` | RetroChimera + Retro* on 50 compounds against the inventory; RDKit ease score recorded | `synthesis_scores.csv` |
| 13 | `MD` | `run_md_batch.py` | `md/` | 100 ns OpenMM (amber14 + OpenFF 2.2.1) on 20 compounds; ligand RMSD stability | `md_scores.csv` |
| 14 | `MMPBSA` | `run_mmpbsa_batch.py` | `mmpbsa/` | AmberTools GB (igb=5, 0.15 M) on the same 20 | `mmpbsa_scores.csv` |
| 15 | `FINAL_RANK` | `final_rank.py` | `final/` | percentile-rank each criterion, weighted sum (see T3), Pareto front, report 50 | `ranked_candidates.csv`, `pareto_candidates.csv` |

GPU steps (Boltz, retrosynthesis, MD) are chained so they never run together. Nextflow also writes `timeline.html`, `report.html`, `trace.txt`, `dag.html` into `results/`.

**Gates:**

1. `ml/locked_test_metrics.json` shows useful enrichment on the locked scaffold test before you trust step 5.
2. `ml/model_comparison.csv`: the production model must be in the `production_eligible` list.
3. Spot-check false negatives among known actives (for example MRTX1133 analogs) before screening billions.
4. Retrosynthesis routes are model output; a chemist must check every route.

### A.5 Extension: 18-allele panel on the final leads

```bash
scripts/run_track_a.sh alleles
```

1. `dock_kras_alleles.py --all` (`--prepare` builds receptors only, `--dock` docks only) aligns one GDP monomer per allele (`config/kras_allele_panel.tsv`) onto 7RPZ, docks the leads in the same box, and rescores with GNINA 1.3 CNN scoring (required for this extension; set `GNINA_EXE` or put `gnina` on PATH).
2. `fetch_kras_msa.py` fetches one ColabFold MSA per allele (**sends the KRAS sequence to the public ColabFold server**; skip if that is not acceptable and supply your own A3M).
3. `boltz_allele_panel.py` runs Boltz-2 affinity for every lead on every allele (leads x 18 jobs).

### A.6 Extension: repurposing screen of approved drugs

```bash
export FDA_DRUG_PKL=/path/to/approved_drugs.pkl     # columns: drug_id, drug_name, smile
export FDA_DRUG_CSV=/path/to/extra_drugs.csv        # optional; columns: name, smiles
cd tracks/a_small_molecule_discovery
python scripts/dock_kras_alleles.py --prepare       # builds results/allele_dock/receptors/6MBU.pdb and the allele receptors
python scripts/fetch_kras_msa.py                    # in the boltz2 env; contacts the ColabFold server
python scripts/boltz_g12d_fda.py                    # writes results/g12d_fda_boltz_msa/g12d_fda_boltz.csv
# keep the rows you care about, for example boltz_affinity_probability > 0.5
python - <<'PY'
import pandas as pd
t = pd.read_csv("results/g12d_fda_boltz_msa/g12d_fda_boltz.csv")
t[t.boltz_affinity_probability > 0.5].to_csv("results/g12d_fda_boltz_msa/g12d_passed_gt_0.5.csv", index=False)
PY
python scripts/boltz_fda_passed_alleles.py          # kept drugs x 17 other alleles
```

---

## Track A2. Lightweight ensemble-docking triage

Use it for a tractable library (up to about 1e5 molecules): the final leads of Track A, a vendor subset, or your own designs.

```bash
scripts/run_track_a2.sh --library my_library.csv     # columns: id,smiles
```

Scripts, in order (`tracks/a2_ensemble_docking_triage/scripts`):

| # | Script | What happens |
|---|---|---|
| 00 | `00_check_environment.py` | verifies obabel, Vina and Python modules |
| 01 | `01_download_structures.py` | downloads the PDB panel from RCSB |
| 02 | `02_build_structure_panel.py` | aligns to 7RPZ; builds resistance-mutant stress models; derives the box from crystallographic 6IC; keeps GDP/GTP analog and Mg |
| 03 | `03_prepare_ligands.py` | standardizes your SMILES, 3D conformers, adds MRTX1133 as positive control |
| 04 | `04_prepare_receptors.py` | receptor PDBQTs |
| 05 | `05_dock_panel.py` | docks every ligand into every receptor |
| 06 | `06_rank_hits.py` | G12D affinity, WT margin, resistance penalties; rank-based score |
| 07 | `07_admet_ai.py` | ADMET-AI v2 (optional, failure tolerated) |
| 08 | `08_make_summary.py` | tables, selectivity-versus-resistance plot, `SUMMARY.md` |
| 09 | `09_export_top_poses.py` | exports top poses |

MD and MM/GBSA are intentionally not automatic here; see `tracks/a2_ensemble_docking_triage/docs/METHODS_AND_CAVEATS.md`.

---

## Track B. ADAR-recruiting guide RNAs

Entry point: `scripts/run_track_b.sh <mode>`. Behind it: `tracks/b_adar_rna_editing/run_all.sh`.

### B.1 Reference setup

```bash
scripts/run_track_b.sh setup
```

1. `scripts/fetch_ensembl_reference.sh` downloads Ensembl 113 cDNA and GTF, strips FASTA headers to bare `ENST` ids, indexes with samtools.
2. `scripts/setup_real_run.py` (via `./run_all.sh setup-real`) cuts the KRAS-4B CDS (ENST00000311936) from the transcriptome using the GTF CDS and stop-codon features, validates it (ATG start, one stop, human KRAS N-terminus, codon 12 = Gly), builds the transcript-to-gene map, and downloads UniProt P55265 to write p150, p110 and the 716-1226 construct FASTAs with catalytic-residue checks (H910, E912, C966, C1036).

### B.2 Variants and preflight

Edit `config/kras_pdac_variants.csv` (`variant_id,hgvs_c`). Then:

```bash
scripts/run_track_b.sh preflight
```

17 checks: fresh output directory, free disk, variant table, CDS sanity, **every HGVS reference base matches the CDS**, transcriptome and gene map, Bowtie2, ADAR FASTA, Boltz and weights, GPU, RNAfold/RNAcofold. **Gate:** zero FAIL.

### B.3 Design and screening

```bash
scripts/run_track_b.sh plan      # design stages only, CPU
scripts/run_track_b.sh full      # + Bowtie2 off-target + Boltz-2 (GPU)
```

Always point `out:` in the config at a new directory; non-empty directories are refused.

| # | Stage | What happens | Output (under `out/`) |
|---|---|---|---|
| 1 | variant annotation | each HGVS change checked against the CDS; A-to-I reversibility call | `discovery/` |
| 2 | guide design | guides of 33, 40 and 59 nt with an orphan C opposite the target A | `discovery/02_guides.csv` |
| 3 | bystanders | every other A in the footprint, with codon consequence | `discovery/` |
| 4 | chemistry transfer | position-specific RESTORE 2.0 patterns mapped onto guides with provenance keys (needs the publisher reference below) | `restore2_published/` |
| 5 | local structure | ViennaRNA fold and cofold of guide and target | `discovery/` |
| 6 | off-target | Bowtie2 local alignment of each guide against the sense transcriptome; hits tagged `gene_name` and `is_on_target_gene`; KRAS isoforms are not counted as off-target | `offtarget/` |
| 7 | structure | Boltz-2 YAML for ADAR1 716-1226 + RNA duplex; ensemble samples reported separately | `boltz_input/`, `boltz_models/` |
| 8 | research ranking | experiment-gated; no therapeutic ranking without experimental input | `ranking/` |

`RUN_STATUS.json` is the stage ledger (`COMPLETED`, `PLANNED_NOT_RUN`, `SKIPPED_INPUT_REQUIRED`, `INPUT_REQUIRED_NOT_EXECUTED`, `FAILED`). Completed means the executor finished, not that the biology is valid.

### B.4 RESTORE 2.0 chemistry reference (optional, user supplied)

1. Download the Supplementary Data and Source Data workbooks of the RESTORE 2.0 article (DOI 10.1038/s41467-025-64434-7) from the publisher into `tracks/b_adar_rna_editing/source_material/` (git-ignored). `scripts/restore2_download.py --url URL --output FILE --expect xlsx` fetches a URL you supply and verifies its signature and hash.
2. `./run_all.sh publisher --supplement source_material/<supplement>.xlsx --source-data source_material/<source_data>.xlsx --out reference/publisher`.
3. Re-run `plan`; the transfer stage now fills `restore2_published/`.

Transferred chemistry is **not validated on KRAS**.

### B.5 Stages that need your own validated inputs

| Stage | You must supply | Template |
|---|---|---|
| somatic calling | tumor/normal FASTQ or BAM, GRCh38, known sites, VEP cache | `config/full_project.yaml` |
| force field + MD | validated Amber `.lib/.off/.frcmod`, stereochemistry, validation JSON | `examples/input_templates/` |
| chemistry-aware ML | curated measured editing + chemistry positions + independent target groups | `examples/input_templates/measured_*.csv` |
| experimental read counts | treated and control BAM, genome FASTA, site list | `examples/input_templates/rna_editing_sites.csv` |

Without them each stage abstains and says why.

---

## Track C. Peptide binders (RFD4-Proteina)

```bash
scripts/bootstrap.sh c                      # clone + install RFD4-Proteina, add gemmi/biotite/pyyaml to its venv
scripts/run_track_c.sh prepare              # CPU: build targets and the design spec
scripts/run_track_c.sh smoke                # GPU: one design, one allele
scripts/run_track_c.sh design               # GPU: 8 designs per allele x 18 alleles
N=4 BATCH=2 ALLELES=G12D,G12C scripts/run_track_c.sh design
```

1. `prepare_targets.py` takes 7RPZ chain A, drops residue 1 (disconnected because residue 2 has no coordinates), keeps G12D as in the crystal, restores Gly at 12 for alleles away from codon 12, grafts an ideal side chain at the mutated residue, and writes `targets/KRAS_<allele>.cif` plus `specs/kras_mutant_peptides.yaml` (hotspots, H-bond condition, contig `8-16, /0, A3-169`).
2. `run_kras_peptide_design.sh` calls `rfd4-generate` with that spec. Metrics are off by default; enable AlphaFold3 / OpenFold3 refolding by installing the oracle and removing the `~metrics` flag. `bin/apptainer` is a Docker stand-in for hosts that cannot run rootless Apptainer.

Outputs land in `tracks/c_peptide_binders/runs/` (git-ignored).

---

## Step Z. Close the loop in the lab

Everything above produces hypotheses. See Figure 8 in the README for the experiments that turn them into evidence, and `docs/SCIENTIFIC_CAVEATS.md` for what each computed number can and cannot support.
