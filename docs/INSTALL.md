# Install

## What you need

| Item | Needed for | Notes |
|---|---|---|
| Linux x86_64 | everything | macOS can run parts of Track B and A2; GPU stages need Linux + NVIDIA |
| conda / mamba (Miniforge recommended) | all tracks | [github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge) |
| git, curl, gzip, awk | all | |
| Nextflow | Track A | `curl -s https://get.nextflow.io \| bash && mv nextflow ~/bin/` (Java 17 or newer) |
| NVIDIA GPU + recent driver | Boltz-2, MD, RetroChimera, Track C | CUDA 12.8 wheels are installed by default; Blackwell (CUDA 13) works with `use_kernels: false` / `--no_kernels` |
| About 32 CPU threads (the defaults) | Track A at full scale | lower `docking.workers`, `ml.inference_workers`, `library.sample_*` in `config.yaml` for smaller machines |
| Disk | Track A | about 330 GB for compressed Enamine REAL, plus room for results; Track B needs about 20 GB free (preflight checks) |

## One command

```bash
scripts/bootstrap.sh          # everything
scripts/bootstrap.sh a        # or one track: a | a2 | b | c
scripts/bootstrap.sh check    # report only
```

Environments created (names matter; the wrappers look them up):

| Env | File | Used by |
|---|---|---|
| `pdac-dd` | `tracks/a_small_molecule_discovery/environment.yml` | Track A (RDKit, Vina, Meeko, scikit-learn, PyTorch) |
| `pdac-md` | `tracks/a_small_molecule_discovery/envs/md.yml` | OpenMM, PDBFixer, OpenFF, AmberTools |
| `admet-ai` | `.../envs/admet.yml` | ADMET-AI (`admet_predict`) |
| `boltz2` | `.../envs/boltz2.yml` | Boltz-2 |
| `amber` | `.../envs/amber.yml` | MM/GBSA (falls back to `pdac-md`) |
| `openfe` | `.../envs/openfe.yml` | OpenFE absolute binding free energy, the alchemical check on MM/GBSA (off by default) |
| `lbz` | `.../envs/lbz.yml` | multi-core `lbzip2` for streaming Enamine shards |
| `retrochimera` | `.../envs/retrochimera.yml` | retrosynthesis (manual step below) |
| `krasg12d_e2e` | `tracks/a2_ensemble_docking_triage/environment.yml` | Track A2 |
| `kras-adar` | `tracks/b_adar_rna_editing/environment-kras-adar.yml` | Track B |
| RFD4 `.venv` | created by RFD4-Proteina `setup.sh` | Track C |

On an NVIDIA GPU, pin `cuda-version` in the `openfe` environment to the CUDA version your driver supports (`nvidia-smi`, top right); a newer one fails with `CUDA_ERROR_UNSUPPORTED_PTX_VERSION` when OpenMM compiles its kernels.

`scripts/configure_local.py` resolves `admet_predict`, `boltz`, `vina`, Meeko, AmberTools, `lbzip2` and every stage interpreter into `config.local.yaml`, so you do not need to activate environments by hand. Use `PDAC_DD_ENV`, `PDAC_MD_ENV`, `BOLTZ_ENV`, `AMBER_ENV`, `ADMET_ENV`, `RETROCHIMERA_ENV` or `LBZIP2_ENV` if your environments have other names.

## RetroChimera (manual)

RetroChimera and Syntheseus are installed from source and need checkpoints that are not redistributable here.

```bash
conda activate retrochimera
git clone https://github.com/microsoft/retrochimera.git external/retrochimera
git clone https://github.com/microsoft/syntheseus.git external/syntheseus
pip install -e external/syntheseus -e external/retrochimera     # follow upstream notes for torch-geometric wheels
# download the USPTO-50K checkpoint per upstream instructions, then
export RETROCHIMERA_CKPT=/path/to/retrochimera_checkpoints
scripts/configure_local.py
```

The purchasable building-block inventory (`data/purchasable/retrostar_inventory.smi` in `config.yaml`) comes from the Syntheseus Retro* benchmark set; see `docs/DATA_ACCESS.md`. To skip retrosynthesis, set `retrosynthesis.enabled: false`.

## Optional binaries

| Tool | When | How |
|---|---|---|
| GNINA 1.3 | 18-allele CNN rescoring (`run_track_a.sh alleles`) | download the release binary; `export GNINA_EXE=/path/to/gnina` |
| Ensembl VEP, GATK | Track B somatic calling from FASTQ/BAM | not required for guide design from a variant table; see `tracks/b_adar_rna_editing/README.md` |
| AlphaFold3 / OpenFold3 | Track C refolding metrics | `bash scripts/setup_search_envs.sh` inside RFD4-Proteina |

## Verify

```bash
scripts/run_all.sh smoke                    # Track A smoke test + Track B unit tests
scripts/run_track_a.sh toy                  # library QC on a toy file
scripts/fetch_public_inputs.sh a            # downloads 7RPZ and writes a docking-ready receptor
scripts/run_track_b.sh preflight            # after `run_track_b.sh setup`
```

## Troubleshooting

| Symptom | Fix |
|---|---|
| `conda env 'pdac-dd' not found` | run `scripts/bootstrap.sh a`; or create an env with that name |
| Meeko cannot match a residue template | the helper retries with histidine tie-break templates; check `structures/*.meeko.log` |
| Boltz aborts on Blackwell GPUs | keep `boltz.use_kernels: false`; wrappers pass `--no_kernels` |
| Nextflow tasks cannot find Vina or RDKit | `scripts/run_track_a.sh` exports `PDAC_DD_BIN`; do not launch `nextflow run` by hand without it |
| `INPUT_REQUIRED_NOT_EXECUTED` in Track B | a required input is missing by design; read `RUN_STATUS.json` |
| Track B refuses the output directory | choose a new `out:`; existing outputs are never reused |
