> **Part of [Four_therapeutic_strategies_to_target_Pancreatic_Cancer](../../README.md).** Use `scripts/run_track_a2.sh` from the repository root. Put your own library in `inputs/library.csv` (columns `id,smiles`; see `inputs/library.template.csv`). No compounds are shipped.

# KRAS G12D / PDAC resistance-robust E2E screen

This is a **single-command, reproducible first-pass computational screening pipeline** for the **noncovalent MRTX1133-like KRAS-G12D Switch-II-pocket track**.

It does **not** pretend conventional docking can model zoldonrasib-like cyclophilin-A ternary/covalent chemistry. That requires a separate ternary-complex/covalent workflow.

## What the single run does

1. Downloads the experimental PDB panel from RCSB.
2. Aligns all structures to 7RPZ.
3. Builds additional resistance-mutant stress-test models from 7RPZ.
4. Defines a common docking box from crystallographic MRTX1133 (CCD **6IC**).
5. Keeps GDP/GTP analog and Mg while removing the bound screening ligand.
6. Standardizes user SMILES, builds 3D conformers, and adds crystallographic MRTX1133 as a positive control.
7. Docks every ligand against:
   - G12D GDP state,
   - G12D active-like state,
   - WT KRAS counter-screen,
   - experimental resistance structures,
   - modeled resistance variants.
8. Calculates G12D affinity heuristic, WT selectivity margin, and resistance-robustness penalties.
9. Runs **ADMET-AI v2** on user compounds.
10. Produces ranked CSVs and a selectivity-vs-resistance plot.

## Install once

```bash
conda env create -f environment.yml
conda activate krasg12d_e2e
```

Apple Silicon is supported by the chosen CPU-first toolchain, subject to availability of the listed conda packages. If a package resolver fails, use a current Miniforge/Mambaforge installation and `mamba env create -f environment.yml`.

## Add your library

Edit `inputs/library.csv`:

```csv
id,smiles
cmpd_001,CCO...
cmpd_002,CCN...
```

Do not put billions of compounds directly into this local ensemble-docking stage. Prefilter ultra-large libraries first (physicochemical/ML/pharmacophore/faster docking), then feed a tractable set here.

## One command

```bash
./run_all.sh
```

Key outputs:

- `results/tables/docking_long.csv`
- `results/tables/ranked_hits.csv`
- `results/tables/final_candidates.csv`
- `results/tables/admet_ai.csv`
- `results/figures/selectivity_vs_resistance.png`
- `results/SUMMARY.md`

## Scientific interpretation

The final score is intentionally **rank-based**, because Vina scores are not physical binding free energies. The purpose is triage, not proof of affinity.

The hard requirements before a compound becomes a real lead remain:

**biochemical G12D binding/inhibition → WT KRAS/HRAS/NRAS selectivity → pERK suppression → G12D-selective cell killing → resistance-variant cells → 3D PDAC organoid/CAF context → PK/tumor exposure → in vivo efficacy/safety.**

### Why no fully automatic MD/MMGBSA is enabled by default
KRAS structures contain nucleotide and Mg, and production-quality MD requires validated parameters/protonation for the nucleotide, Mg coordination, protein, and each ligand. Automatically forcing a generic ligand force field onto GDP/GppCp/Mg would look “complete” but can be scientifically wrong. This package therefore makes the docking/selectivity/resistance/ADMET stages fully one-click and leaves MD/FEP as a **post-ranking validation stage** after the user selects a chemically small shortlist and validates topology/parameterization. This is deliberate, not an omission hidden from the user.

For top hits, recommended follow-up is explicit-solvent MD with a validated KRAS nucleotide/Mg setup and ligand parameters, then relative FEP where a congeneric series permits it; do not use MM/GBSA alone as the lead-selection criterion.
