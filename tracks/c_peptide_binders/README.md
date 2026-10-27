> **Part of [two-strategies-to-target-kras-mutants](../../README.md).** Use `scripts/run_track_c.sh` from the repository root.

# Track C: de novo peptide binders for 18 KRAS alleles

Designs 8 to 16 residue peptides against wild-type KRAS and 17 oncogenic or resistance-linked variants with
[RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina) (RFdiffusion3 / Proteina family). This track is a
design engine only: it proposes sequences and structures, it does not rank them against experiment.

## Alleles

G12D, WT, G12C, G12V, G12R, G12A, G13D, A146T, G12S, G12F, G13C, G13R, Q61H, Q61L, Q61R, Q61K, A146V, K117N.
Codons 12, 13 and 61 define the P-loop and switch-II hot-spots; A146 and K117 sit at the nucleotide site.

## Files

| File | Role |
|---|---|
| `prepare_targets.py` | Downloads the public KRAS structure, builds one mutant target per allele, writes the design spec |
| `specs/kras_mutant_peptides.yaml` | Design template (hot-spot residues, peptide length range) |
| `run_kras_peptide_design.sh` | Orchestrates prepare + design; env vars `ALLELES`, `N`, `BATCH`, `PEPTIDE`, `SMOKE` |
| `bin/apptainer` | Small shim so containerized tools resolve on hosts without a system install |

## Run

```bash
scripts/bootstrap.sh c                 # clones RFD4-Proteina into external/ and builds its venv
scripts/run_track_c.sh prepare         # CPU: targets and spec only
scripts/run_track_c.sh smoke           # GPU: one design, one allele
N=4 BATCH=2 ALLELES=G12D,G12C scripts/run_track_c.sh design
```

RFD4-Proteina is an external dependency, never vendored. Set `RFD4_REPO` to use an existing clone and
`RFD4_GIT_URL` to use a fork.

## Interpretation

Structure-confidence metrics (AlphaFold3 / OpenFold3 refolding) are disabled by default because those tools are not
installed by the bootstrap. Treat every output as a hypothesis for binding assays, not as a validated binder.
