> **Part of [two-strategies-to-target-kras-mutants](../../README.md).** Use `scripts/run_track_b.sh` from the repository root. Reference files (Ensembl GTF and transcriptome) are downloaded by `scripts/fetch_public_inputs.sh b` and are not shipped. The RESTORE 2.0 supplementary material is not redistributed; obtain it from the publication (see `docs/DATA_ACCESS.md` at the repository root). The notes below describe the software audit history of this track.

# v1.4 Adversarial Audit Release (2026-10-08)

**Scope:** source review, automated regressions, chemistry linkage invariants, and output isolation. This is a software audit, not experimental or GPU/MD verification.

**Fixes compared with v1.3**

1. Bowtie2 index is now derived under the unique per-run output directory. A configured index outside that directory is rejected to prevent stale or shared index reuse.
2. `./run_all.sh smoke` invokes all `test_*.py` unit and regression modules, including the release-audit tests previously omitted.
3. Chemistry-aware ML requires every internal internucleotide linkage to be a recognized PO/PS state; unknown or missing values cause explicit rejection.
4. The terminal nucleotide has no 3-prime internucleotide bond and is excluded from backbone-linkage feature denominators.
5. New adversarial tests exercise these invariants. 21 Python unit/regression tests passed on the package source.

**Unverified externally:** GATK/Ensembl VEP execution on actual tumor-normal data, Boltz GPU folding/complex inference, parameterized modified MD trajectories, calibrated KRAS-specific editing ML and biological editing/tumor responses. No claim of guaranteed correctness or biological efficacy is made.

# v1.2 strict reproducibility audit (2026-10-08)

This is a targeted bug-fix release based on a direct source-code audit. **Always select a fresh `out` path for each run.** Existing output directories containing files are rejected, preventing stale VCF, off-target, or fitted ML artifacts from a prior run from being interpreted as newly computed. `chemistry_model.joblib` is accepted only if the current run's `ml_report.json` explicitly marks a qualifying model. The positional chemistry feature encoder rejects a fictitious 3-prime linkage after the terminal nucleotide. One legacy synthetic unit fixture was corrected to have a proper terminal linkage. A new regression test checks stale-output rejection.

None of these engineering fixes establish KRAS G12D editing efficiency or validate force fields. The implementation still requires independent biological validation and verified modified chemical structures.

# KRAS G12D ADAR RNA Correction: E2E Computational Research Package

**Version 1.1, final source-code audit, 2026-10-08**

This research package is a reproducible computational workflow for prioritizing ADAR-recruiting antisense oligonucleotides against verified KRAS variants in pancreatic ductal adenocarcinoma (PDAC). It is not evidence of therapeutic activity, an ADAR editing-rate predictor validated on KRAS, or a turnkey synthesis protocol. The user supplies biological input data, external model weights, validated chemical parameters and eventually experimental readouts. No patient sample or experimental measurement is bundled.

## 1. Before you run

Run from the extracted package root. First test on the bundled artificial dataset:

```bash
./run_all.sh smoke
./run_all.sh full-e2e --config config/synthetic_integration.yaml
```

To run with real input, edit `config/full_project.yaml` and run:

```bash
./run_all.sh doctor
./run_all.sh full-e2e --config config/full_project.yaml
./run_all.sh full-e2e --config config/full_project.yaml --execute-heavy
```

**The default real-project config has no variant data.** Running it without data writes `RUN_STATUS.json` with `INPUT_REQUIRED_NOT_EXECUTED`. The demo files `config/variants.csv` and `config/synthetic_for_integration.csv` are explicitly rejected unless `discovery.allow_synthetic_panel: true` is set for testing. Do not enable that switch for patient projects.

`--execute-heavy` is an instruction to invoke installed external tools, not a substitute for correctly supplied input files. Run external steps individually first. Linux with CUDA hardware is recommended for Boltz-2. Other modules run on macOS or Linux with the corresponding dependencies.

## 2. Install strategy

The core sequence and basic design tools can run with Python 3.11, PyYAML and the bundled Python source. For a practical CPU environment:

```bash
micromamba create -n kras-adar-core -c conda-forge -c bioconda python=3.11 pyyaml numpy pandas scikit-learn joblib biopython viennarna bowtie2 samtools bcftools -y
micromamba activate kras-adar-core
```

For advanced chemistry and simulation, see `environment-full.yml` and install the official tools separately where necessary. In particular, setup for Boltz-2 (models, MSA and GPUs), GATK, Ensembl VEP, AmberTools/tleap, custom residue parameters, OpenMM and MDTraj must be verified on the destination computer. An optional `artifact-tool` dependency is needed for the publisher workbook importer. Run `./run_all.sh doctor` to see installed dependencies. Environment installation alone does not validate the scientific workflow.

## 3. Phases, code, inputs and genuine outputs

| Phase | Executor | Input | Output and boundary |
| --- | --- | --- | --- |
| Somatic variant evidence | `full.alignment`, `full.somatic` | matched tumor and normal FASTQ/BAM, GRCh38 and VEP | GATK somatic VCF, filtered calls, transcript-specific KRAS variant panel; purity, CNV/LOH and phase require independent evidence |
| RNA-editability filtering | `scripts/adar_pipeline.py` | version-verified KRAS CDS and coding HGVS variants | direct A-to-I codon reversibility; not an editing-rate estimate |
| Guide generation | `scripts/adar_pipeline.py` | correctable KRAS variants | canonical guide sequences, A-C mismatch and bystander candidates |
| Experimental chemistry provenance | `scripts/ingest_publisher_workbooks.py`, `scripts/transfer_publisher_chemistry.py`, `full.chemistry` | six publication files, especially supplementary and Source Data workbooks | position-specific RESTORE 2.0 pattern transfers to KRAS with source keys; **not validated on KRAS** |
| Local RNA structure | `scripts/advanced.py` | guide/target canonical RNA | ViennaRNA outputs, not modification-aware free energies |
| Sequence complementarity | `full.offtarget` | sense-oriented human transcriptome and guides | Bowtie2 local binding opportunities, not demonstrated off-target editing |
| Protein and RNA complex | `full.structure` | exact ADAR protein isoform FASTA, MSA, canonical RNA | Boltz-2 YAML, predicted ensemble, per-model geometry, **not catalytic mechanism or chemically modified RNA** |
| Chemical topology | `full.parameterize` | **atomically correct** model, verified Amber .off/.lib/.frcmod and independent validation manifest | solvated prmtop/inpcrd with approximate NaCl concentration; **only** if all validations pass |
| Molecular dynamics | `full.dynamics` | verified solvated topology and coordinate files, matched evidence checksums | OpenMM restrained then unrestrained NPT equilibration, production trajectories and conformational summaries; **not measured editing** |
| Chemistry-aware ML | `full.learning` | curated measured editing + positional chemistry + independent target groups | grouped CV and external holdout evaluation, otherwise abstain; **unvalidated for KRAS until actual KRAS evaluation** |
| Experimental read counts | `full.experimental` | matched treated/control RNA BAM, sites, genome FASTA | base fractions with heterozygous allele warnings, not mutant-specific editing efficiency |
| Research prioritization | `full.prioritize` | verified guide metadata and optionally experimental measures | engineering triage; experiment-gated research ranking only |

Important details for `full-e2e`:

* It chooses **one unambiguous mutation source**, never a demo panel silently. With verified somatic output and VEP cache, the actual variant table is fed to design on execute mode.
* Force-field source files, residue stereochemistry and structure atom mapping must be externally validated. The supplied booleans alone are insufficient. Modified models additionally check hashes of the PDB, `.lib/.off`, `.frcmod`, Amber topology and coordinates. `examples/input_templates/parameter_validation.json` deliberately has false validation flags.
* If the output of a parameterization stage is eligible and `md.auto_jobs_from_forcefield: true`, it can pass into MD, **but** modified systems require a separate post-build checksum validation JSON. The workflow does not make this file up.
* Boltz MSA must be supplied or generated with `--use-msa-server`. Every diffusion sample is reported separately. A canonical predicted protein-RNA pose is not a validated active enzyme-substrate state; the catalytic Zn2+ center may be absent.
* ML uses grouped cross-validation and can evaluate truly independent target groups using `./run_all.sh ml validate-external ...`. Groups overlapping the training source are rejected. Even a favorable holdout result is not a KRAS experimental editing guarantee.
* `RUN_STATUS.json` is the authoritative stage ledger. Possible values include `INPUT_REQUIRED_NOT_EXECUTED`, `PLANNED_NOT_RUN`, `SKIPPED_INPUT_REQUIRED`, `COMPLETED`, and `FAILED`. A stage marked completed means its executor finished, not that the biological result is valid.

## 4. Required input fields

The main schema is `config/full_project.yaml`. The following templates are in `examples/input_templates/`: `build_specs.csv`, `md_jobs.csv`, `parameter_validation.json`, `measured_editing.csv`, `measured_chemistry_positions.csv`, `functional_assays.csv`, `rna_editing_sites.csv`.

Critical source requirements:

1. **Variant identity.** Genomic assembly, actual somatic variant, sample/patient provenance, a version-locked coding transcript and matched coding reference sequence. Multimutant haplotypes must be phased before being represented as a semicolon-separated cis allele. No tumor lineage or patient state is inferred from a mutation alone.
2. **ADAR.** Precisely identified ADAR1 p110, p150 or ADAR2 isoform, domain boundaries and localization. ADAR2 structural data alone do not validate ADAR1 activity in PDAC cells.
3. **RESTORE 2.0 chemistry.** Published chemistry transfers preserve source provenance but may include 2'-O-methyl, 2'-fluoro, 2'-deoxy, PO, PS, special bases or conjugates. For this publisher dataset, 104 original records are stored; 34 passed the conservative standard parser and 70 need separate chemical interpretation. Avoid a generic, invented chemistry pattern.
4. **Stereochemistry and simulation.** Rp/Sp for every PS linkage, correct protonation, atom types and partial charges, force-field verification, solvated explicit water/ions, QC across replicas. The standard OpenMM stage cannot predict A-to-I reaction chemistry, cellular delivery or in-vivo effectiveness.
5. **Editing quantification.** At KRAS c.35G>A, edited mutant mRNA can look like WT mRNA in standard A/G counting. Heterozygous WT expression confounds the apparent G fraction. Independent phasing, allele-specific assays and protein signaling measurements are required.

**Key target caveat:** A guide C opposite mutant KRAS G12D A forms A-C mismatch on the mutant but a perfect G-C pair with WT KRAS RNA. This does not ensure mutant-selective occupancy. All candidates require WT transcript counterassessment and bystander/off-target evaluation.

## 5. Suggested run sequence

```bash
./run_all.sh doctor
./run_all.sh smoke
./run_all.sh full-e2e --config config/full_project.yaml  # truthful plan and missing-input report
# Supply real/verified data paths and run the desired external executors
./run_all.sh full-e2e --config config/full_project.yaml --execute-heavy
# Individual stages for detailed debugging
./run_all.sh protein-rna --help
./run_all.sh forcefield --help
./run_all.sh dynamics --help
./run_all.sh ml --help
./run_all.sh audit
```

For all failures, inspect `RUN_STATUS.json` and the stage-specific log. Never use a simulated, synthetic or out-of-domain estimate as a therapy efficacy claim or a synthesis-ready chemical specification.

## 6. Core sources and external software

* Pfeiffer et al. 2025, RESTORE 2.0, Nature Communications, DOI 10.1038/s41467-025-64434-7. The original supplementary PDFs/XLSX supplied by the user are archived under `source_material/`, with a parsed chemical reference under `reference/publisher/`.
* Boltz official documentation: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
* GATK Mutect2 and FilterMutectCalls: https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2 and https://gatk.broadinstitute.org/hc/en-us/articles/5257888603803-FilterMutectCalls
* Ensembl VEP: https://github.com/Ensembl/ensembl-vep
* ViennaRNA: https://github.com/ViennaRNA/ViennaRNA
* Bowtie2: https://github.com/BenLangmead/bowtie2
* OpenMM: https://docs.openmm.org/latest/
* AmberTools: https://ambermd.org/AmberTools.php

For specific fixes, tests and incomplete validation obligations see `docs/FINAL_AUDIT_v1_1.md`. Old project development notes are preserved in `docs/LEGACY_HISTORY_v1_0.md`, but any statement there inconsistent with this README is superseded.

## 7. Prepare checksum scaffolds for custom chemistry

After assembling independently reviewed 3D PDB and Amber residue libraries, generate a **nonvalidated** manifest with exact file digests:

```bash
./run_all.sh evidence parameters --pdb /path/to/modified_duplex.pdb --libs /path/to/custom.lib --frcmods /path/to/custom.frcmod --out /path/to/parameter_validation.json
# After validated tleap output has been created:
./run_all.sh evidence topology --topology /path/to/system.prmtop --coordinates /path/to/system.inpcrd --out /path/to/md_validation.json
```

Both outputs intentionally contain `validated: false`. Appropriate external chemical/structural validation is required before updating the boolean fields to true. Only then do the parameter and MD executors accept modified oligonucleotide inputs. File digests prevent accidental use of stale validation after a topology or parameter file has been edited. A checksum is not a force-field validation.

## v1.3 integrated release audit

This version corrects planned-VEP vs produced-VCF stage transitions, rejects absent annotated VCFs, strengthens competing variant-source checking, and normalizes chemistry linkage labels. It adds dedicated regression tests and a reproducible audit report at `docs/RELEASE_AUDIT_v1_3.md`. The full suite and the synthetic integration path passed. Actual modified ADAR catalysis, experimentally calibrated KRAS editing efficiency, and PDAC therapeutic efficacy remain outside the claims of these software tests.
