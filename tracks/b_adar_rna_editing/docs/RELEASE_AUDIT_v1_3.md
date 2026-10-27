# KRAS ADAR E2E v1.3: Integrated release audit

## Audit scope

Examined the actual v1.2 Python and Bash code, orchestrator control flow, current regression tests, manifest normalization, software-only structural and molecular dynamics execution contracts, ML target-group validation gates, and the packaged published RESTORE 2.0 references.

## Confirmed fixes in v1.3

1. **Somatic VEP stage transition:** The planning-only mode no longer attempts to extract a KRAS panel from a VEP VCF that has not been produced. A VEP output is used only when actually created.
2. **User-supplied VCF validation:** A configured `somatic.annotated_vcf` that does not exist generates an explicit error rather than silently falling back to a demonstration panel.
3. **Competing mutation sources:** A user-provided discovery CSV cannot simultaneously compete with a configured raw VCF queued for VEP annotation.
4. **Chemical-annotation normalization:** Whitespace and case normalization for standard backbone labels and PS recognition are consistent in chemistry manifest validation.
5. **Regression coverage:** New tests specifically check the above data-flow gates.

## Verification on supplied synthetic fixtures

- Python modules: compilation succeeded.
- Bash scripts: `bash -n` succeeded.
- Unit and regression suite: 18 tests passed.
- Synthetic full-e2e orchestration: exit 0, 11 recorded stages (5 completed, 6 appropriately skipped for missing scientific input).
- Synthetic artifacts: 18 guide records, 21 chemistry-transfer constructs.
- No real KRAS ADAR experiments, ML editing efficacy, modified force-field dynamics, or therapeutic efficacy were computed in these tests.

## Important limitations that remain intentionally unclaimed

- No claim of complete software correctness over every external dependency and file format.
- Boltz canonical RNA surrogate complexes do not provide modified-oligonucleotide structures or ADAR base-flipping/catalytic rates.
- RESTORE 2.0 positional pattern transfer is not a validated editing-efficiency predictor for KRAS G12D.
- Stereo-random PS chemistry cannot be represented by a single arbitrary stereoisomer.
- The E2E suite does not independently validate manually supplied force fields, parameterization proofs, biological editing labels, delivery, PDAC responses or toxicity.
- Mature tumor-normal calling and off-target analysis require real reference genomes, indexed sequencing files, aligner models and external-tool version testing on the user's compute setup.
- The current RNA-seq read-count module must not be interpreted as allele-specific correction without phasing or orthogonal allele-resolving measurement.
- ML target-group CV on historical non-KRAS editing experiments does not establish external validity for RESTORE-modified KRAS guides.

## Release rule

v1.3 is the audited code release. Use a fresh `config.out` for every execution and verify every missing-input/conditional stage in `RUN_STATUS.json`. `FINISHED_WITH_CONDITIONAL_MODULES` is not equivalent to completion of a scientific efficacy study.
