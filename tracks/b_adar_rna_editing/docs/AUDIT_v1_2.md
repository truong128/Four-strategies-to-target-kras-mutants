# Final targeted audit: v1.2

## Verified issues corrected
1. **Stale artifact reuse**. Previous workflow allowed writing into a previously populated `out` directory. This could retain an old ML model, annotated call set, or RNA off-target alignment file and misattribute it to a new run. Now existing nonempty output is rejected. Use a unique output directory per run.
2. **ML model evidence gate**. An existing model file by itself is insufficient. Downstream ML steps require the current run's report to assert the expected qualifying status.
3. **Terminal linkage validity**. Position-level chemistry encoding now forbids a PS or PO internucleotide linkage after the last (3-prime) guide residue.
4. **Legacy synthetic fixture** fixed to represent a physical terminal oligonucleotide.
5. **Off-target provenance**. Ranking only uses Bowtie2 results if the transcriptome screen was explicitly requested and actually executed in the current run.

## Remaining limits
- Published RESTORE modifications were established in other target contexts, not validated in KRAS G12D pancreatic cancer.
- Boltz outputs on canonical RNA do not establish a chemically modified ADAR1 catalytic state.
- Missing bespoke force-field parameters must be validated outside of the pipeline and tied to hash evidence.
- Chemistry-aware training requires accurately matched editing assay labels and independent target validation.
- Somatic calling and structure/MD external executables were not run on real patient or wet-lab data in this audit.
- Therefore this is a **source-code correctness improvement**, not a validated therapeutic predictor.
