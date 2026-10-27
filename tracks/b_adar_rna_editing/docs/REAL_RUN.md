# Real-run setup (KRAS variants in PDAC)

## What the setup does

- `scripts/setup_real_run.py` builds the KRAS CDS (strand-aware, with stop codon, validated: 188 aa, codon 12 = GGT), a transcript-to-gene map from the GTF, and the ADAR1 FASTAs (UniProt P55265; p150, p110 = isoform 5 residues 296-1226, DRBM3+deaminase 716-1226), each with provenance JSON.
- Off-target screening is gene-aware: hits carry `gene_name` and `is_on_target_gene`; KRAS isoform hits are excluded from the off-target count, and `offtarget_gene_summary.csv` lists the genes hit. A screened guide with zero hits is reported as 0 (it was previously `NOT_SCREENED`).
- Boltz predict now stages only the `*.yaml` files into `boltz_models/boltz_yaml`; previously `structure_index.csv` and the disclaimer JSON in the input folder made Boltz abort ("Unable to parse filetype .csv").
- `bowtie2-build` and the screen use `--threads`.
- Boltz stage: `structure.variants` and `structure.max_guides` filters; `--no_kernels` is on by default (required on Blackwell GPUs).
- `preflight` checks, before any heavy step: empty output directory, free disk, variant table, CDS sanity, HGVS reference bases against the CDS, transcriptome, gene map, bowtie2, ADAR FASTA, boltz and weights, GPU, RNAfold/RNAcofold, and tleap/OpenMM or GATK/VEP when those stages are configured.
- Inputs: `config/kras_pdac_variants.csv` and `config/real_run_pdac_kras.yaml`.

## Transcript choice

Use ENST00000311936 (KRAS-4B, MANE, 188 aa). ENST00000256078 is KRAS-4A (189 aa). Codons 12, 13 and 61 are identical in both. No FASTQ, VEP, or GRCh38 alignment is needed for guide design from a variant table.

## Which variants ADAR can correct

ADAR deaminates A to I (read as G). Only G>A mutants are directly correctable:

| Variant | Change | Correctable |
|---|---|---|
| G12D | c.35G>A | yes |
| G12S | c.34G>A | yes |
| G13D | c.38G>A | yes |
| G12V, G12C | G>T | no |
| G12R | G>C | no |
| Q61H | A>C / A>T | no |

The non-correctable variants are included so the pipeline reports them as non-correctable rather than silently skipping them. The guide places an orphan C opposite the target A. The wild-type G-C pairing means mutant-vs-wild-type selectivity is not established by this pipeline.

## Commands

```bash
cd KRAS_ADAR_E2E_v14_adversarial_audit
./run_all.sh smoke                                   # 31 tests
./run_all.sh setup-real                              # CDS, gene map, ADAR1 FASTAs (already generated)
./run_all.sh preflight --config config/real_run_pdac_kras.yaml
./run_all.sh full-e2e  --config config/real_run_pdac_kras.yaml --execute-heavy
```

The output directory (`results/real_run_01`) must not already exist.

## Still not runnable without your inputs

- MD and force-field stages need validated parameters for the modified residues and PS stereochemistry.
- The ML stage needs curated measured editing data and abstains otherwise.
- Boltz-2 models use canonical RNA with no Zn2+ and no chemical modifications, so they give geometry only.
- With `use_msa_server: true`, the public ADAR1 sequence is sent to the ColabFold server. Set it to false and supply `structure.msa` to avoid that.
- `FILE_INTEGRITY.json` and `RELEASE_SHA256.json` predate these edits and are stale.
