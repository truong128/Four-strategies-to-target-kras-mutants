# KRAS-ADAR computational E2E candidate-discovery package (v0.1 research prototype)

## Scope and honesty

This is a **working one-command mutation-to-guide candidate generator**, not a validated end-to-end ADAR therapeutic design and not the previously discussed fully parameterized ADAR-RNA MD suite. Outputs explicitly mark unperformed analyses as BLOCKED/NOT_RUN. Do not use guide sequences for synthesis or animal/human administration without confirming transcript/variant identities, independent alignment, chemical patterns, and experimental performance.

## Workflow implemented
1. Download and cache Ensembl canonical KRAS CDS using Ensembl REST (or provide a locally verified coding-sequence FASTA).
2. Validate CDS triplet/frame and codon-12 Gly reference; validate HGVS coding SNVs against FASTA.
3. Compute mutant CDS/protein and check **direct A-to-I functional reversal** of each altered amino acid.
4. Generate antisense guide sequences with target A-C mismatch for different guide lengths and offsets.
5. Tabulate all target-region adenosines as possible bystander sites and WT versus mutant sequence controls.
6. Emit a conservative chemical-modification **placeholder manifest** (unmodified control; exploratory terminal 2′-OMe/2′-F/PS). This is NOT a verified RESTORE 2.0 modification positional recipe.
7. Optionally evaluate unmodified **single-strand** guide secondary structures using installed `RNAfold` (NOT guide-target duplex energy, NOT modification-aware).
8. Generate a BLOCKED molecular simulation input manifest. No unparameterized chemistry is simulated.
9. Save auditable tabular results and completion status.

## Quick start: macOS/Linux, Python 3.9+ (no mandatory pip dependencies)

```bash
unzip KRAS_ADAR_E2E.zip
cd KRAS_ADAR_E2E
./run_all.sh smoke
./run_all.sh discover
```

`discover` needs internet to retrieve Ensembl CDS on first execution. Optional ViennaRNA: `conda install -c bioconda viennarna`; then run `./run_all.sh discover --fold`.

Offline mode with a **verified complete KRAS canonical CDS**:

```bash
./run_all.sh discover --cds-fasta /path/to/verified_KRAS_CDS.fasta
```

Customize `config/variants.csv` using HGVS `c.35G>A` format. Semicolon-separated mutations in one row mean **cis on one allele**, which must be confirmed with phased patient data. Indels, fusions, trans mutations, isoform-specific coordinates and RNA splicing mutations are not supported yet.

## Files

- `results/run/01_variants.csv`  -  only direct coding SNV correction classification
- `02_guides.csv`  -  antisense guide candidates
- `03_bystander_sites.csv`  -  enumerated A positions, NOT predicted editing efficiencies
- `04_chemistry_manifest.csv`  -  reference/exploratory chemistry, not validated synthesis specification
- `05_wt_comparison.csv`  -  sequence comparison, NOT proven allele specificity
- `06_secondary_structure.csv`  -  optional RNAfold output (single guide strand)
- `07_offtarget_status.json`  -  explicit off-target status, not a negative screen
- `08_md_input_manifest.csv`  -  blocks molecular dynamics until molecular inputs/parameters exist
- `run_manifest.json`, `SUMMARY.md`  -  audit trail/limitations

## Mandatory scientific and implementation gates (NOT IMPLEMENTED)

**Genomics:** Full PDAC variant discovery requires tumor/normal WGS/WES/RNA-seq analysis; VCF normalization; somatic evidence and clonality; copy number; purity, phase/cis/trans; transcript version lock, WT expression, RNA abundance. This tool starts at a **variant panel**, not raw FASTQ and does not discover mutations in patient data. Clinical or research sequence provenance must be independently validated.

**Guide design:** Pseudosequence complementarity alone does not predict editing. Target and guide accessibility, ADAR1 p110/p150 vs ADAR2, isoform localization/expression, loop/bulge engineering, secondary interactions, RNA secondary structures, mutant vs WT discrimination, and effective concentrations remain unmeasured.

**Modification:** Exact published RESTORE 2.0 patterns must be extracted from source supplementary tables or explicit chemistry files with positional, linkage and stereo information. `2'-OMe + 2'-F + PS` is a family, not one unambiguously defined oligo. Never label the generated placeholder as RESTORE 2.0 itself.

**Off-target:** To claim transcriptome-scale coverage one must run a documented alignment/search against relevant human transcript isoforms (and assess guide mismatches, bulges, RNA editing context, native editing baseline, and gene consequence). This package intentionally makes no such claim.

**Structure/MD:** The prior Amber duplex workflow (Amber RNA OL3 / custom parameters, ~150mM NaCl, 310.15K, explicit solvent) can be reused *after* 3D modified oligos, accurate charges, PS Rp/Sp ensembles, residues/force field parameter validation, and protein-RNA complex preparation. RNA duplex MD does not measure ADAR editing, and 100ns classical MD does not demonstrate A-base flipping nor chemistry. ADAR-RNA modeling demands catalytic-state structures and stereochemical reproducibility.

**Biology:** KRAS reversion at RNA does not equal tumor suppression. Need cellular editing measured from RNA sequencing and matched DNA controls, guide localization/uptake, WT allele monitoring, bystander analysis, ADAR dependency, protein-level KRAS-GTP/pERK measurements, multiple PDAC cellular models, tumor microenvironment delivery and safety/toxicity studies. Assess RNA turnover, mutant allele dosage, resistance alleles and heterogeneity.

**Model calibration:** No measured editing efficiencies supplied means no ML training, probability, accuracy or performance claims. Future data splits should be by target locus/sequence context and chemistry rather than random guide duplicates.

## Sources and reproducibility

- RESTORE 2.0: https://www.nature.com/articles/s41467-025-64434-7
- G-U wobble/circular CLUSTER: https://www.nature.com/articles/s41587-024-02313-0
- Endogenous guide review: https://academic.oup.com/nar/article/53/13/gkaf672/8203451
- Ensembl REST sequence API: https://rest.ensembl.org/documentation/info/sequence_id

These references are for design context; this package does not implement the published models exactly.

## v0.2 extension: executable steps 7-12

Run `./run_all.sh e2e` (online Ensembl reference), or the two stages separately:

```bash
./run_all.sh discover --cds-fasta /path/to/VERIFIED_KRAS_CDS.fa
./run_all.sh advanced --input results/run --out results/run \
  --transcriptome /path/to/human_transcripts.fa --max-mismatches 2
```

`advanced` produces:

- `07_chemistry_constructs.csv`, `07_chemistry_positions.csv`: concrete per-position **exploratory** sugar/backbone tables. They are **NOT exact, experimentally validated RESTORE 2.0 recipes**, and contain no unstated Rp/Sp stereochemistry.
- `08_cofold.csv`: if `RNAcofold` exists, duplex MFE computed for **unmodified** RNA only; otherwise unavailable status.
- `09_transcriptome_hits.csv` and `09_offtarget_limits.json`: actual full-length reverse-complement local scanning for up to `--max-mismatches` substitutions on both orientations, with pigeonhole exact-seed acceleration. This is only a sequence-complementarity risk screen, not transcriptome-wide ADAR editing or gapped/seed-only interaction prediction. FASTA choice, alternate isoforms and expression matter.
- `10_md_readiness.csv`: gates simulation on independently verified coordinates, topology, forcefield/stereochemistry and template. It does not create structural models or fit force fields.
- `11_editing_ml.csv`, `11_model_status.json`: requires measured experimental editing fractions + independent target groups; group CV and baseline comparison; no experimental data = ABSTAIN. Even a model that beats baseline is NOT validated for KRAS without external target-specific testing.
- `12_candidate_evidence.csv`: evidence ledger with measured unknowns, no invented therapeutic performance or numerical overall rank.

For optional experimental data fitting (only when real labels are available):

```bash
python3 -m pip install -r requirements-ml.txt
./run_all.sh advanced --measured-editing /path/to/measurements.csv
```

Labeled CSV columns: `guide_rna_5to3,target_rna_5to3,editing_fraction,study_or_target_group` (fraction in [0,1]). Requires >=30 observations and >=3 genuinely distinct groups. The simplistic sequence features are an **exploratory baseline**, not a publication-grade predictor.

For OPTIONAL Amber run, once you have validated, pre-equilibrated solvated periodic systems, prepare a CSV whose fields include:
`construct_id,structure_pdb,topology_prmtop,coordinates_inpcrd,chemistry_verified,stereochemistry_verified,parameterization_validated,template_validated,equilibrated_explicit_solvent_ions`.
Set each flag to `YES` only if independently confirmed. Execute:

```bash
python3 scripts/md_runner.py --manifest /path/to/validated_md.csv
# inspect generated command.json and production.in; then:
python3 scripts/md_runner.py --manifest /path/to/validated_md.csv --execute --engine pmemd.cuda
```

This executes short classical MD of **supplied, validated** structures; it does NOT build ADAR protein-RNA complexes, simulate catalysis or demonstrate base flipping. For chemically modified duplexes, validated Amber forcefields and stereochemistry are prerequisite inputs. Mac users without Amber CUDA can select `sander`/`pmemd`.

### Important unfinished scientific work
Patient FASTQ→VCF variant calling/phasing; published exact chemical positional recipes; automated parameterization of arbitrary modified RNA; validated ADAR-RNA complex modeling; editing outcome external validation; PDAC-specific guide uptake and tumor-response assays are not delivered. They cannot be honestly claimed completed just because the pipeline generates files.

## v0.3 RESTORE 2.0 source-verified integration (2026-10-08)
**Critical reproducibility correction:** the RESTORE 2.0 publication describes a chemical editing platform, not a demonstrated turnkey editing-prediction software. Previous v0.2 terminal chemistry constructs are *illustrative*, NOT RESTORE 2.0.

Publication: Pfeiffer et al., Nat Commun 16:8849 (2025), DOI 10.1038/s41467-025-64434-7, PMCID PMC12494997.
Supplementary Data 1 and Source Data must be downloaded from the official publication and audited before an oligonucleotide is marked as an exact published reference. As of this build, publication workbook bytes could not be retrieved in the build environment (network access blocked). **No author-specific modified positions or numeric editing fractions have been fabricated, inferred or preloaded.**

### Integrating source data
1. Download **Supplementary Data 1** and **Source Data** from the paper. You can use `python3 scripts/restore2_download.py --url 'ACTUAL_OFFICIAL_SUPPLEMENT_URL' --output reference/restore2_supp.xlsx --expect xlsx` if direct browser download is inconvenient. Never use a guessed link.
2. Review workbook tabs; normalize published oligonucleotide position records into `config/restore2_reference_TEMPLATE.csv` format using literal positions/chemistry from the publication. Evidence `PUBLISHED_EXACT` is allowed only if each row can be traced to a specific original sheet/row; attach workbook using `--publisher-supplement` for SHA256 provenance. A sequence in a spreadsheet without a positional chemistry map is insufficient. Link every inter-nucleotide phosphate bond (3' outgoing linkage for base i); terminal base must be `NONE`.
3. Add exactly one `is_orphan_C=YES` at the guide nucleotide opposite the intended edited A, *only if confirmed for that published construct*. Important: source editing-site orientation must match the KRAS guide 5'→3' orientation.
4. Run:
```
python3 scripts/restore2_reference.py \
  --reference-csv config/restore2_reference_filled.csv \
  --publisher-supplement reference/restore2_supp.xlsx \
  --publisher-source-data reference/restore2_source.xlsx \
  --guides-csv results/run/02_guides.csv --out results/run/restore2
```
5. Successful imported patterns generate `restore2_transferred_chemistry.csv`, not experimental predictions. Designs with different length or editing-site position are rejected rather than extrapolated.

### Scientific hard stops
- Published source editing efficiencies should **not** be silently attached to transferred KRAS guides. Match reference IDs and assays to Source Data, then separately perform grouped benchmark analyses.
- `RNAcofold` default parameters are unmodified-RNA approximations. Actual 2'-OMe, 2'-F, DNA and PS/PO chemistries require appropriate thermodynamics/atomistic parameters.
- A full transcriptomic match screen is NOT equivalent to off-target RNA editing measurement; include known RNA isoforms and endogenous editing controls.
- PS_RANDOM is a mixture of stereoisomers and cannot be represented by one uniquely specified atomistic topology. For MD sample explicitly specified Rp/Sp microstates or report the limitation.
- MD of duplex and ADAR-RNA structural models do not predict editing yields. Neither molecular docking nor geometry establishes catalytic transition-state accessibility.
- Do not infer PDAC efficacy from in-vitro editing; allele selectivity, RNA delivery, molecular rescue and in-vivo data remain experimentally unverified.

### Acceptance criteria
A truly publication-reproduced RESTORE 2.0 construct needs all of: original source XLSX + SHA256, precise source row, target, guide sequence, each sugar substitution, each bond, editing center and assay outcome linked without manual ambiguity. `restore2_reference.py` rejects incomplete positional maps and mismatching sequence/base labels. ML development is blocked without proper measured labels and independent validation. This is deliberately a **source-auditable research pipeline**, not a proven therapeutic end-to-end design.

### One-command orchestration with verified chemistry
Once a complete normalized reference CSV is present:
```
RESTORE2_REF_CSV=config/restore2_reference_filled.csv ./run_all.sh e2e --cds-fasta /path/to/verified_KRAS_CDS.fa
```
Without `RESTORE2_REF_CSV`, E2E runs only its original exploratory-chemistry path and DOES NOT claim exact RESTORE 2.0 integration. A successful mutation-to-guide run is not evidence of clinical target validation.

## v0.4  -  Publisher-workbook ingestion (uploaded 8 Oct 2026)

This package **contains publisher-supplied workbook evidence** and also retains all six original files under `source_material/`. Run:

```bash
./run_all.sh publisher --supplement source_material/41467_2025_64434_MOESM3_ESM.xlsx --source-data source_material/41467_2025_64434_MOESM6_ESM.xlsx --out reference/publisher
./run_all.sh e2e
```

The first command uses `artifact_tool` (`pip install artifact-tool` if unavailable). It outputs `publisher_constructs.csv`, `publisher_residues.csv`, `publisher_parser_issues.csv`, `source_data_assay_index.csv`, and `publisher_audit.json`. The second command transfers *only* published chemistry templates matching both guide length and orphan-C position, producing `results/run/restore2_published/`.

**Audit results of supplied workbooks:** 104 sequence records, 34 conservative standard-alphabet records, 70 requiring special-base/conjugate/linker interpretation (conjugate, linker, LNA/terminal annotation, etc.); 7 source-data worksheets and 163 assay index rows, of which 96 have shared/missing version-to-target disambiguation. Editing outcomes are NOT joined to variants by ID alone, and no unvalidated ML labels are generated. 2′-F/2′-OMe/PS transfer is not evidence of efficacy in KRAS. PDB generation, chemical custom residues and MD require validated molecular structures. Multiple PS stereocenters are an ensemble, not one stereoisomer.

**Critical biological assumptions:** Individual somatic SNVs from c.35G>A are not patient VCF calls; RNA allele expression and cis/trans configurations require patient material. A guide may bind wild-type KRAS despite an A-C mismatch in mutant, and WT can acquire unwanted edits at bystander adenosines. Off-target sequence searches are not an ADAR-editing specificity assay. Efficacy in liver SERPINA1 is not evidence of efficacy/delivery in pancreatic cancer. The oncology-specific safety of decreasing oncogenic KRAS activity should be experimentally verified. Do not interpret RNAfold/RNAcofold energy as modified-guide free energy without calibrated parameters.

**Peer-review caution:** Source reviewer file explicitly questioned therapeutic-index/safety of 2′-F/PS, PK/delivery and effects on RNA/protein binding. The manuscript supplementary RNAfold figure reports no clear correspondence of minimum-energy hairpins with editing success. Treat these as mechanistic uncertainties, not validated ranking predictors.

## v0.5: ADAR pretrained-model compatibility bridge (external baseline)

A public **AdarEdit** repository (https://github.com/Scientific-Computing-Lab/AdarEdit) includes pre-trained graph-model checkpoints. Its published model analyzes natural endogenous editing in Alu-derived substrates and **does not validate the effect of 2'-OMe, 2'-F or PS/PO chemical modifications on RESTORE 2.0 guides**. This project must not report AdarEdit results as KRAS G12D measured or predicted editing percentage.

- `./run_all.sh adaredit-fetch`: fetch the upstream repository with Git when network access is available (not included or vendored here), store checkout revision. Install it separately using the upstream README and test its own published checkpoints first. This environment cannot download the repository, so no model checkpoints are bundled and **no checkpoint inference has been run**.
- `./run_all.sh adaredit-import --guides results/run/02_guides.csv --predictions external_native_site_scores.csv --out results/run --model-checkpoint /absolute/path/to/checkpoint --model-context Combined --model-version UPSTREAM_GIT_SHA --prediction-type native_site_editability_probability`: import *genuine separately computed* model scores with integrity checks, provenance and **OUT_OF_DOMAIN** labels. Input requires `guide_id,model_output` columns and full guide coverage. This is an evidence bridge, not an inference engine.

The previous Step 11 grouped-ML experimental training option remains exploratory and requires correctly paired measured outcomes. Original RESTORE 2.0 workbook measurements should be audited and grouped by **target, construct, assay, cell type and experimental conditions** before model fitting. Do not merge by ON version number alone. Prospective validation in KRAS-mutant PDAC is mandatory before calibrated editing efficiency claims.

## v0.6 Scientific scope and reproducibility audit

Run `./run_all.sh audit` after `e2e`. Outputs `results/audit/{AUDIT_REPORT.md,audit.json,stage_matrix.csv,tool_inventory.csv}`. `./run_all.sh validate-inputs --cds <verified_CDS.fa>` rejects incompatible KRAS CDS and HGVS coding variant positions. This is **not** a clinical somatic variant caller: use matched tumor/normal VCF, GRCh38 VEP and transcript-version validation upstream. The included `config/variants.csv` is a synthetic mutation panel, not patient observations.

### Tools reviewed and recommended, with limitations

* [Ensembl VEP](https://github.com/Ensembl/ensembl-vep) - input VCF consequence annotation with versioned offline cache. Must keep GRCh38, transcript version, and allele phase.
* [ViennaRNA](https://github.com/ViennaRNA/ViennaRNA) - RNAfold/RNAcofold reference; standard RNA thermodynamics do **not** represent mixed 2'-F/2'-OMe/PS chemistry.
* [AdarEdit](https://github.com/Scientific-Computing-Lab/AdarEdit) - published graph models/weights for natural Alu-based A-to-I editability, **out of distribution** for RESTORE modified guides. The repository is not automatically invoked in v0.6.
* [EditPredict](https://github.com/wjd198605/EditPredict) - independent pretrained natural ADAR baseline, not a calibrated chemically modified guide percentage.
* [JACUSA2](https://github.com/dieterich-lab/JACUSA2) and [REDItools2](https://github.com/BioinfoUNIBA/REDItools2) - processing experimental RNA-seq BAM data to quantify edit sites, requiring DNA/SNP/artifact and matched-control filtering.
* [AmberTools](https://ambermd.org/AmberTools.php) - requires verified modified residues, topology, PS stereochemistry and ADAR protein-RNA structure. No automatic modeling of reactive deamination.

**Critical corrected issue:** sense-oriented transcriptome FASTA is searched only in the biological transcript orientation. The earlier v0.5 scan of reverse-complement strands is not appropriate as evidence of expressed RNA binding. Genomic screening is a separate task; full-length ≤2-mismatch screening still does not cover partial duplexes, bulges, or transcript abundance.

### Gate statuses

`CODED` means script exists, not experimentally validated. `CONDITIONAL` requires external data/executables. `GATED` requires verified input model and measured data. `NOT_AVAILABLE` means it is not implemented; the audit never substitutes placeholder values.

### Minimum acceptance tests before therapeutic conclusions

1. Matched tumor/normal confirmation of actual mutant KRAS allele, heterozygosity, copy number, allele phase and expression.
2. Allele-specific RNA editing measurement with appropriate controls and full transcriptome off-target profiling; ADAR1 isoform/context profiling.
3. ADAR-RNA structure and modified forcefield accuracy checks against independent published or measured data.
4. Independent experimental KRAS G12D guide editing rates, KRAS signaling suppression, PDAC cellular efficacy, delivery and toxicity.

Neither AdarEdit output nor RNAcofold MFE nor all-atom RMSD should be interpreted as editing percentage or pancreatic tumor therapeutic benefit.

### Executable external-tool commands

```bash
./run_all.sh external setup-conda
./run_all.sh external vep input/variants.GRCh38.vcf /path/to/vep_cache results/annotated.vcf
./run_all.sh external samtools-check input/treated.bam
./run_all.sh external adaredit-clone external/AdarEdit
```
These are actual wrappers, not fake inference. `vep` requires a predownloaded current, assembly-consistent GRCh38 cache. `adaredit-clone` clones source and prints its commit but does not run a non-generalizable model on RESTORE guides. Do not mix liver-trained or Alu-trained model output with validated chemical editing percentages. JACUSA2 and REDItools2 are documented analysis choices, not invoked without properly designed treated/control BAM cohorts.

## v0.7 independent structural/ML readiness audit (2026-10-08)

Run `./run_all.sh structural-audit` after `./run_all.sh e2e` has produced `results/run/02_guides.csv`. To test without network on the bundled synthetic example, run:

```bash
./run_all.sh structural-audit --guides results/synthetic/02_guides.csv --transferred results/synthetic/restore2_published/transferred_per_residue.csv
```

This builds `results/structural_ml_preflight/structural_tasks.csv`, `SCIENTIFIC_READINESS.json`, and `READINESS_REPORT.md`. It verifies residue-level chemistry map integrity, flags PS stereo-ensembles, checks claimed MD evidence against actual files, and rejects editing-specific ML status where curated chemistry-aware, group-independent experimental labels are missing. It does **not** build validated modified residue force fields, perform ADAR docking or catalytic MD, or predict measured KRAS editing efficiency. Neither RESTORE Supplementary Data nor a pretrained Alu-specific model proves such capabilities. The new preflight module is a real executable evidence gate, not a replacement for those scientific tasks.

---

# v1.0 Full-source execution package (2026-10-08)

**User contract:** All source-code stages are implemented and runnable with their declared dependencies and correct supplied inputs. This is not a claim that a validated KRAS G12D ADAR editing model, complete modified-chemistry 3D coordinates, or a biologically effective therapy currently exists. Heavy third-party executables are intentionally **not** run or pretrained on this computer.

## One-command orchestration

```bash
# Micromamba or conda installation on YOUR machine (Linux with CUDA recommended for Boltz):
micromamba create -f environment-full.yml
micromamba activate kras-adar-full

# Check optional dependencies and what is missing:
./run_all.sh doctor

# Copy and edit all user input paths in config/full_project.yaml
./run_all.sh full-e2e --config config/full_project.yaml

# When FASTQ/BAM, verified ADAR isoform FASTA, transcriptome, modified topology,
# validated chemistry force fields and measured editing data are actually present:
./run_all.sh full-e2e --config config/full_project.yaml --execute-heavy
```

The non-heavy run **creates Boltz/Amber/Bowtie plans but never claims these have executed**. `RUN_STATUS.json` labels each module `COMPLETED`, `PLANNED_NOT_RUN`, `SKIPPED_INPUT_REQUIRED` or `FAILED`. Heavy execution requires the actual programs installed and the matching input files.

### Download reproducible references (optional)

```bash
./run_all.sh references gencode --release 48 --out inputs/references --execute
./run_all.sh references protein --accession YOUR_VERIFIED_UNIPROT_ADAR1_ISOFORM_ID --out inputs/references --execute
```

Always check protein isoform and GENCODE/GRCh38 version consistency; reference download alone does not establish ADAR expression in PDAC.

### Standalone executors (source is in `full/`)

| Module | Command | Input required | Actual output/claim |
|---|---|---|---|
| Raw sequencing | `./run_all.sh align ... --execute` | tumor/normal paired FASTQ, GRCh38 genome, sample IDs | bwa-mem2+samtools+GATK duplicate-marked BAM (optional BQSR) |
| Somatic mutation | `./run_all.sh somatic mutect2 ... --execute` | indexed matched BAMs, references, artifact/QC inputs | Mutect2 + FilterMutectCalls filtered VCF |
| Effect annotation | `./run_all.sh somatic vep ... --execute` | VEP cache, GRCh38 VCF | annotated VCF |
| KRAS panel | `./run_all.sh somatic panel ...` | VEP `CSQ` VCF | HGVS c. SNV table; no cis/trans grouping |
| Candidate generation | `./run_all.sh discover ...` | verified KRAS CDS, allele variants | guide sequences and bystander table |
| RESTORE 2.0 | `./run_all.sh publisher` then `transfer` | source XLSX included | position-level 2' chemistry pattern with provenance |
| Chemistry validation | `python3 -m full.chemistry ...` | published transfer CSV | construct-level/position chemistry manifest |
| RNAfold/cofold | `./run_all.sh advanced ...` | ViennaRNA executable | UNMODIFIED-RNA approximation |
| RNA off-target | `./run_all.sh rna-offtarget index/screen ... --execute` | strand-normalized human transcriptome | Bowtie2 local hybridization candidate matches |
| Protein/RNA modeling | `./run_all.sh protein-rna prepare/predict/geometry ...` | verified ADAR isoform FASTA, canonical RNA guides, GPU+Boltz | canonical protein-RNA structure models, geometry, not catalysis |
| Custom chemistry | `./run_all.sh forcefield ... --execute` | custom `.lib/.off`, `.frcmod`, validated atomistic PDB, PS Rp/Sp proof | Amber topology + coordinates, solvated water/ions at specified molarity |
| MD | `./run_all.sh dynamics simulate/analyze ... --execute` | verified solvated Amber prmtop/inpcrd, OpenMM | per-replicate DCD and conformational measurements |
| Chemistry-aware ML | `./run_all.sh ml train/infer ...` | curated guide-level measured editing rates and per-position chemistry | group-held-out predictive metrics, OOD predictions |
| RNA-seq validation | `./run_all.sh rna-assay ...` | matched RNA BAM controls and known edited sites | nucleotide proportions, NOT allele-specific editing fraction |
| Evidence ranking | `./run_all.sh ranking ...` | verified measurements for ranking; basic guide metadata | engineering triage; measured in-vitro rank only when evidence available |

### Required file formats

`forcefield.build_specs_csv`: `construct_id,pdb,modified,libs,frcmods,validation_json,salt_pairs` where `libs` and `frcmods` are semicolon-separated absolute or package-relative paths; `salt_pairs` can be left blank to estimate pairs from first-stage water count for 0.15 M bulk NaCl.

`md.jobs_csv`: `construct_id,topology_prmtop,coordinates_inpcrd,modified,validation_json`. Modified topologies require `chemistry_validation.json` containing `validated`, `stereochemistry_resolved`, `charge_and_torsion_checked`, `atom_mapping_checked` all `true`. This is a user-supplied independent scientific validation record, not proof synthesized by script.

`machine_learning.curated_labels_csv`: `construct_id,guide_rna_5to3,target_rna_5to3,guide_index_0based,target_group,assay_context,editing_fraction`. `curated_chemistry_positions_csv`: `construct_id,position_1based,base,sugar,link_to_next_3prime`. Training demands enough distinct target groups and strictly grouped cross-validation. Original RESTORE source index is **not automatically treated as unambiguous measured labels** because sample ID alone does not distinguish assay/target; curate from paper and experiments first.

`ranking.experimental_csv`: `guide_id,mutant_edit_fraction,wt_bystander_fraction,max_offtarget_edit_fraction,KRAS_GTP_reduction_fraction,cell_viability_effect_fraction,assay_batch,biological_replicates`.

`sites.csv` for RNA data: `site_id,contig,position_1based,transcript_strand,genomic_reference_base`. **For KRAS G12D heterozygotes, a G read can come from WT transcript OR A-to-I edited mutant RNA; an A/G read ratio is not directly an editing rate.** Use phasing/matched DNA, ADAR perturbation, proteomics and functional assays.

### Essential scientific interpretation and computational boundaries

1. **Allele discrimination cannot be assumed.** On KRAS G12D mutant mRNA, a C opposite mutant A gives an A-C mismatch, but the same C forms a canonical G-C pair with WT KRAS at the site. WT hybridization can be strong. WT bystanders and ADAR recruitment must be measured independently.
2. **Boltz-2 canonical RNA is a structural reference only**. It cannot reliably simulate the complete stereo-random PS backbone, altered ribose chemistries or ADAR Zn chemistry as a black box. Your own atom-mapped modified structure and validated force-field libraries are needed.
3. **Classical MD does not model enzymatic conversion A→I** or editing frequency. Actual catalytic free-energy barriers require validated enhanced-sampling/QM/MM protocols and catalytic geometry, not trajectory RMSD.
4. **RESTORE labels must be curated** across target, assay, construct and cell conditions. An ML model can be trained with the source code provided, but independent PDAC/KRAS validation cannot be created by programming alone.
5. **Biological delivery remains separate**: RESTORE in vivo mouse liver / hepatocyte results do not validate pancreatic tumor delivery; off-target RNA-seq and ADAR1 expression need tumor-specific data.
6. **G12D direct correction does not eliminate other mutations** or prove KRAS-GTP suppression, tumor regression, toxicity or therapeutic index.

Versioned public sources: [RESTORE 2.0](https://doi.org/10.1038/s41467-025-64434-7), [Boltz2](https://github.com/jwohlwend/boltz), [OpenMM](https://docs.openmm.org/latest/), [Ensembl VEP](https://github.com/Ensembl/ensembl-vep), [ViennaRNA](https://www.tbi.univie.ac.at/RNA/), [Bowtie2](https://github.com/BenLangmead/bowtie2), [GATK Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2). These tools have different scopes and must not be treated as one validated clinical model.
