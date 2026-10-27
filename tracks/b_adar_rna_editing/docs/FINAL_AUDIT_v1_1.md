# KRAS ADAR E2E source-code audit, v1.1

## Confirmed code corrections

1. Safety and input identity: Real-project config no longer points to demonstration HGVS variants. Synthetic reference data require an explicit synthetic-only configuration. A missing patient mutation source generates an `INPUT_REQUIRED_NOT_EXECUTED` report rather than a false run.
2. Tumor-normal workflow: Updated GATK Mutect2 arguments, added optional F1R2 orientation-bias modeling and user-supplied contamination table support. Extractor now retains only `PASS` VCF calls and checks VEP ALT-allele agreement.
3. Variant transfer: Connected filtered somatic calls through VEP to the mutation-to-guide executor on actual execution when no competing manually provided panel exists.
4. RESTORE chemistry: Published per-position transfer requires correct bases, contiguous positions and valid 3' terminal linkage. It cannot claim successful editing from chemical pattern matching.
5. Boltz-2: A protein MSA is required unless the server is enabled or single-sequence usage was explicitly selected. Each diffusion sample has its own geometry entry; no silent choice of one predicted pose. No affinity or catalytic score is fabricated.
6. Modified parameterization: The actual structure PDB, residue libraries and parameter files must match SHA-256 digests in independent validation records. Amber RNA OL3, DNA OL15 and protein ff14SB reference force fields are loaded, but custom residues still need independent QC.
7. OpenMM: Added a staged restrained then unrestrained NPT equilibration before production, with independent replica seeds. Modified MD is blocked if validated topology or coordinate hashes do not match evidence.
8. E2E MD handoff: When a valid Amber build succeeds and automatic MD handoff is enabled, job construction is automated. It never invents the separate post-build modified-topology validation. A conformational summary is generated from completed trajectories.
9. Editing ML: Added consistency checks for nucleotide-specific chemistry, duplicate positions and identical target sequences assigned inconsistent holdout groups. External target-group evaluation is available; without truly independent data, claims remain abstained.
10. RNA-seq and ranking: Optional experimental read count stage can run with matched treated/control BAMs, but reports base fractions only. G/A allele-mixture caveat remains prominent in final prioritization.

## Code complete is not science validated

* Source code for external tools is provided, but the package does not include proprietary experimental inputs, ADAR validation data for KRAS G12D, Boltz weights, calibrated modification-specific Amber parameters, actual catalytic transition-state information or pancreatic tumor delivery evidence.
* Entire-protein predictions of ADAR1 with guide-target RNA can be unreliable. Use experimentally supported ADAR domains and verified structural templates, protein-protein cofactors and metal sites as indicated by a defined study design.
* RNA structure free energy and 100 ns RMSD do not predict editing yield. Caspase, interferon, off-target transcriptome effects, WT protein function and PDAC pathway adaptation require independent data.
* The biological model should be benchmarked on known ADAR substrates and a separate KRAS G12D cellular dataset before prioritizing synthesis.
* Protein restoration is not equivalent to anticancer efficacy, especially with KRAS secondary mutations and non-KRAS dependencies.

## Publication dataset provenance

The original six supplier-provided files are preserved. Publisher workbook parsing reports 104 constructs, 34 conservative standard-chemistry records and 70 complex or special records. This is a parser scope result, not an estimate of useful guides. Source Data assay indices have incomplete target-to-construct joins and are not assigned false ML labels.

## Independent acceptance sequence

1. Molecular identity: Verify reference genome, transcript, variant calls, phase, allele expression and mutant/WT dosage.
2. Chemistry: Confirm every base, covalent connection, chirality and force-field parameter against external references or QM/empirical checks.
3. Structural method: Benchmark ADAR protein-RNA topology and catalytic geometry against experimental controls, including Zn and RNA orientation.
4. Model validity: Demonstrate independent held-out target performance; prospective matched KRAS validation and phenotype required.
5. Delivery and efficacy: In-vitro PDAC cell and organoid activity, mutant allele suppression, viability, WT/bystander/off-target, immune activation and delivery; then animal studies only following safety and ethics approval.

All modifications to source were tested against synthetic datasets and CLI plans only. No Boltz-2 GPU inference, custom-residue MD simulation or true KRAS editing ML validation was performed in the packaging environment.
