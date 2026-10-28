# Scientific caveats

These are properties of the *method*. They are repeated in outputs where relevant.

## What a score means

| Quantity | It is | It is not |
|---|---|---|
| Vina docking score | a fast rank heuristic | a binding free energy |
| ML probability (Track A) | calibrated likelihood of being a top-1% Vina scorer on a scaffold-held-out test | likelihood of real binding |
| Boltz-2 affinity probability / value | a model prediction of binder likelihood and log10 IC50 (uM) | a measurement; it was not trained on your chemotypes or alleles |
| OpenMM ligand RMSD | pose stability over the simulated window | affinity or residence time |
| MM/GBSA | a single, noisy end-point estimate (generalized Born, not PB) | a free-energy difference you can compare across series |
| OpenFE absolute binding free energy | an alchemical estimate with an uncertainty, only as good as the pose, force field, cofactors and sampling | an experimental affinity; it also needs protonation and cofactors to match the real pocket |
| RetroChimera route | a model suggestion that can hallucinate reactions | a synthesis plan; a chemist must check it |
| ADMET-AI | a screening-level prediction | tissue exposure; PDAC stroma makes PK a first-class question |
| ViennaRNA energy | unmodified RNA thermodynamics | folding of chemically modified guides |
| Boltz-2 ADAR1 geometry | a plausible duplex pose with canonical RNA | catalytic mechanism; Zn2+ and modifications are absent |
| Bowtie2 hit | sequence complementarity | demonstrated off-target editing |
| Chemistry transfer | RESTORE 2.0 patterns mapped by position with provenance | validated on KRAS |

## Track A

1. **Wild-type counter-screen.** If the WT receptor is an apo structure or a different construct, the WT margin is confounded by conformation. Use a matched receptor (PDB 9O0N, or a D12G back-mutation of 7RPZ).
2. **Modeled mutants** are single substitutions modeled on 7RPZ with a quick side-chain relaxation. They are stress tests, labeled `MODEL_*`.
3. **Covalent chemistry.** Vina cannot form the G12C covalent bond; poses are noncovalent. The method also does not model RAS(ON) cyclophilin-A tri-complex inhibitors.
4. **Bypass biology** (EGFR/FGFR feedback, PI3K, YAP/MYC, EMT, autophagy, CAFs) is out of scope for a binding-site screen.
5. **Virtual actives** are defined by Vina rank, so the ML learns to imitate Vina, not biology. Enrichment on the locked test is a statement about that proxy.
6. **Library prefilters** use Enamine's precomputed columns; compounds with missing columns are handled as documented in the scripts.

## Track B

1. **Only G>A variants are directly reversible** (G12D, G12S, G13D). Others are reported as not correctable.
2. **Selectivity.** A guide C opposite a mutant A forms a perfect G-C pair with wild-type G, so mutant-versus-wild-type selectivity is not guaranteed and not computed. Every candidate needs a WT-transcript counter-assessment.
3. **Quantifying editing.** Edited mutant mRNA looks like WT in A/G counting; heterozygous WT expression confounds the fraction. Use phasing and allele-specific assays.
4. **No stage predicts editing efficiency.** The ML stage abstains without curated, independent, KRAS-relevant labels.
5. **ADAR isoform matters.** ADAR2 structural data do not validate ADAR1 activity in PDAC cells.

## Track C

Designs are generated against a rigid single-chain target with one substitution. Refolding metrics are off by default. Peptides need binding assays long before they need anything else.

## The rule

Computed does not mean validated. Figure 8 lists the experiments that convert each hypothesis into evidence.
