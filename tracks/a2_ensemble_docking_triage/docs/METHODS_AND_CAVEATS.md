# Methods and caveats

## Primary scientific objective
Optimize **KRAS G12D selectivity + resistance robustness + developability**, rather than docking score alone.

## Structural logic
- **7RPZ** anchors the crystallographic MRTX1133 pose and GDP-state pocket.
- **7T47** tests an active-like nucleotide state.
- **6MBU** is a WT counter-screen.
- **9U8J / 9U95 / 9U99 / 9U97** provide experimentally observed secondary-mutant structural contexts.
- Additional validated resistance variants are modeled from 7RPZ and treated as uncertainty-bearing stress tests.

## Resistance liabilities included
Validated 2024 MRTX1133 resistance screen: V9E, V9W, V9Q, G13P, T58Y, R68G, Y96W, Q99L. H95/Y96 substitutions are included because H95 is a KRAS-paralog selectivity handle and Y96 substitutions can directly disrupt MRTX1133 activity.

## Critical caveats
1. A modeled mutant is not evidence that the mutant will adopt that exact structure.
2. Vina scores are not delta-G values and should not be compared as if they were experimental affinities.
3. WT counter-screen docking is only an early selectivity filter; biochemical/paralog experiments are mandatory.
4. PDAC resistance is not only on-target. EGFR/FGFR feedback, PI3K-AKT-mTOR, YAP/MYC, EMT/lineage plasticity, autophagy/macropinocytosis, CAFs and drug exposure all matter.
5. The dense PDAC microenvironment makes tumor exposure a first-class property; cell-free affinity cannot substitute for tissue PK.
6. This workflow is for the noncovalent Switch-II-pocket route. Cyclophilin-A ternary/covalent RAS(ON) inhibitors need a different physical model.
