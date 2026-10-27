#!/usr/bin/env python3
"""Generate the method tables (CSV + one Markdown file) from the shipped configuration.

All content is method metadata: stages, parameters, tools, access rules. The receptor and allele
panels are read straight from the TSV files in tracks/. No result or dataset is used.

    python docs/make_tables.py       # writes docs/tables/T*.csv and docs/tables/TABLES.md
"""
from __future__ import annotations

import csv
from pathlib import Path

import yaml

DOCS = Path(__file__).resolve().parent
REPO = DOCS.parent
OUT = DOCS / "tables"
OUT.mkdir(exist_ok=True)
A = REPO / "tracks" / "a_small_molecule_discovery"
B = REPO / "tracks" / "b_adar_rna_editing"

TABLES: list[tuple[str, str, list[str], list[list[str]], str]] = []


def add(tid: str, title: str, header: list[str], rows: list[list[str]], note: str = "") -> None:
    TABLES.append((tid, title, header, rows, note))


def read_tsv(path: Path) -> tuple[list[str], list[list[str]]]:
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


cfg = yaml.safe_load((A / "config.yaml").read_text())

# ------------------------------------------------------------------------------------------------
add("T1", "Tracks, purpose and entry points",
    ["Track", "Modality", "Question it answers", "Engine", "Entry point"],
    [["A", "Small molecules", "Which purchasable molecules in a 13.6 billion library are worth buying for KRAS G12D?",
      "Nextflow DSL2 + Python", "scripts/run_track_a.sh screen --library PATH"],
     ["A2", "Small molecules (light)", "Rank up to about 1e5 molecules against a resistance-aware receptor ensemble in one command",
      "Python scripts 00 to 09", "scripts/run_track_a2.sh"],
     ["B", "ADAR RNA editing", "Which ADAR-recruiting guide RNAs can reverse a G>A KRAS driver, and what do they hit elsewhere?",
      "Python + Bowtie2 + ViennaRNA + Boltz-2", "scripts/run_track_b.sh plan | full"],
     ["C", "Peptides", "Which short peptides does a generative model propose for each of 18 KRAS allele targets?",
      "RFD4-Proteina (external)", "scripts/run_track_c.sh design"]])

stage_rows = [
    ["1", "Target and box", "fetch_structures.py", "PDBFixer, Meeko", "RCSB 7RPZ", "receptor PDBQT, box", "box center (1.714, 4.927, -23.164), size 24 A"],
    ["2", "Library streaming and filter", "prepare_library.py / sample_library.py", "RDKit, lbzip2", "Enamine .cxsmiles.bz2 shards", "standardized or sampled SMILES",
     "MW <= %s, sLogP <= %s, TPSA <= %s, RotB <= %s" % (cfg["library"]["max_mw"], cfg["library"]["max_clogp"], cfg["library"]["max_tpsa"], cfg["library"]["max_rotatable_bonds"])],
    ["3", "Training sample", "sample_library.py", "reservoir sampling", "library", "training_sample.smi", "n = %d, seed %d" % (cfg["library"]["training_sample_n"], cfg["library"]["reservoir_seed"])],
    ["4", "Training docking", "dock_vina_batch.py", "AutoDock Vina 1.2", "sample + receptor", "training_docking.csv", "exhaustiveness %d, %d workers" % (cfg["docking"]["exhaustiveness"], cfg["docking"]["workers"])],
    ["5", "Classifier benchmark", "benchmark_16_models.py, tabular_nets.py", "scikit-learn, PyTorch", "docking labels", "model bundle, comparison, locked-test metrics",
     "top %s%% = active; scaffold split; 20 models" % cfg["docking"]["active_percent"]],
    ["6", "Ultra-large inference", "predict_ultralarge.py", "joblib, RDKit", "library + bundle", "ml_hits.csv.gz", "keep top %d; conformal epsilon %s" % (cfg["ml"]["screen_top_k"], cfg["ml"]["conformal_epsilon"])],
    ["7", "Redocking", "dock_vina_batch.py --redock", "AutoDock Vina", "ML hits", "redocked_hits.csv", "exhaustiveness %d" % cfg["docking"]["redock_exhaustiveness"]],
    ["8", "ADMET filter", "run_admet.py", "ADMET-AI 2", "redocked hits", "admet_filtered.csv", "hERG, AMES, hepatotoxicity < 0.5"],
    ["9", "Diversity clustering", "cluster_hits.py", "RDKit Butina-style", "filtered hits", "diverse_hits.csv", "Tanimoto cutoff %s" % cfg["clustering"]["tanimoto_cutoff"]],
    ["10", "Resistance panel", "build_resistance_panel.py, dock_resistance_panel.py", "Biopython, PDBFixer, Meeko, Vina", "diverse hits + panel TSV", "resistance_scores.csv",
     "up to %d compounds; 7 experimental + 9 modeled receptors" % cfg["resistance"]["max_input"]],
    ["11", "Boltz-2 complex and affinity", "run_boltz_batch.py", "Boltz-2", "short list + protein", "boltz_scores.csv", "up to %d compounds; no custom kernels on Blackwell" % cfg["boltz"]["max_input"]],
    ["12", "Retrosynthesis", "synthetic_feasibility.py", "RetroChimera + Syntheseus Retro*", "best docked", "synthesis_scores.csv", "%d compounds, USPTO-50K checkpoint" % cfg["retrosynthesis"]["max_input"]],
    ["13", "Molecular dynamics", "run_md_batch.py, md_openmm.py", "OpenMM, amber14, OpenFF 2.2.1", "top docked + protein", "md_scores.csv", "%d compounds, %d ns" % (cfg["md"]["max_short"], cfg["md"]["short_ns"])],
    ["14", "MM/GBSA", "run_mmpbsa_batch.py, mmgbsa_from_md.py", "AmberTools", "MD trajectories", "mmpbsa_scores.csv", "igb=5, 0.15 M salt; %d compounds" % cfg["mmpbsa"]["max_input"]],
    ["15", "Final ranking", "final_rank.py", "pandas", "all tables", "ranked and Pareto CSVs", "13 weighted criteria + Pareto front; report %d" % cfg["ranking"]["final_n"]],
]
add("T2", "Track A stages (Nextflow processes in main.nf)",
    ["#", "Stage", "Script", "Tool", "Input", "Output", "Key parameter"], stage_rows,
    "Parameters are defaults from config.yaml; they cap compute and are not outcomes.")

add("T3", "Track A ranking weights (config.yaml, section ranking)",
    ["Criterion", "Weight", "Direction"],
    [[k, str(v), d] for (k, v), d in zip(cfg["ranking"]["weights"].items(),
     ["higher", "lower", "higher", "higher", "lower", "higher", "lower", "higher", "lower", "higher", "lower", "higher", "higher"])],
    "Each criterion is converted to a percentile rank before weighting. Missing criteria are skipped, not imputed.")

h, rows = read_tsv(A / "config" / "pdb_panel.tsv")
add("T4", "Resistance receptor panel (config/pdb_panel.tsv)", h, rows,
    "MODEL_* rows are single substitutions modeled on 7RPZ and are stress tests, not experimental structures.")

h, rows = read_tsv(A / "config" / "kras_allele_panel.tsv")
add("T5", "KRAS allele panel (config/kras_allele_panel.tsv)", h, rows,
    "Alleles without a deposited GDP structure with an intact switch II are built as one substitution on wild-type 6MBU.")

add("T6", "Track B stages",
    ["#", "Stage", "Executor", "Input", "Output", "Boundary"],
    [["1", "Reference setup", "scripts/setup_real_run.py", "Ensembl GTF + cDNA, UniProt P55265", "KRAS CDS, gene map, ADAR1 FASTAs", "CDS validated: ATG start, one stop, codon 12 = Gly"],
     ["2", "Variant check", "scripts/adar_pipeline.py", "HGVS c. table", "reversibility call", "reference base must match the CDS"],
     ["3", "Guide design", "scripts/adar_pipeline.py", "correctable variants", "guides at 33, 40, 59 nt", "A-C mismatch design; no editing-rate estimate"],
     ["4", "Chemistry transfer", "scripts/transfer_publisher_chemistry.py", "RESTORE 2.0 workbooks (user downloaded)", "positional chemistry per guide", "transferred, not validated on KRAS"],
     ["5", "Local RNA structure", "scripts/advanced.py", "guide and target", "ViennaRNA fold and cofold", "not modification-aware"],
     ["6", "Off-target screen", "full/offtarget.py", "transcriptome + guides", "Bowtie2 local hits, per gene", "binding opportunity, not editing"],
     ["7", "Structure", "full/structure.py", "ADAR1 716-1226 + RNA", "Boltz-2 ensemble geometry", "canonical RNA, no Zn2+, no modifications"],
     ["8", "Parameterization and MD", "full/parameterize.py, full/dynamics.py", "validated Amber parameters", "NPT trajectories", "refused without hash-bound validation"],
     ["9", "Chemistry-aware ML", "full/learning.py", "curated measured editing", "grouped CV, external holdout", "abstains without labeled data"],
     ["10", "Experimental read-out", "full/experimental.py", "treated and control BAM", "base fractions at sites", "heterozygous allele warnings"],
     ["11", "Ranking", "full/prioritize.py", "all evidence", "research ranking", "experiment-gated"]],
    "Every run needs a NEW output directory; existing outputs are refused to prevent stale-file reuse.")

add("T7", "Which KRAS variants ADAR can reverse directly",
    ["Variant", "cDNA change", "Direct A-to-I reversal", "Reason"],
    [["G12D", "c.35G>A", "Yes", "mutant A is the editing target"],
     ["G12S", "c.34G>A", "Yes", "mutant A is the editing target"],
     ["G13D", "c.38G>A", "Yes", "mutant A is the editing target"],
     ["G12V", "c.35G>T", "No", "G>T cannot be undone by A-to-I"],
     ["G12C", "c.34G>T", "No", "G>T cannot be undone by A-to-I"],
     ["G12R", "c.34G>C", "No", "G>C cannot be undone by A-to-I"],
     ["Q61H", "c.183A>C / c.183A>T", "No", "mutant base is not A"]],
    "Non-correctable variants are reported as such, never silently skipped. Mutant versus wild-type selectivity is not established by the pipeline.")

add("T8", "Track C allele targets and pocket definitions",
    ["Codon", "Alleles", "Pocket", "Hotspot rule"],
    [["12", "G12D, G12C, G12V, G12R, G12A, G12S, G12F (+ WT)", "P-loop / switch II", "mutant side chain + two nearest side chains"],
     ["13", "G13D, G13C, G13R", "P-loop / switch II", "mutant side chain + two nearest side chains"],
     ["61", "Q61H, Q61L, Q61R, Q61K", "P-loop / switch II", "mutant side chain + two nearest side chains"],
     ["117", "K117N", "G4 motif, nucleotide base side", "local nucleotide-site residues"],
     ["146", "A146T, A146V", "G5 motif, nucleotide cleft", "residues pointing into the cleft"]],
    "Backbone: 7RPZ chain A residues 3 to 169 (residue 2 has no coordinates). Designed chain: 8 to 16 residues.")

add("T9", "Software stack (versions used during development; pin your own with the environment files)",
    ["Component", "Version tested", "Used by", "Environment"],
    [["Python", "3.11 / 3.12", "all tracks", "pdac-dd, pdac-md, kras-adar"],
     ["Nextflow", "26.04", "Track A", "system"],
     ["RDKit", "2026.03", "A, A2", "pdac-dd, pdac-md"],
     ["AutoDock Vina", "1.2.7", "A, A2", "pdac-dd"],
     ["Meeko", "0.8.0", "A, A2", "pdac-dd"],
     ["scikit-learn / XGBoost / LightGBM / CatBoost", "1.9 / 3.4 / 4.7 / 1.2", "A (classifiers)", "pdac-dd"],
     ["PyTorch", "2.11 (CUDA 12.8 wheels)", "A (neural nets, Boltz-2, RetroChimera)", "pdac-dd, boltz2, retrochimera"],
     ["ADMET-AI", "2.0", "A, A2", "admet-ai"],
     ["Boltz-2", "2.2", "A, B, extensions", "boltz2, kras-adar"],
     ["OpenMM / PDBFixer", "8.6 / 1.12", "A (MD, receptors)", "pdac-md"],
     ["OpenFF toolkit / forcefields", "0.18 / 2026.01", "A (ligand parameters)", "pdac-md"],
     ["AmberTools", "26.0", "A (MM/GBSA)", "pdac-md or amber"],
     ["RetroChimera / Syntheseus", "0.1 dev / 0.9", "A (retrosynthesis)", "retrochimera"],
     ["ViennaRNA", "2.4.7", "B", "kras-adar"],
     ["Bowtie2", "2.5.5", "B", "kras-adar"],
     ["samtools / bcftools", "current conda-forge", "B", "kras-adar"],
     ["RFD4-Proteina", "upstream main", "C", "external/rfd4-proteina (.venv)"],
     ["GNINA", "1.3", "A (allele CNN rescoring, optional)", "binary on PATH or GNINA_EXE"]])

add("T10", "Inputs that are NOT in this repository and how to obtain them",
    ["Resource", "Why excluded", "How to obtain", "Used by"],
    [["Enamine REAL (13.6 B)", "vendor license, about 330 GB compressed", "register at enamine.net and download the .cxsmiles.bz2 shards", "A"],
     ["Purchasable building-block inventory", "third-party data", "Syntheseus Retro* benchmark set (Figshare 10.6084/m9.figshare.25376728)", "A (retrosynthesis)"],
     ["PDB structures (7RPZ, 6MBU, ...)", "public, regenerated", "scripts/fetch_public_inputs.sh (RCSB)", "A, A2, C"],
     ["Ensembl GRCh38 release 113 cDNA + GTF", "public, large", "scripts/fetch_public_inputs.sh b", "B"],
     ["UniProt P55265 (ADAR1)", "public, regenerated", "downloaded by setup_real_run.py", "B"],
     ["RESTORE 2.0 supplementary workbooks", "publisher copyright", "download from the publisher, ingest with scripts/ingest_publisher_workbooks.py", "B"],
     ["Boltz-2 weights", "large, licensed by Boltz", "downloaded by Boltz on first use", "A, B"],
     ["RetroChimera checkpoints", "large, third party", "see upstream RetroChimera instructions", "A"],
     ["RFD4-Proteina weights", "large, third party", "scripts/bootstrap.sh c (Hugging Face release, about 1.7 GB)", "C"],
     ["FDA drug SMILES tables", "user supplied", "FDA_DRUG_PKL / FDA_DRUG_CSV environment variables", "A (repurposing scripts)"],
     ["Tumor-normal or RNA-seq data", "patient or experimental data", "your own; never commit", "B (optional)"]])

add("T11", "Compute guide (what drives cost, not timings)",
    ["Stage group", "Resource", "Scales with", "Config knob"],
    [["Streaming + sampling", "CPU cores, disk read", "library size", "sample_file_parallel, sample_decompress_threads"],
     ["Vina docking", "CPU cores", "molecules x exhaustiveness", "docking.workers, exhaustiveness"],
     ["Classifier training", "CPU, optional GPU for neural nets", "rows x models", "ml.*_max_rows, dl_epochs"],
     ["Library inference", "CPU cores, RAM", "library size", "inference_workers, inference_chunk_size"],
     ["Boltz-2", "GPU", "compounds x samples", "boltz.max_input, diffusion_samples"],
     ["OpenMM MD", "GPU", "compounds x ns", "md.max_short, short_ns"],
     ["Track B off-target", "CPU cores, RAM", "guides x transcriptome", "offtarget.threads"],
     ["Track C design", "GPU", "alleles x replicates", "N, BATCH, ALLELES"]])

add("T12", "Known limitations built into the method",
    ["Limitation", "Consequence", "Mitigation in the repository"],
    [["Docking scores are not binding free energies", "ranks are heuristics", "rank-based scoring; MD and MM/GBSA on a short list only"],
     ["WT counter-screen receptor can be apo", "WT margin is confounded by conformation", "documented; use a matched WT receptor (9O0N or a D12G back-mutation of 7RPZ)"],
     ["Modeled mutants are unrelaxed single substitutions", "stress tests only", "labeled MODEL_* everywhere, never reported as experimental"],
     ["Vina cannot model covalent chemistry (G12C)", "pose is noncovalent", "documented in the allele panel notes"],
     ["Boltz-2 geometry is not catalysis", "no editing claim from structure", "Track B reports geometry only"],
     ["ViennaRNA is not modification-aware", "chemistry effects on folding are unknown", "labeled in outputs"],
     ["Chemistry transfer is from RESTORE 2.0, not KRAS", "unvalidated on KRAS", "provenance keys on every transferred pattern"],
     ["Editing ML needs curated labels", "no predictor without data", "the stage abstains explicitly"],
     ["Everything is computed", "no biological claim", "Figure 8 and docs/SCIENTIFIC_CAVEATS.md"]])


def md_table(header: list[str], rows: list[list[str]]) -> str:
    esc = lambda s: str(s).replace("|", "\\|")
    lines = ["| " + " | ".join(esc(h) for h in header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows]
    return "\n".join(lines)


def main() -> None:
    md = ["# Tables", "", "Generated by `docs/make_tables.py`. Method metadata only; no results.", ""]
    for tid, title, header, rows, note in TABLES:
        with open(OUT / f"{tid}.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(header)
            w.writerows(rows)
        md += [f"## {tid}. {title}", "", md_table(header, rows), ""]
        if note:
            md += [f"*{note}*", ""]
    (OUT / "TABLES.md").write_text("\n".join(md))
    print(f"wrote {len(TABLES)} tables to {OUT}")


if __name__ == "__main__":
    main()
