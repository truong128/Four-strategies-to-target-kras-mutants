#!/usr/bin/env python3
"""Draw the method schematics used in the README and docs.

Every number shown here is a DESIGN PARAMETER taken from the shipped configuration files
(config.yaml, real_run_pdac_kras.yaml, ...). No figure is built from experimental or computed
results, so this script needs no data and can run anywhere:

    python docs/figures/make_figures.py          # writes PNG (300 dpi) and SVG next to this file
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

OUT = Path(__file__).resolve().parent

NAVY, TEAL, ORANGE, RED = "#1F3A5F", "#2A9D8F", "#E9A23B", "#C0392B"
PURPLE, GREY, LIGHT, WHITE = "#6C5B9E", "#7F8C8D", "#ECEFF1", "#FFFFFF"
SOFT_TEAL, SOFT_ORANGE, SOFT_PURPLE, SOFT_NAVY = "#D5EEEA", "#FBEBD0", "#E4DFF1", "#D9E2EF"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "savefig.dpi": 300,
    "svg.fonttype": "none", "axes.linewidth": 0.7,
})


def canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, text, fc=SOFT_NAVY, ec=NAVY, fs=8, tc="black", bold=False, lw=1.0, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", linespacing=1.25)


def arrow(ax, x0, y0, x1, y1, color=GREY, lw=1.3, style="-|>"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle=style, lw=lw, color=color, shrinkA=0, shrinkB=0))


def title(ax, text, sub=None, y=None):
    top = ax.get_ylim()[1] if y is None else y
    ax.text(1, top - 2.2, text, fontsize=13, fontweight="bold", color=NAVY, va="top")
    if sub:
        ax.text(1, top - 6.4, sub, fontsize=8.5, color=GREY, va="top")


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def chain(ax, y, items, x0=14, x1=97, h=11, fc=SOFT_NAVY, ec=NAVY, fs=7.4, gap=2.2, colors=None, start=None):
    n = len(items)
    w = (x1 - x0 - gap * (n - 1)) / n
    for i, text in enumerate(items):
        x = x0 + i * (w + gap)
        f, e = (colors[i] if colors else (fc, ec))
        box(ax, x, y, w, h, text, fc=f, ec=e, fs=fs)
        if start is not None:
            ax.text(x + 0.9, y + h - 1.3, str(start + i), fontsize=7.5, fontweight="bold", color=e, va="top")
        if i:
            arrow(ax, x - gap + 0.15, y + h / 2, x - 0.15, y + h / 2)


# ---------------------------------------------------------------------------------------------
def fig1_overview():
    fig, ax = canvas(10.5, 9.4)
    title(ax, "Four therapeutic strategies against KRAS-mutant pancreatic cancer",
          "One driver mutation, four computational routes. Method only: no data, no results.")
    box(ax, 1, 33, 12, 39, "KRAS\nG12D\n(c.35G>A)\n\nPDAC\ndriver", fc=NAVY, ec=NAVY, tc="white", bold=True, fs=10.5)

    lanes = [
        (72, "1  DRUG DISCOVERY    small molecules (Track A)", TEAL, SOFT_TEAL,
         ["Enamine REAL\n13.6B streamed,\nnever unpacked", "ML triage of the full\nlibrary, Vina redock\nof the top 50,000", "ADMET-AI, resistance\npanel, Boltz-2", "MD, MM/GBSA\n(+ OpenFE), retro,\nPareto top 50"]),
        (54, "2  DRUG REPURPOSING    approved drugs (Track A+)", ORANGE, SOFT_ORANGE,
         ["Table of approved\nsmall molecules", "Boltz-2 affinity on\nKRAS G12D (MSA)", "Keep affinity\nprobability above 0.5", "Re-score on the 17\nother KRAS alleles"]),
        (36, "3  ANTISENSE    ADAR RNA editing (Track B)", PURPLE, SOFT_PURPLE,
         ["KRAS-4B CDS +\nHGVS variant", "A-to-I correctability\n(G>A changes only)", "Guide design, bystander\nmap, chemistry transfer", "Bowtie2 off-target,\nViennaRNA, Boltz-2\nADAR1 geometry"]),
        (18, "4  PEPTIDE BINDERS    18 KRAS alleles (Track C)", RED, "#F6D9D5",
         ["7RPZ backbone,\n18 allele side chains", "Hotspot and H-bond\nconditions (YAML)", "RFD4-Proteina\n8-16 aa designs", "Optional AF3 /\nOpenFold3 refolding"]),
    ]
    for y, label, ec, fc, items in lanes:
        ax.text(15, y + 14.2, label, fontsize=9.6, fontweight="bold", color=ec)
        chain(ax, y, items, fc=fc, ec=ec, fs=8.6, x0=15, x1=99, h=12.5)
        arrow(ax, 13.2, 52.5, 14.7, y + 6, color=ec, lw=1.0)

    box(ax, 1, 1.5, 98, 9.5,
        "EXPERIMENTAL VALIDATION (outside this repository)\nbinding / editing  >  WT and RAS-isoform selectivity  >  pERK, G12D-selective killing  >  organoids + CAFs  >  PK  >  in vivo",
        fc=LIGHT, ec=GREY, fs=8.4, ls="--")
    save(fig, "fig1_overview")


# ---------------------------------------------------------------------------------------------
def fig2_funnel():
    fig, ax = canvas(14, 7.6)
    title(ax, "Track A funnel: design caps at every stage",
          "Caps come from config.yaml. They bound the compute at each step; they are not measured outcomes.")
    stages = [
        ("Enamine REAL library", "13.6 billion (streamed)", 92, NAVY),
        ("Property pre-filter (Enamine columns)", "MW 450, sLogP 5, TPSA 140, RotB 10", 84, NAVY),
        ("Representative training sample", "200,000 compounds", 70, TEAL),
        ("Vina labels (exhaustiveness 4)", "top 1% = virtual actives", 62, TEAL),
        ("Classifier screen of the full library", "keep top 50,000 by probability", 52, TEAL),
        ("Vina redock (exhaustiveness 8)", "50,000 compounds", 46, ORANGE),
        ("ADMET-AI hard filters", "hERG, AMES, hepatotoxicity < 0.5", 40, ORANGE),
        ("Diversity clustering", "Tanimoto cutoff 0.55", 34, ORANGE),
        ("Resistance panel + Boltz-2", "up to 200 compounds", 28, PURPLE),
        ("Retrosynthesis (RetroChimera + Retro*)", "50 compounds", 22, PURPLE),
        ("100 ns OpenMM MD + MM/GBSA (+ OpenFE)", "20 compounds; OpenFE optional", 16, PURPLE),
        ("Final ranking (Pareto + rank aggregation)", "50 reported", 11, RED),
    ]
    top, step = 81.0, 6.4
    cx = 50
    for i, (name, cap, width, col) in enumerate(stages):
        y = top - i * step
        hgt = step - 1.2
        ax.text(1, y + hgt / 2, f"{i + 1:>2}.  {name}", ha="left", va="center", fontsize=8.2, color="black")
        bw = 0.30 * width
        ax.add_patch(Rectangle((cx - bw / 2, y), bw, hgt, fc=col, ec="white", lw=1.0, alpha=0.95))
        ax.text(70, y + hgt / 2, cap, ha="left", va="center", fontsize=8.2, color="black")
    ax.text(70, 88.0, "Design cap", fontsize=9, fontweight="bold", color=NAVY)
    ax.text(1, 88.0, "Stage", fontsize=9, fontweight="bold", color=NAVY)
    ax.text(1, 2.5, "Colors: teal = machine-learning triage, orange = physical and ADMET filters, purple = shortlist-only expensive modeling.",
            fontsize=7.6, color=GREY)
    save(fig, "fig2_trackA_funnel")


def fig3_ml_design():
    fig, ax = canvas(15, 7.2)
    title(ax, "Track A machine-learning design",
          "Model selection never touches the locked test set. Scaffold-separated partitions replace temporal ones.")
    box(ax, 1, 62, 14, 20, "Vina scores on\n200,000 sampled\nmolecules\n\ntop 1% = active", fc=SOFT_TEAL, ec=TEAL, fs=8)
    box(ax, 19, 62, 16, 20, "Scaffold split\ntrain 0.65\nvalidation 0.15\ncalibration 0.10\nlocked test 0.10", fc=SOFT_NAVY, ec=NAVY, fs=8)
    box(ax, 39, 62, 18, 20, "20 classifiers\n16 classical (BRFSS set)\n+ FT-Transformer x2\n+ RealMLP\n+ Tabular ResNet", fc=SOFT_NAVY, ec=NAVY, fs=8)
    box(ax, 61, 62, 17, 20, "Select on validation\nbest single model\nvs top-3 soft vote\n(metric: average\nprecision)", fc=SOFT_NAVY, ec=NAVY, fs=8)
    box(ax, 82, 62, 16, 20, "Isotonic calibration\n+ class-conditional\nsplit conformal\n(epsilon 0.15)", fc=SOFT_ORANGE, ec=ORANGE, fs=8)
    for a, b in [(15, 19), (35, 39), (57, 61), (78, 82)]:
        arrow(ax, a + 0.2, 72, b - 0.2, 72)
    box(ax, 19, 30, 38, 18, "Features: ECFP4 (radius 2, 2048 bits) + 10 descriptors\nSMOTE intentionally NOT used on fingerprint bits\nEnrichment and recall reported because actives are about 1%", fc=LIGHT, ec=GREY, fs=8)
    box(ax, 61, 30, 17, 18, "ONE locked-test\nevaluation", fc=SOFT_PURPLE, ec=PURPLE, fs=8.5, bold=True)
    box(ax, 82, 30, 16, 18, "Chunked inference\n(9 production-eligible\nmodels, 4,000 per chunk,\n26 workers + 2 readers)", fc=SOFT_TEAL, ec=TEAL, fs=7.6)
    arrow(ax, 89.5, 62, 89.5, 48.2)
    arrow(ax, 78.2, 39, 81.8, 39)
    arrow(ax, 69.5, 62, 69.5, 48.2)
    ax.text(1, 14, "Slow benchmark models are row-capped and reported but cannot become production models.", fontsize=8, color=GREY)
    save(fig, "fig3_trackA_ml_design")


def fig4_resistance_panel():
    fig, ax = canvas(15, 8.2)
    title(ax, "Resistance-aware ensemble docking (shortlist stage)",
          "The same switch-II box is applied to every aligned receptor. Modeled variants are stress tests, never experimental structures.")
    ax.text(1, 82, "Receptor panel (config/pdb_panel.tsv)", fontsize=10, fontweight="bold", color=NAVY)
    groups = [
        ("Primary G12D", ["7RPZ\nGDP", "7T47\nGTP-like"], TEAL, SOFT_TEAL),
        ("WT control", ["6MBU\nWT GDP"], NAVY, SOFT_NAVY),
        ("Experimental resistance", ["9U8J\nY96D", "9U95\nQ99L", "9U99\nE62K", "9U97\nG13D"], ORANGE, SOFT_ORANGE),
        ("Modeled variants (labeled)", ["V9E", "V9W", "V9Q", "G13P", "T58Y", "R68G", "Y96W", "H95Q", "H95L"], PURPLE, SOFT_PURPLE),
    ]
    x = 1
    for name, items, ec, fc in groups:
        bw_, gp_ = 5.0, 0.8
        w = bw_ * len(items) + gp_ * (len(items) - 1)
        ax.text(x, 77.5, name, fontsize=8, fontweight="bold", color=ec)
        for i, it in enumerate(items):
            box(ax, x + i * (bw_ + gp_), 60, bw_, 14, it, fc=fc, ec=ec, fs=7.2 if len(items) < 9 else 6.8)
        x += w + 2.0
    ax.text(1, 52, "Score components (per compound)", fontsize=10, fontweight="bold", color=NAVY)
    comp = [
        ("G12D affinity", "mean over G12D receptors\n(GDP and GTP-like)", 0.45, TEAL),
        ("WT selectivity margin", "WT minus G12D docking score\n(kcal/mol), pass if >= 1.0", 0.25, NAVY),
        ("Resistance robustness", "worst resistance score minus\nG12D score, pass if <= 1.5", 0.30, PURPLE),
    ]
    x = 1
    for name, desc, wgt, col in comp:
        wd = 96 * wgt
        ax.add_patch(Rectangle((x, 36), wd - 0.6, 10, fc=col, ec="white"))
        ax.text(x + wd / 2 - 0.3, 41, f"{name}\nweight {wgt:.2f}", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
        ax.text(x + wd / 2 - 0.3, 31, desc, ha="center", va="top", fontsize=8, color="black")
        x += wd
    ax.text(1, 14, "Each component is converted to a percentile rank, then combined with the weights above; Vina scores are not binding free energies.\n"
                   "Caveat built into the method: a WT margin is only meaningful against a matched WT receptor (see docs/SCIENTIFIC_CAVEATS.md).",
            fontsize=8, color=GREY, va="top")
    save(fig, "fig4_trackA_resistance_panel")


def fig5_trackB():
    fig, ax = canvas(15, 8.4)
    title(ax, "Track B: from a variant table to ranked ADAR-recruiting guides",
          "Stages 1-13 run in order; each writes to a fresh output directory. Preflight refuses stale outputs and unverified inputs.")
    rows = [
        (66, "DESIGN (CPU only)", PURPLE, SOFT_PURPLE,
         ["Variant table\n(HGVS c.)", "Check reference base\nagainst validated CDS", "Is it G>A?\n(A-to-I reversible)", "Guide layouts\n33 / 40 / 59 nt,\norphan C opposite A", "Bystander A\nmap + codon\nconsequence"]),
        (40, "SPECIFICITY AND CHEMISTRY", TEAL, SOFT_TEAL,
         ["RESTORE 2.0 positional\nchemistry transfer\n(provenance-keyed)", "Bowtie2 local\nbinding vs sense\ntranscriptome", "Gene-aware filter\n(KRAS isoforms are\non target)", "ViennaRNA fold +\ncofold of guide\nand target", ""]),
        (14, "STRUCTURE AND TRIAGE (GPU)", ORANGE, SOFT_ORANGE,
         ["Boltz-2: ADAR1\n716-1226 + RNA duplex\n(geometry only)", "Optional MD\n(needs validated\nmodified-residue\nparameters)", "Chemistry-aware ML\n(abstains without\ncurated labels)", "Experiment-gated\nresearch ranking", ""]),
    ]
    counter = 1
    for y, label, ec, fc, items in rows:
        ax.text(1, y + 17.5, label, fontsize=9.5, fontweight="bold", color=ec)
        items = [i for i in items if i]
        chain(ax, y, items, x0=1, x1=1 + 96 * len(items) / 5, h=14, fc=fc, ec=ec, fs=7.6, start=counter)
        counter += len(items)
    ax.text(1, 3.2, "Hard boundaries: Boltz geometry is not catalytic evidence, ViennaRNA is not modification-aware, and no stage estimates editing efficiency.",
            fontsize=8, color=GREY)
    save(fig, "fig5_trackB_flow")


def fig6_adar_concept():
    fig, ax = canvas(15, 8.6)
    title(ax, "Which KRAS variants can ADAR correct?",
          "ADAR deaminates adenosine to inosine, read as guanosine. Only G>A driver changes can be reversed directly.")
    hdr = ["Variant", "cDNA change", "Direct A-to-I reversal"]
    rows = [
        ("G12D", "c.35G>A", "Yes", TEAL), ("G12S", "c.34G>A", "Yes", TEAL), ("G13D", "c.38G>A", "Yes", TEAL),
        ("G12V", "c.35G>T", "No", RED), ("G12C", "c.34G>T", "No", RED), ("G12R", "c.34G>C", "No", RED),
        ("Q61H", "c.183A>C / c.183A>T", "No", RED),
    ]
    x = [1, 12, 33]
    for i, h in enumerate(hdr):
        ax.text(x[i], 80, h, fontsize=9, fontweight="bold", color=NAVY)
    for r, (v, c, ok, col) in enumerate(rows):
        y = 73 - r * 6.6
        ax.add_patch(Rectangle((0.5, y - 2.2), 46, 5.8, fc=LIGHT if r % 2 == 0 else WHITE, ec="none"))
        ax.text(x[0], y, v, fontsize=9, va="center")
        ax.text(x[1], y, c, fontsize=9, va="center")
        ax.text(x[2], y, ok, fontsize=9, va="center", color=col, fontweight="bold")
    ax.text(1, 21, "The pipeline reports non-correctable variants explicitly\ninstead of skipping them silently.", fontsize=8, color=GREY, va="top")

    ax.text(52, 80, "Guide design principle (schematic)", fontsize=10, fontweight="bold", color=NAVY)
    target = list("GGUAGUUGGAGCUGUAGGCGUAGGCAAGAG")
    n = len(target)
    x0, dx = 52, 1.5
    ax.text(x0 - 1, 66, "mRNA 5'", fontsize=7.5, ha="right", va="center", color=GREY)
    ax.text(x0 + n * dx + 0.2, 66, "3'", fontsize=7.5, va="center", color=GREY)
    ax.text(x0 - 1, 56, "guide 3'", fontsize=7.5, ha="right", va="center", color=GREY)
    ax.text(x0 + n * dx + 0.2, 56, "5'", fontsize=7.5, va="center", color=GREY)
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    target_idx = 15
    for i, base in enumerate(target):
        xi = x0 + i * dx + dx / 2
        hot = i == target_idx
        ax.text(xi, 66, base, ha="center", va="center", fontsize=8.5, fontweight="bold" if hot else "normal",
                color=RED if hot else "black", family="DejaVu Sans Mono")
        pair = "C" if hot else comp[base]
        ax.text(xi, 56, pair, ha="center", va="center", fontsize=8.5, fontweight="bold" if hot else "normal",
                color=RED if hot else NAVY, family="DejaVu Sans Mono")
        ax.plot([xi, xi], [58.2, 63.8], color=RED if hot else GREY, lw=1.0, ls=":" if hot else "-")
    ax.annotate("target adenosine\n(edited A to I)", xy=(x0 + target_idx * dx + dx / 2, 68), xytext=(x0 + 38, 78),
                fontsize=8, color=RED, ha="center", arrowprops=dict(arrowstyle="-|>", color=RED, lw=1))
    ax.annotate("orphan C opposite A\n(A-C mismatch recruits ADAR)", xy=(x0 + target_idx * dx + dx / 2, 54), xytext=(x0 + 12, 42),
                fontsize=8, color=RED, ha="center", arrowprops=dict(arrowstyle="-|>", color=RED, lw=1))
    ax.text(52, 28, "Illustrative sequence only. Real guides come from the validated KRAS CDS.\n"
                    "Layouts: 33, 40 and 59 nt (RESTORE 2.0 constructs). Every other A in the\n"
                    "footprint is a bystander and is mapped for codon consequence.", fontsize=8, color=GREY, va="top")
    save(fig, "fig6_adar_principle")


def fig7_trackC():
    fig, ax = canvas(15, 6.4)
    title(ax, "Track C: short peptide binders for 18 KRAS allele targets",
          "Backbone from 7RPZ chain A (residues 1-169). Mutant side chains are grafted, so each target carries one substitution.")
    chain(ax, 56, ["PDB 7RPZ\n(G12D, MRTX1133\nremoved)", "Graft side chain\nat codon 12, 13, 61,\n117 or 146", "Condition spec:\nhotspots, H-bond\ndonor / acceptor", "RFD4-Proteina\n(RFdiffusion3 /\nProteina)", "Complexes with\n8-16 aa peptide", "Optional refolding\nmetrics (AF3 /\nOpenFold3 oracle)"],
          x0=1, x1=98, h=22, fc=SOFT_ORANGE, ec=ORANGE, fs=8)
    alleles = "G12D  WT  G12C  G12V  G12R  G12A  G13D  A146T  G12S  G12F  G13C  G13R  Q61H  Q61L  Q61R  Q61K  A146V  K117N"
    ax.text(1, 40, "Allele panel", fontsize=9, fontweight="bold", color=NAVY)
    ax.text(1, 34, alleles, fontsize=8.5, family="DejaVu Sans Mono")
    ax.text(1, 24, "Codons 12, 13 and 61 use the P-loop / switch-II pocket. K117 and A146 use the local nucleotide-site residues.\n"
                   "rfd4-proteina is an external dependency (cloned by scripts/bootstrap.sh c), not vendored here.", fontsize=8, color=GREY, va="top")
    save(fig, "fig7_trackC_peptides")


def fig8_evidence_ladder():
    fig, ax = canvas(15, 8.2)
    title(ax, "Claim boundaries: what a computed result does and does not mean",
          "This repository stops at the first rung. Every higher rung needs a wet-lab experiment.")
    mol = ["Computed prioritization\n(this repository)", "Biochemical binding, GTP\nexchange inhibition", "WT KRAS / HRAS / NRAS\nselectivity", "pERK suppression",
           "G12D-selective cell killing", "Resistance-variant cells", "3D PDAC organoids + CAFs", "PK and tumor exposure", "In vivo efficacy and safety"]
    rna = ["Computed guide ranking\n(this repository)", "Synthesis + in vitro editing\nof KRAS RNA", "Editing in G12D cell lines\n(on-target, bystanders)", "Transcriptome-wide\noff-target editing",
           "Protein restoration,\npERK suppression", "Selective killing of\nmutant cells", "Organoids + CAFs", "Delivery (LNP / AAV)", "In vivo efficacy and safety"]
    n = len(mol)
    for i in range(n):
        y = 77 - i * 8.2
        fc_m = SOFT_TEAL if i == 0 else LIGHT
        ec_m = TEAL if i == 0 else GREY
        fc_r = SOFT_PURPLE if i == 0 else LIGHT
        ec_r = PURPLE if i == 0 else GREY
        off = i * 1.2
        box(ax, 2 + off, y, 36, 6.8, mol[i], fc=fc_m, ec=ec_m, fs=7.6, bold=(i == 0), ls="-" if i == 0 else "--")
        box(ax, 54 + off, y, 36, 6.8, rna[i], fc=fc_r, ec=ec_r, fs=7.6, bold=(i == 0), ls="-" if i == 0 else "--")
    ax.text(2, 86.5, "Small molecules and peptides", fontsize=9.5, fontweight="bold", color=TEAL)
    ax.text(54, 86.5, "ADAR RNA editing", fontsize=9.5, fontweight="bold", color=PURPLE)
    save(fig, "fig8_evidence_ladder")


if __name__ == "__main__":
    for fn in (fig1_overview, fig2_funnel, fig3_ml_design, fig4_resistance_panel, fig5_trackB,
               fig6_adar_concept, fig7_trackC, fig8_evidence_ladder):
        fn()
