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

NAVY, TEAL, ORANGE, RED = "#1F3A5F", "#2A9D8F", "#D98E1F", "#C0392B"
PURPLE, GREY, LIGHT, WHITE = "#6C5B9E", "#5F6B73", "#ECEFF1", "#FFFFFF"
SOFT_TEAL, SOFT_ORANGE, SOFT_PURPLE, SOFT_NAVY = "#D5EEEA", "#FBEBD0", "#E4DFF1", "#D9E2EF"
SOFT_RED = "#F6D9D5"
COLORS = {"teal": (SOFT_TEAL, TEAL), "orange": (SOFT_ORANGE, ORANGE), "purple": (SOFT_PURPLE, PURPLE),
          "red": (SOFT_RED, RED), "navy": (SOFT_NAVY, NAVY), "grey": (LIGHT, GREY)}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "savefig.dpi": 300,
    "svg.fonttype": "none", "axes.linewidth": 0.7,
})


def canvas(w, h):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")
    return fig, ax


def box(ax, x, y, w, h, text, fc=SOFT_NAVY, ec=NAVY, fs=10, tc="black", bold=False, lw=1.2, ls="-", align="center"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.0,rounding_size=0.8",
                                fc=fc, ec=ec, lw=lw, ls=ls))
    if align == "center":
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc,
                fontweight="bold" if bold else "normal", linespacing=1.3)
    else:
        ax.text(x + 1.6, y + h / 2, text, ha="left", va="center", fontsize=fs, color=tc,
                fontweight="bold" if bold else "normal", linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, color=GREY, lw=1.3, style="-|>"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle=style, lw=lw, color=color, shrinkA=0, shrinkB=0))


def title(ax, text, sub=None):
    ax.text(1, 98.5, text, fontsize=15, fontweight="bold", color=NAVY, va="top")
    if sub:
        ax.text(1, 93.2, sub, fontsize=10, color=GREY, va="top")


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def bands(ax, rows, y_top, pitch, h, left_w=27, x0=1, x1=99, numbered=True, arrows=True):
    """Vertical flow, one band per step. rows = [(head, body, color_key)]."""
    for i, (head, body, key) in enumerate(rows):
        y = y_top - i * pitch
        fc, ec = COLORS[key]
        label = f"{i + 1}   {head}" if numbered else head
        box(ax, x0, y, left_w, h, label, fc=ec, ec=ec, tc="white", bold=True, fs=11, align="left")
        box(ax, x0 + left_w + 1, y, x1 - x0 - left_w - 1, h, body, fc=fc, ec=ec, fs=10, align="left")
        if arrows and i:
            arrow(ax, x0 + left_w / 2, y + pitch - 0.3, x0 + left_w / 2, y + h + 0.3)


def swatches(ax, items, y):
    """Legend of (label, color_key) pairs."""
    x = 1
    for label, key in items:
        fc, ec = COLORS[key]
        ax.add_patch(FancyBboxPatch((x, y), 2.4, 3.0, boxstyle="round,pad=0,rounding_size=0.4", fc=fc, ec=ec, lw=1.2))
        ax.text(x + 3.4, y + 1.5, label, fontsize=10, va="center")
        x += 5.5 + 0.95 * len(label)


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
    title(ax, "Four strategies against KRAS-mutant pancreatic cancer",
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
    fig, ax = canvas(10.5, 8.0)
    title(ax, "Track A funnel: design caps at every stage",
          "Caps come from config.yaml. They bound the compute at each step; they are not measured outcomes.")
    stages = [
        ("Enamine REAL library", "13.6 billion (streamed)", 92, NAVY),
        ("Property pre-filter", "MW 450, sLogP 5,\nTPSA 140, RotB 10", 84, NAVY),
        ("Representative training sample", "200,000 compounds", 70, TEAL),
        ("Vina labels (exhaustiveness 4)", "top 1% = virtual actives", 62, TEAL),
        ("Classifier screen, full library", "keep top 50,000", 52, TEAL),
        ("Vina redock (exhaustiveness 8)", "50,000 compounds", 46, ORANGE),
        ("ADMET-AI hard filters", "hERG, AMES, hepato-\ntoxicity < 0.5", 40, ORANGE),
        ("Diversity clustering", "Tanimoto cutoff 0.55", 34, ORANGE),
        ("Resistance panel + Boltz-2", "up to 200 compounds", 28, PURPLE),
        ("Retrosynthesis", "50 compounds", 22, PURPLE),
        ("MD + MM/GBSA (+ OpenFE)", "20 compounds, 100 ns;\nOpenFE optional", 16, PURPLE),
        ("Final ranking", "50 reported", 11, RED),
    ]
    top, step = 81.0, 5.9
    cx = 55
    ax.text(1, 87, "Stage", fontsize=11, fontweight="bold", color=NAVY)
    ax.text(72, 87, "Design cap", fontsize=11, fontweight="bold", color=NAVY)
    for i, (name, cap, width, col) in enumerate(stages):
        y = top - i * step
        hgt = step - 1.1
        ax.text(1, y + hgt / 2, f"{i + 1:>2}.  {name}", ha="left", va="center", fontsize=10, color="black")
        bw = 0.26 * width
        ax.add_patch(Rectangle((cx - bw / 2, y), bw, hgt, fc=col, ec="white", lw=1.0, alpha=0.95))
        ax.text(72, y + hgt / 2, cap, ha="left", va="center", fontsize=10, color="black", linespacing=1.2)
    swatches(ax, [("ML triage", "teal"), ("physical + ADMET filters", "orange"), ("shortlist-only modeling", "purple")], 3.0)
    save(fig, "fig2_trackA_funnel")


def fig3_ml_design():
    fig, ax = canvas(10.5, 9.0)
    title(ax, "Track A machine-learning design",
          "Model selection never touches the locked test set. Scaffold-separated partitions replace temporal ones.")
    rows = [
        ("Labels", "Vina scores on 200,000 sampled molecules; the top 1% are the actives", "teal"),
        ("Scaffold split", "train 0.65, validation 0.15, calibration 0.10, locked test 0.10", "navy"),
        ("20 classifiers", "16 classical models + FT-Transformer x2 + RealMLP + Tabular ResNet", "navy"),
        ("Selection", "On validation only: best single model vs top-3 soft vote (average precision)", "navy"),
        ("Calibration", "Isotonic calibration + class-conditional split conformal (epsilon 0.15)", "orange"),
        ("Locked test", "ONE evaluation, after every choice is frozen", "purple"),
        ("Inference", "9 production-eligible models; 4,000 molecules per chunk;\n26 scoring workers + 2 readers", "teal"),
    ]
    bands(ax, rows, y_top=80, pitch=9.4, h=7.8, left_w=26)
    ax.text(1, 9.5, "Features: ECFP4 (radius 2, 2048 bits) + 10 descriptors. SMOTE is intentionally not used on fingerprint bits.\n"
                     "Actives are about 1%, so enrichment and recall are reported next to accuracy.\n"
                     "Slow benchmark models are row-capped and reported, but cannot become production models.",
            fontsize=9.5, color=GREY, va="center", linespacing=1.5)
    save(fig, "fig3_trackA_ml_design")


def fig4_resistance_panel():
    fig, ax = canvas(10.5, 9.8)
    title(ax, "Resistance-aware ensemble docking (shortlist stage)",
          "The same switch-II box is applied to every aligned receptor. Modeled variants are stress tests.")
    ax.text(1, 87, "Receptor panel (config/pdb_panel.tsv)", fontsize=11.5, fontweight="bold", color=NAVY)
    groups = [
        ("Primary G12D", ["7RPZ\nGDP", "7T47\nGTP-like"], "teal"),
        ("WT control", ["6MBU\nWT GDP"], "navy"),
        ("Experimental\nresistance", ["9U8J\nY96D", "9U95\nQ99L", "9U99\nE62K", "9U97\nG13D"], "orange"),
        ("Modeled variants\n(labeled MODEL_*)", ["V9E", "V9W", "V9Q", "G13P", "T58Y", "R68G", "Y96W", "H95Q", "H95L"], "purple"),
    ]
    y = 77.5
    for name, items, key in groups:
        fc, ec = COLORS[key]
        box(ax, 1, y, 24, 7.2, name, fc=ec, ec=ec, tc="white", bold=True, fs=10.5, align="left")
        n, gap = len(items), 0.9
        w = min(12.5, (99 - 26.5 - gap * (n - 1)) / n)
        for k, it in enumerate(items):
            box(ax, 26.5 + k * (w + gap), y, w, 7.2, it, fc=fc, ec=ec, fs=10 if n < 9 else 9.5)
        y -= 9.4
    ax.text(1, 43.5, "Score components (per compound)", fontsize=11.5, fontweight="bold", color=NAVY)
    comps = [
        ("G12D affinity   0.45", "Mean docking score over the G12D receptors (GDP and GTP-like)", "teal"),
        ("WT selectivity   0.25", "WT minus G12D docking score in kcal/mol; passes if >= 1.0", "navy"),
        ("Resistance robustness   0.30", "Worst resistance score minus G12D score; passes if <= 1.5", "purple"),
    ]
    y = 35.5
    for name, desc, key in comps:
        fc, ec = COLORS[key]
        box(ax, 1, y, 33, 7.2, name, fc=ec, ec=ec, tc="white", bold=True, fs=10.5, align="left")
        box(ax, 35, y, 64, 7.2, desc, fc=fc, ec=ec, fs=10, align="left")
        y -= 9.4
    ax.text(1, 7.0, "Each component is converted to a percentile rank, then combined with these weights.\n"
                    "Vina scores are not binding free energies. A WT margin is only meaningful against a matched WT\n"
                    "receptor (see docs/SCIENTIFIC_CAVEATS.md).",
            fontsize=9.5, color=GREY, va="center", linespacing=1.5)
    save(fig, "fig4_trackA_resistance_panel")


def fig5_trackB():
    fig, ax = canvas(10.5, 10.4)
    title(ax, "Track B: from a variant table to ranked ADAR-recruiting guides",
          "Thirteen stages run in order; each writes to a fresh directory. Preflight refuses stale or unverified inputs.")
    rows = [
        ("Variant table", "HGVS c. notation, one variant per row", "purple"),
        ("Reference check", "Reference base is checked against the validated CDS", "purple"),
        ("Correctability", "Is it G>A? Only then is A-to-I editing a direct reversal", "purple"),
        ("Guide layouts", "33 / 40 / 59 nt, with an orphan C opposite the target A", "purple"),
        ("Bystander map", "Every other A in the footprint, with codon consequence", "purple"),
        ("Chemistry transfer", "RESTORE 2.0 positional transfer, provenance-keyed", "teal"),
        ("Off-target screen", "Bowtie2 local binding of the guide against the sense transcriptome", "teal"),
        ("Gene-aware filter", "KRAS isoforms count as on target, not off target", "teal"),
        ("RNA structure", "ViennaRNA fold and cofold of guide and target", "teal"),
        ("ADAR1 geometry", "Boltz-2: ADAR1 716-1226 + RNA duplex (geometry only)", "orange"),
        ("Optional MD", "Needs validated modified-residue parameters", "orange"),
        ("Chemistry-aware ML", "Abstains without curated labels", "orange"),
        ("Research ranking", "Experiment-gated", "orange"),
    ]
    bands(ax, rows, y_top=84, pitch=5.35, h=4.4, left_w=30)
    swatches(ax, [("design (CPU)", "purple"), ("specificity and chemistry", "teal"), ("structure and triage (GPU)", "orange")], 9.0)
    ax.text(1, 3.6, "Hard boundaries: Boltz geometry is not catalytic evidence, ViennaRNA is not modification-aware,\nand no stage estimates editing efficiency.",
            fontsize=9.5, color=GREY, va="center", linespacing=1.5)
    save(fig, "fig5_trackB_flow")


def fig6_adar_concept():
    fig, ax = canvas(10.5, 10.4)
    title(ax, "Which KRAS variants can ADAR correct?",
          "ADAR deaminates adenosine to inosine, read as guanosine. Only G>A driver changes can be reversed directly.")
    ax.text(2, 86, "Variant", fontsize=10.5, fontweight="bold", color=NAVY)
    ax.text(16, 86, "cDNA change", fontsize=10.5, fontweight="bold", color=NAVY)
    ax.text(40, 86, "Direct A-to-I reversal", fontsize=10.5, fontweight="bold", color=NAVY)
    rows = [
        ("G12D", "c.35G>A", "Yes", TEAL), ("G12S", "c.34G>A", "Yes", TEAL), ("G13D", "c.38G>A", "Yes", TEAL),
        ("G12V", "c.35G>T", "No", RED), ("G12C", "c.34G>T", "No", RED), ("G12R", "c.34G>C", "No", RED),
        ("Q61H", "c.183A>C / c.183A>T", "No", RED),
    ]
    for r, (v, c, ok, col) in enumerate(rows):
        y = 80 - r * 5.2
        ax.add_patch(Rectangle((1, y - 2.3), 62, 4.8, fc=LIGHT if r % 2 == 0 else WHITE, ec="none"))
        ax.text(2, y, v, fontsize=10.5, va="center")
        ax.text(16, y, c, fontsize=10.5, va="center")
        ax.text(40, y, ok, fontsize=10.5, va="center", color=col, fontweight="bold")
    ax.text(66, 63, "The pipeline reports\nnon-correctable variants\nexplicitly instead of\nskipping them silently.", fontsize=10, color=GREY, va="center", linespacing=1.4)

    ax.text(1, 42, "Guide design principle (schematic)", fontsize=11.5, fontweight="bold", color=NAVY)
    target = list("GGUAGUUGGAGCUGUAGGCGUAGGCAAGAG")
    n = len(target)
    x0, dx = 9, 2.9
    ax.text(x0 - 0.8, 31, "mRNA 5'", fontsize=9, ha="right", va="center", color=GREY)
    ax.text(x0 + n * dx + 0.5, 31, "3'", fontsize=9, va="center", color=GREY)
    ax.text(x0 - 0.8, 21, "guide 3'", fontsize=9, ha="right", va="center", color=GREY)
    ax.text(x0 + n * dx + 0.5, 21, "5'", fontsize=9, va="center", color=GREY)
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    ti = 15
    for i, base in enumerate(target):
        xi = x0 + i * dx + dx / 2
        hot = i == ti
        ax.text(xi, 31, base, ha="center", va="center", fontsize=11, fontweight="bold" if hot else "normal",
                color=RED if hot else "black", family="DejaVu Sans Mono")
        pair = "C" if hot else comp[base]
        ax.text(xi, 21, pair, ha="center", va="center", fontsize=11, fontweight="bold" if hot else "normal",
                color=RED if hot else NAVY, family="DejaVu Sans Mono")
        ax.plot([xi, xi], [23.4, 28.6], color=RED if hot else GREY, lw=1.1, ls=":" if hot else "-")
    xt = x0 + ti * dx + dx / 2
    ax.annotate("target adenosine (edited A to I)", xy=(xt, 33.2), xytext=(xt + 8, 38.5),
                fontsize=10, color=RED, ha="left", arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.2))
    ax.annotate("orphan C opposite A (A-C mismatch recruits ADAR)", xy=(xt, 18.8), xytext=(xt - 40, 11.5),
                fontsize=10, color=RED, ha="left", arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.2))
    ax.text(1, 4, "Illustrative sequence only. Real guides come from the validated KRAS CDS. Layouts: 33, 40 and 59 nt\n"
                  "(RESTORE 2.0 constructs). Every other A in the footprint is a bystander, mapped for codon consequence.",
            fontsize=9.5, color=GREY, va="center", linespacing=1.5)
    save(fig, "fig6_adar_principle")


def fig7_trackC():
    fig, ax = canvas(10.5, 7.8)
    title(ax, "Track C: short peptide binders for 18 KRAS allele targets",
          "Backbone from 7RPZ chain A (residues 1-169). Mutant side chains are grafted: one substitution per target.")
    steps = ["1   PDB 7RPZ\n(G12D, MRTX1133 removed)", "2   Graft side chain at\ncodon 12, 13, 61, 117 or 146", "3   Condition spec:\nhotspots, H-bond\ndonor / acceptor",
             "4   RFD4-Proteina\n(RFdiffusion3 / Proteina)", "5   Complexes with an\n8-16 aa peptide", "6   Optional refolding metrics\n(AF3 / OpenFold3 oracle)"]
    w, gap = 30.0, 4.0
    for k, text in enumerate(steps):
        r, c = divmod(k, 3)
        x = 1 + c * (w + gap)
        y = 70 - r * 24
        box(ax, x, y, w, 16, text, fc=SOFT_ORANGE, ec=ORANGE, fs=10)
        if c:
            arrow(ax, x - gap + 0.2, y + 8, x - 0.2, y + 8)
    ax.text(1, 34, "Allele panel", fontsize=11.5, fontweight="bold", color=NAVY)
    ax.text(1, 27.5, "G12D   WT   G12C   G12V   G12R   G12A   G13D   A146T   G12S", fontsize=11, family="DejaVu Sans Mono")
    ax.text(1, 21.5, "G12F   G13C  G13R   Q61H   Q61L   Q61R   Q61K   A146V   K117N", fontsize=11, family="DejaVu Sans Mono")
    ax.text(1, 9, "Codons 12, 13 and 61 use the P-loop / switch-II pocket. K117 and A146 use the local nucleotide-site\n"
                  "residues. rfd4-proteina is an external dependency (cloned by scripts/bootstrap.sh c), not vendored here.",
            fontsize=9.5, color=GREY, va="center", linespacing=1.5)
    save(fig, "fig7_trackC_peptides")


def fig8_evidence_ladder():
    fig, ax = canvas(10.5, 8.4)
    title(ax, "What a computed result does and does not mean",
          "This repository stops at the first rung. Every higher rung needs a wet-lab experiment.")
    mol = ["Computed prioritization\n(this repository)", "Biochemical binding,\nGTP exchange inhibition", "WT KRAS / HRAS / NRAS\nselectivity", "pERK suppression",
           "G12D-selective cell killing", "Resistance-variant cells", "3D PDAC organoids + CAFs", "PK and tumor exposure", "In vivo efficacy and safety"]
    rna = ["Computed guide ranking\n(this repository)", "Synthesis + in vitro editing\nof KRAS RNA", "Editing in G12D cell lines\n(on-target, bystanders)", "Transcriptome-wide\noff-target editing",
           "Protein restoration,\npERK suppression", "Selective killing of\nmutant cells", "Organoids + CAFs", "Delivery (LNP / AAV)", "In vivo efficacy and safety"]
    n = len(mol)
    w, off = 40.0, 0.85
    ax.text(1, 87, "Small molecules and peptides", fontsize=11.5, fontweight="bold", color=TEAL)
    ax.text(52, 87, "ADAR RNA editing", fontsize=11.5, fontweight="bold", color=PURPLE)
    for i in range(n):
        y = 78 - i * 8.6
        first = i == 0
        box(ax, 1 + i * off, y, w, 7.2, mol[i], fc=SOFT_TEAL if first else LIGHT, ec=TEAL if first else GREY, fs=10,
            bold=first, ls="-" if first else "--")
        box(ax, 52 + i * off, y, w, 7.2, rna[i], fc=SOFT_PURPLE if first else LIGHT, ec=PURPLE if first else GREY, fs=10,
            bold=first, ls="-" if first else "--")
    save(fig, "fig8_evidence_ladder")


if __name__ == "__main__":
    for fn in (fig1_overview, fig2_funnel, fig3_ml_design, fig4_resistance_panel, fig5_trackB,
               fig6_adar_concept, fig7_trackC, fig8_evidence_ladder):
        fn()
