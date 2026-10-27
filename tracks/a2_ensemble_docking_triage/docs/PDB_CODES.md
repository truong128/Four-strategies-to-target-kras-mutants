# KRAS structural panel used by the pipeline

## Experimental structures used directly

| PDB | What it is | Role |
|---|---|---|
| **7RPZ** | KRAS G12D + GDP + MRTX1133 (6IC), 1.30 Å | Primary/reference GDP-state complex |
| **7T47** | KRAS G12D + GppCp + MRTX1133, 1.27 Å | Active-like nucleotide-state target |
| **6MBU** | WT KRAS + GDP + Mg, 1.45 Å | WT counter-screen |
| **9U8J** | KRAS G12D/Y96D + MRTX1133, 1.34 Å | Experimental resistance structure |
| **9U95** | KRAS G12D/Q99L + MRTX1133, 1.75 Å | Experimental resistance structure |
| **9U99** | KRAS G12D/E62K + MRTX1133, 2.50 Å | Experimental resistance structure |
| **9U97** | KRAS G12D/G13D + MRTX1133, 1.79 Å | Experimental secondary-mutant complex |

The runner downloads them from `https://files.rcsb.org/download/<PDB>.pdb` automatically.

## Resistance variants modeled from 7RPZ
The 2024 saturation/site-directed screen identified **V9E, V9W, V9Q, G13P, T58Y, R68G, Y96W, Q99L** as substantial MRTX1133-resistance variants. Q99L has an experimental complex above; the remaining listed variants are modeled from 7RPZ for *computational stress testing only*. H95Q/H95L are additionally modeled because H95 is a mechanistic KRAS-paralog selectivity handle and H95/Y96 substitutions have literature support as resistance liabilities.

**Do not describe a modeled mutant as an experimental structure.** The source type is explicitly recorded in `config/pdb_panel.tsv`.

## Useful additional experimental nucleotide-only structures
- **9U90** - KRAS G12D/E62K + GDP, 1.20 Å.
- **9U8W** - KRAS G12D/Q99L + GDP, 1.50 Å.
- **9NF2** - KRAS G12D + MRTX1133 + GMPPNP, 1.70 Å.

These are useful for follow-up mechanistic/state analysis but are not needed for the default panel.
