# Data access, licensing and what is not included

This repository contains **no datasets, results, models or third-party files**. Everything below is obtained by you, from the source, under the source's terms.

## Public inputs fetched by scripts

| Input | Source | Script | Terms |
|---|---|---|---|
| PDB entries 7RPZ, 7T47, 6MBU, 9U8J, 9U95, 9U99, 9U97 and the allele crystals in `config/kras_allele_panel.tsv` | RCSB PDB | `fetch_structures.py`, `build_resistance_panel.py`, `dock_kras_alleles.py` | PDB open data (CC0) |
| Reference ligand CCD 6IC | RCSB ligand service | `fetch_structures.py` | PDB open data |
| Ensembl release 113 cDNA + GTF (GRCh38) | ftp.ensembl.org | `tracks/b_adar_rna_editing/scripts/fetch_ensembl_reference.sh` | Ensembl terms of use |
| ADAR1 UniProt P55265 | rest.uniprot.org | `tracks/b_adar_rna_editing/scripts/setup_real_run.py` | CC BY 4.0 |
| ColabFold MSAs for KRAS alleles | public MMseqs2 server | `fetch_kras_msa.py`, Boltz `--use_msa_server` | sends the KRAS or ADAR1 sequence to a third party; avoid for confidential sequences and supply your own A3M instead |

## Inputs you must obtain yourself

| Input | Why not included | Where |
|---|---|---|
| **Enamine REAL** (13.6 B) | vendor license; redistribution not allowed | register at [enamine.net](https://enamine.net), download the `.cxsmiles.bz2` shards (about 330 GB), keep them compressed |
| Purchasable building-block inventory | third-party data | Syntheseus Retro* benchmark set, Figshare 10.6084/m9.figshare.25376728 |
| RESTORE 2.0 Supplementary Data and Source Data workbooks | publisher copyright | the Nature Communications article, DOI 10.1038/s41467-025-64434-7 |
| Boltz-2 weights | large, upstream license | downloaded automatically by Boltz into `~/.boltz` on first use |
| RetroChimera USPTO-50K checkpoint | large, upstream license | upstream RetroChimera instructions |
| RFD4-Proteina code and weights | upstream license | cloned and downloaded by `scripts/bootstrap.sh c` ([RosettaCommons/RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina)); override the URL with `RFD4_GIT_URL` |
| Approved-drug SMILES tables | your choice of source | `FDA_DRUG_PKL`, `FDA_DRUG_CSV` |
| Tumor/normal sequencing, RNA-seq BAMs, cell-line data | patient or experimental data | yours; never commit |
| GNINA binary | upstream release | github.com/gnina/gnina |

## Software licenses to check before you redistribute outputs (indicative only)

AutoDock Vina (Apache-2.0), Meeko (LGPL-2.1), RDKit (BSD-3), OpenMM (MIT/LGPL), AmberTools (mixed; free for academic use, see AmberTools license), ADMET-AI (MIT), Boltz (MIT), RetroChimera and Syntheseus (MIT), ViennaRNA (custom academic license), Bowtie2 (GPL-3), RFD4-Proteina (see upstream LICENSE), GNINA (GPL-2/Apache-2). Verify current terms upstream.
