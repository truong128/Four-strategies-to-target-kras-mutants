#!/usr/bin/env bash
# Verified external command wrappers: run only when user supplies valid human assets.
set -euo pipefail
cmd="${1:-}"; shift || true
case "$cmd" in
 setup-conda)
   if ! command -v micromamba >/dev/null && ! command -v conda >/dev/null; then echo 'Install micromamba/conda first.' >&2; exit 3; fi
   backend="$(command -v micromamba || command -v conda)"
   "$backend" create -y -n kras_adar_tools -c conda-forge -c bioconda python=3.11 viennarna bowtie2 minimap2 samtools bcftools bedtools pysam biopython pandas openpyxl scikit-learn
   ;;
 vep)
   # Correct reference build and transcript cache must be supplied by user.
   if [[ $# -ne 3 ]]; then echo 'Usage: ./run_all.sh external vep input.vcf /path/to/vep_cache output.vcf' >&2; exit 2; fi
   input="$1"; cache="$2"; output="$3"
   [[ -f "$input" && -d "$cache" ]] || { echo 'Missing input VCF or verified VEP cache' >&2; exit 3; }
   command -v vep >/dev/null || { echo 'Ensembl VEP not installed.' >&2; exit 3; }
   vep --offline --cache --dir_cache "$cache" --assembly GRCh38 --species homo_sapiens --input_file "$input" --output_file "$output" --format vcf --vcf --force_overwrite --symbol --hgvs --canonical --mane --variant_class
   ;;
 samtools-check)
   [[ $# -eq 1 && -f "$1" ]] || { echo 'Usage: ./run_all.sh external samtools-check <input.bam>' >&2; exit 2; }
   samtools quickcheck -v "$1";samtools flagstat "$1"
   ;;
 adaredit-clone)
   [[ $# -eq 1 ]] || { echo 'Usage: ./run_all.sh external adaredit-clone <destination>' >&2; exit 2; }
   command -v git >/dev/null || exit 3
   [[ ! -e "$1" ]] || { echo 'Destination exists; refusing to overwrite' >&2; exit 3; }
   git clone https://github.com/Scientific-Computing-Lab/AdarEdit.git "$1"
   git -C "$1" rev-parse HEAD
   echo 'NOTE: model is out-of-domain for modified RESTORE guide efficiency; no automated KRAS predictions generated.'
   ;;
 *) echo 'Usage: ./run_all.sh external {setup-conda|vep|samtools-check|adaredit-clone}' >&2; exit 2 ;;
esac
