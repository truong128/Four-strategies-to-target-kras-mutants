#!/usr/bin/env bash
# Download Ensembl release 113 (GRCh38) cDNA + GTF and make the transcriptome FASTA that
# setup-real, the Bowtie2 off-target screen and the CDS cutter expect:
#   * sense-oriented transcripts,
#   * one header token per record (ENSTxxxxxxxxxxx, version and description removed),
#   * a samtools .fai index next to the FASTA.
# Needs: curl, gzip, awk, samtools. Output: reference/ensembl/ (git-ignored).
set -euo pipefail
cd "$(dirname "$0")/.."
REL="${ENSEMBL_RELEASE:-113}"
OUT="reference/ensembl"
BASE="https://ftp.ensembl.org/pub/release-${REL}"
mkdir -p "$OUT"
[[ -s "$OUT/Homo_sapiens.GRCh38.cdna.all.fa.gz" ]] || \
  curl -fsSL --retry 3 -o "$OUT/Homo_sapiens.GRCh38.cdna.all.fa.gz" "$BASE/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
[[ -s "$OUT/Homo_sapiens.GRCh38.${REL}.gtf.gz" ]] || \
  curl -fsSL --retry 3 -o "$OUT/Homo_sapiens.GRCh38.${REL}.gtf.gz" "$BASE/gtf/homo_sapiens/Homo_sapiens.GRCh38.${REL}.gtf.gz"
gzip -dc "$OUT/Homo_sapiens.GRCh38.${REL}.gtf.gz" > "$OUT/Homo_sapiens.GRCh38.${REL}.gtf"
if [[ "$REL" != "113" ]]; then ln -sf "Homo_sapiens.GRCh38.${REL}.gtf" "$OUT/Homo_sapiens.GRCh38.113.gtf"; fi
gzip -dc "$OUT/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
  | awk '/^>/ { split(substr($1,2), a, "."); print ">" a[1]; next } { print }' \
  > "$OUT/Homo_sapiens.GRCh38.cdna.sense.fa"
samtools faidx "$OUT/Homo_sapiens.GRCh38.cdna.sense.fa"
echo "Ensembl release ${REL} reference ready in $OUT"
