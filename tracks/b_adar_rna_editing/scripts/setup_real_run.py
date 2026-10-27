#!/usr/bin/env python3
"""Build the reference inputs for a real KRAS ADAR run from an Ensembl GTF plus a sense-oriented transcriptome FASTA.

Run scripts/fetch_ensembl_reference.sh first (or pass --gtf and --transcriptome for your own release).

* KRAS CDS: cut from a sense-oriented transcriptome FASTA using the CDS and stop-codon
  features of the matching GTF, so no Ensembl REST call is needed and the exact release
  used for off-target screening is also the release used for guide design.
* Transcript to gene map: lets the off-target screen separate KRAS on-target hits from
  paralogs (NRAS, HRAS) and unrelated genes.
* ADAR1 protein: UniProt P55265 full length (p150), p110 and a DRBM3 plus deaminase
  construct. Domain numbering is checked against catalytic residues before any file is used.
"""
import argparse, hashlib, json, re, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TX = 'ENST00000311936'
KRAS_PROTEIN_START = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'
CODONS = {a + b + c: 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i]
          for i, (a, b, c) in enumerate((a, b, c) for a in 'TCAG' for b in 'TCAG' for c in 'TCAG')}
ADAR1_CHECKS = {910: 'H', 912: 'E', 966: 'C', 1036: 'C'}
ADAR1_P110_START = 296
ADAR1_CONSTRUCT = (716, 1226)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def translate(seq):
    return ''.join(CODONS.get(seq[i:i + 3], 'X') for i in range(0, len(seq) - 2, 3))


def read_fai(fai_path):
    index = {}
    with open(fai_path) as f:
        for line in f:
            name, length, offset, bases, width = line.rstrip('\n').split('\t')[:5]
            index[name] = (int(length), int(offset), int(bases), int(width))
    return index


def fetch_fasta_record(fasta_path, index, name):
    if name not in index:
        raise KeyError(f'{name} not found in transcriptome index')
    length, offset, bases, width = index[name]
    nlines = (length + bases - 1) // bases
    with open(fasta_path, 'rb') as f:
        f.seek(offset)
        raw = f.read(length + nlines * (width - bases) + 2)
    return raw.decode().replace('\n', '').replace('\r', '')[:length].upper()


def parse_gtf_attributes(text):
    return dict(re.findall(r'(\w+) "([^"]*)"', text))


def transcript_features(gtf_path, transcript):
    needle = f'transcript_id "{transcript}"'
    feats = {'exon': [], 'CDS': [], 'stop_codon': []}
    strand = None
    version = ''
    gene = ''
    with open(gtf_path) as f:
        for line in f:
            if needle not in line:
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            kind = fields[2]
            strand = fields[6]
            if kind == 'transcript':
                attrs = parse_gtf_attributes(fields[8])
                version = attrs.get('transcript_version', '')
                gene = attrs.get('gene_name', '')
            if kind in feats:
                feats[kind].append((int(fields[3]), int(fields[4])))
    if not feats['exon'] or not feats['CDS'] or strand not in ('+', '-'):
        raise ValueError(f'{transcript}: exon or CDS features missing from {gtf_path}')
    return feats, strand, version, gene


def cds_slice(feats, strand, transcript_seq):
    minus = strand == '-'
    exons = sorted(feats['exon'], reverse=minus)
    cds = sorted(feats['CDS'], reverse=minus)
    stop = sorted(feats['stop_codon'], reverse=minus)
    first = cds[0][1] if minus else cds[0][0]
    offset = 0
    for start, end in exons:
        if start <= first <= end:
            offset += (end - first) if minus else (first - start)
            break
        offset += end - start + 1
    else:
        raise ValueError('first CDS base is not inside any exon')
    length = sum(e - s + 1 for s, e in cds) + sum(e - s + 1 for s, e in stop)
    return offset, transcript_seq[offset:offset + length]


def validate_kras_cds(cds):
    protein = translate(cds)
    problems = []
    if len(cds) % 3:
        problems.append('length is not a multiple of 3')
    if not cds.startswith('ATG'):
        problems.append('CDS does not start with ATG')
    if not protein.endswith('*') or '*' in protein[:-1]:
        problems.append('CDS must end with exactly one stop codon')
    if not protein.startswith(KRAS_PROTEIN_START[:60]):
        problems.append('N-terminal protein does not match human KRAS')
    if cds[33:36] not in ('GGT', 'GGC', 'GGA', 'GGG'):
        problems.append(f'codon 12 is {cds[33:36]}, expected Gly')
    if problems:
        raise ValueError('KRAS CDS failed validation: ' + '; '.join(problems))
    return protein


def build_kras_cds(gtf, transcriptome, transcript, out_dir):
    index = read_fai(str(transcriptome) + '.fai')
    feats, strand, version, gene = transcript_features(gtf, transcript)
    tx_seq = fetch_fasta_record(transcriptome, index, transcript)
    offset, cds = cds_slice(feats, strand, tx_seq)
    protein = validate_kras_cds(cds)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fasta = out / f'{transcript}_cds.fa'
    fasta.write_text(f'>{transcript}.{version} {gene} CDS cut from local GTF and transcriptome\n{cds}\n')
    meta = {'transcript': transcript, 'transcript_version': version, 'gene': gene, 'strand': strand,
            'cds_offset_in_transcript_0based': offset, 'cds_length_nt': len(cds),
            'protein_length_aa': len(protein) - 1, 'codon12': cds[33:36],
            'cds_sha256': sha256(fasta), 'gtf': str(gtf), 'transcriptome': str(transcriptome),
            'note': 'Codon numbering c.1 starts at the A of ATG. c.34 to c.38 cover G12 and G13.'}
    (out / f'{transcript}_cds.provenance.json').write_text(json.dumps(meta, indent=2) + '\n')
    return fasta, meta


def build_gene_map(gtf, out_path):
    count = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(gtf) as f, open(out_path, 'w') as o:
        o.write('transcript_id\tgene_id\tgene_name\n')
        for line in f:
            if '\ttranscript\t' not in line:
                continue
            attrs = parse_gtf_attributes(line.split('\t', 8)[8])
            tid = attrs.get('transcript_id')
            if not tid:
                continue
            gid = attrs.get('gene_id', '')
            o.write(f'{tid}\t{gid}\t{attrs.get("gene_name") or gid}\n')
            count += 1
    return count


def read_single_fasta(text):
    lines = [x.strip() for x in text.splitlines() if x.strip()]
    return lines[0], ''.join(x for x in lines[1:] if not x.startswith('>')).upper()


def check_adar1(seq):
    if len(seq) != 1226:
        raise ValueError(f'ADAR1 isoform 1 should be 1226 aa, got {len(seq)}')
    wrong = {pos: (want, seq[pos - 1]) for pos, want in ADAR1_CHECKS.items() if seq[pos - 1] != want}
    if wrong:
        raise ValueError(f'ADAR1 catalytic residue numbering mismatch (position: expected, found) {wrong}')
    if seq[ADAR1_P110_START - 1] != 'M':
        raise ValueError('ADAR1 residue 296 is not Met; p110 start does not match UniProt isoform 5')


def build_adar1(out_dir, accession='P55265', source_text=None):
    if source_text is None:
        url = f'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
        with urllib.request.urlopen(url, timeout=60) as response:
            source_text = response.read().decode()
    header, seq = read_single_fasta(source_text)
    check_adar1(seq)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    s, e = ADAR1_CONSTRUCT
    variants = {
        'ADAR1_p150_P55265-1_1-1226.fasta': (1, 1226),
        'ADAR1_p110_P55265-5_296-1226.fasta': (ADAR1_P110_START, 1226),
        f'ADAR1_DRBM3_deaminase_{s}-{e}.fasta': (s, e),
    }
    meta = {'uniprot': accession, 'source_header': header, 'catalytic_residues_checked': ADAR1_CHECKS,
            'numbering': 'UniProt P55265 isoform 1 (p150). p110 is isoform 5, residues 296 to 1226.',
            'domains_uniprot': {'DRBM3': '726-794', 'A_to_I_editase': '886-1221'}, 'files': {}}
    for name, (a, b) in variants.items():
        path = out / name
        path.write_text(f'>ADAR1_{accession}_{a}-{b}\n{seq[a - 1:b]}\n')
        meta['files'][name] = {'residues': f'{a}-{b}', 'length': b - a + 1, 'sha256': sha256(path)}
    (out / 'ADAR1_provenance.json').write_text(json.dumps(meta, indent=2) + '\n')
    return meta


def main():
    # Default inputs are produced by scripts/fetch_ensembl_reference.sh (Ensembl release 113, GRCh38).
    ref = ROOT / 'reference' / 'ensembl'
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--gtf', default=str(ref / 'Homo_sapiens.GRCh38.113.gtf'))
    p.add_argument('--transcriptome', default=str(ref / 'Homo_sapiens.GRCh38.cdna.sense.fa'))
    p.add_argument('--transcript', default=DEFAULT_TX)
    p.add_argument('--out-dir', default=str(ROOT / 'reference'))
    p.add_argument('--skip-adar', action='store_true')
    p.add_argument('--skip-gene-map', action='store_true')
    a = p.parse_args()
    out = Path(a.out_dir)
    fasta, meta = build_kras_cds(a.gtf, a.transcriptome, a.transcript, out / 'kras')
    print(f'KRAS CDS: {fasta} ({meta["cds_length_nt"]} nt, {meta["protein_length_aa"]} aa, version {meta["transcript_version"]})')
    if not a.skip_gene_map:
        n = build_gene_map(a.gtf, out / 'annotation' / 'transcript_gene_map.tsv')
        print(f'gene map: {n} transcripts')
    if not a.skip_adar:
        info = build_adar1(out / 'adar1')
        print('ADAR1 files:', ', '.join(info['files']))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('FATAL:', exc, file=sys.stderr)
        sys.exit(2)
