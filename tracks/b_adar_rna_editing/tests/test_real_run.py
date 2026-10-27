"""Regression tests for the real-run setup: CDS extraction, gene-aware off-target screening, triage and preflight."""
import csv, random, shutil, tempfile, unittest
from pathlib import Path
from types import SimpleNamespace

import yaml

from full.common import GateError
from full import offtarget, prioritize, workflow
from scripts import setup_real_run as setup

ROOT = Path(__file__).resolve().parents[1]


def rc(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def write_fasta_with_fai(path, records, width=60):
    with open(path, 'w') as fa, open(str(path) + '.fai', 'w') as fai:
        offset = 0
        for name, seq in records.items():
            header = f'>{name}\n'
            offset += len(header)
            fa.write(header)
            fai.write(f'{name}\t{len(seq)}\t{offset}\t{width}\t{width + 1}\n')
            for i in range(0, len(seq), width):
                line = seq[i:i + width] + '\n'
                fa.write(line)
                offset += len(line)


def csv_rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


class CdsExtraction(unittest.TestCase):
    def build_gene(self, strand):
        rng = random.Random(7)
        utr5 = ''.join(rng.choice('ACGT') for _ in range(40))
        cds = 'ATG' + 'GCT' * 20 + 'TAA'
        utr3 = ''.join(rng.choice('ACGT') for _ in range(30))
        mrna = utr5 + cds + utr3
        cut1, cut2 = 60, 95
        exons_tx = [(0, cut1), (cut1, cut2), (cut2, len(mrna))]
        base = 1000
        if strand == '+':
            exons = [(base + a, base + b - 1) for a, b in exons_tx]
            cds_tx = (40, 40 + len(cds) - 3)
            stop_tx = (40 + len(cds) - 3, 40 + len(cds))
            def to_genomic(a, b):
                out = []
                for ea, eb in exons_tx:
                    lo, hi = max(a, ea), min(b, eb)
                    if lo < hi:
                        out.append((base + lo, base + hi - 1))
                return out
        else:
            total = len(mrna)
            exons = [(base + total - b, base + total - a - 1) for a, b in exons_tx]
            cds_tx = (40, 40 + len(cds) - 3)
            stop_tx = (40 + len(cds) - 3, 40 + len(cds))
            def to_genomic(a, b):
                out = []
                for ea, eb in exons_tx:
                    lo, hi = max(a, ea), min(b, eb)
                    if lo < hi:
                        out.append((base + total - hi, base + total - lo - 1))
                return out
        feats = {'exon': exons, 'CDS': to_genomic(*cds_tx), 'stop_codon': to_genomic(*stop_tx)}
        return mrna, cds, feats

    def test_plus_and_minus_strand_cds_match_the_true_coding_sequence(self):
        for strand in ('+', '-'):
            mrna, cds, feats = self.build_gene(strand)
            offset, found = setup.cds_slice(feats, strand, mrna)
            self.assertEqual(offset, 40, strand)
            self.assertEqual(found, cds, strand)

    def test_fai_reader_returns_the_exact_record(self):
        with tempfile.TemporaryDirectory() as td:
            records = {'T1': 'ACGT' * 40, 'T2': 'TTGCA' * 33}
            fa = Path(td) / 't.fa'
            write_fasta_with_fai(fa, records)
            index = setup.read_fai(str(fa) + '.fai')
            for name, seq in records.items():
                self.assertEqual(setup.fetch_fasta_record(fa, index, name), seq)

    def test_kras_validation_rejects_a_non_kras_cds(self):
        with self.assertRaisesRegex(ValueError, 'failed validation'):
            setup.validate_kras_cds('ATG' + 'GCT' * 20 + 'TAA')

    def test_adar1_numbering_check_rejects_wrong_sequence(self):
        with self.assertRaises(ValueError):
            setup.check_adar1('M' * 1226)


class OffTargetAndTriage(unittest.TestCase):
    @unittest.skipUnless(shutil.which('bowtie2') and shutil.which('bowtie2-build'), 'bowtie2 not installed')
    def test_hits_are_labelled_by_gene_and_on_target_is_separated(self):
        rng = random.Random(11)
        kras = ''.join(rng.choice('ACGT') for _ in range(300))
        decoy = ''.join(rng.choice('ACGT') for _ in range(300))
        site = kras[100:140]
        guide = rc(site)
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            write_fasta_with_fai(td / 'tx.fa', {'ENST1': kras, 'ENST2': kras, 'ENST3': decoy})
            (td / 'map.tsv').write_text('transcript_id\tgene_id\tgene_name\nENST1\tG1\tKRAS\nENST2\tG1\tKRAS\nENST3\tG3\tOTHER\n')
            (td / 'guides.csv').write_text('guide_id,guide_rna_5to3\ng1,' + guide.replace('T', 'U') + '\n')
            offtarget.build(SimpleNamespace(transcriptome=str(td / 'tx.fa'), out=str(td / 'idx' / 'tx'), execute=True, threads=2))
            args = SimpleNamespace(guides=str(td / 'guides.csv'), index=str(td / 'idx' / 'tx'), out=str(td / 'screen'),
                                   min_coverage=0.7, max_alignments=50, threads=2, gene_map=str(td / 'map.tsv'),
                                   on_target_gene='KRAS', execute=True)
            offtarget.query(args)
            hits = csv_rows(td / 'screen' / 'offtarget_alignments.csv')
            self.assertEqual({h['transcript_id'] for h in hits}, {'ENST1', 'ENST2'})
            self.assertTrue(all(h['gene_name'] == 'KRAS' and h['is_on_target_gene'] == 'YES' for h in hits))
            summary = csv_rows(td / 'screen' / 'offtarget_gene_summary.csv')
            self.assertEqual(summary[0]['off_target_genes'], '0')

    def test_triage_excludes_on_target_hits_and_reports_zero_for_screened_guides(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            guide_cols = 'guide_id,variant_id,wt_guide_target_base,mutant_guide_target_base\n'
            (td / 'g.csv').write_text(guide_cols + 'a,V,G,A\nb,V,G,A\n')
            (td / 'b.csv').write_text('guide_id,is_target\na,NO\nb,NO\n')
            (td / 'o.csv').write_text('guide_id,gene_name,transcript_id,is_on_target_gene\n'
                                      'a,KRAS,T1,YES\na,NRAS,T2,NO\na,NRAS,T3,NO\n')
            args = SimpleNamespace(guides=str(td / 'g.csv'), bystanders=str(td / 'b.csv'), experimental=None,
                                   offtarget_hits=str(td / 'o.csv'), out=str(td / 'r'), min_replicates=3,
                                   min_edit=0.2, max_wt_bystander=0.05, max_offtarget=0.05, min_KRAS_GTP_reduction=0.2)
            prioritize.priority(args)
            rows = {r['guide_id']: r for r in csv_rows(td / 'r' / 'engineering_triage.csv')}
            self.assertEqual(rows['a']['complementary_transcript_hits_not_editing'], '2')
            self.assertEqual(rows['a']['complementary_off_target_genes'], '1')
            self.assertEqual(rows['a']['on_target_gene_transcript_hits'], '1')
            self.assertEqual(rows['b']['complementary_transcript_hits_not_editing'], '0')


class Preflight(unittest.TestCase):
    def run_preflight(self, td, mutate):
        cfg = yaml.safe_load((ROOT / 'config/real_run_pdac_kras.yaml').read_text())
        cfg['out'] = str(Path(td) / 'out')
        mutate(cfg)
        path = Path(td) / 'c.yaml'
        path.write_text(yaml.safe_dump(cfg))
        workflow.preflight(SimpleNamespace(config=str(path), min_free_gb=1.0))

    @unittest.skipUnless((ROOT / 'reference' / 'adar1').is_dir() and (ROOT / 'reference' / 'ensembl').is_dir(),
                         'needs downloaded references: run scripts/fetch_ensembl_reference.sh and ./run_all.sh setup-real')
    def test_bundled_real_config_passes(self):
        with tempfile.TemporaryDirectory() as td:
            self.run_preflight(td, lambda cfg: None)

    def test_wrong_reference_base_in_a_variant_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            bad = Path(td) / 'v.csv'
            bad.write_text('variant_id,hgvs_c\nBAD,c.35A>G\n')
            with self.assertRaises(GateError):
                self.run_preflight(td, lambda cfg: cfg['discovery'].update(variants_csv=str(bad)))

    def test_missing_transcriptome_is_caught(self):
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(GateError):
                self.run_preflight(td, lambda cfg: cfg['offtarget'].update(transcriptome_fasta=str(Path(td) / 'none.fa')))

    def test_nonempty_output_directory_is_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            def mutate(cfg):
                Path(cfg['out']).mkdir(parents=True)
                (Path(cfg['out']) / 'old.txt').write_text('x')
            with self.assertRaises(GateError):
                self.run_preflight(td, mutate)


if __name__ == '__main__':
    unittest.main()
