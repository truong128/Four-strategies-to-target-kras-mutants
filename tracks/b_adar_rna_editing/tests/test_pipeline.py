import sys,tempfile,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import adar_pipeline as a
assert a.rc('GAU'.replace('U','T'))=='ATC'
assert a.GENETIC['GGT']=='G' and a.GENETIC['GAT']=='D'
# Synthetic valid CDS for deterministic offline testing; never claim it is true KRAS sequence.
s='ATG'+'GCT'*10+'GGT'+'GGT'+'GCT'*99
assert s[33:36]=='GGT'
assert a.edit_variant(s,['c.35G>A'])[33:36]=='GAT'
assert a.edit_variant(s,['c.34G>A'])[33:36]=='AGT'
assert a.edit_variant(s,['c.38G>A'])[36:39]=='GAT'
with tempfile.TemporaryDirectory() as t:
 d=Path(t);(d/'seq.fa').write_text('>SYNTHETIC_CDS\n'+s+'\n')
 (d/'variants.csv').write_text('variant_id,hgvs_c\nG12D,c.35G>A\nG12V,c.35G>T\n')
 class A:pass
 z=A();z.out=str(d/'out');z.cds_fasta=str(d/'seq.fa');z.transcript='SYNTHETIC';z.variants=str(d/'variants.csv');z.lengths=[33];z.fold=False;z.transcriptome=None
 a.run(z)
 with (d/'out/01_variants.csv').open() as f: v=list(csv.DictReader(f))
 assert v[0]['directly_correctable']=='YES' and v[1]['directly_correctable']=='NO'
 with (d/'out/02_guides.csv').open() as f:g=list(csv.DictReader(f))
 assert len(g)>0 and all(x['guide_rna_5to3'][int(x['guide_index_0based'])]=='C' for x in g)
 assert (d/'out/08_md_input_manifest.csv').exists()
print('PASS: reference checks, variant edits, editability and guide mismatch on synthetic CDS.')
