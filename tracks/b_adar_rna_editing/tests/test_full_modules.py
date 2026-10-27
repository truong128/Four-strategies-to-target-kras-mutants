import csv,importlib,tempfile,unittest,json
from pathlib import Path
from types import SimpleNamespace
from full.learning import transform
from full.chemistry import normalize
from full.prioritize import priority
from full.somatic import extract_panel
from full.common import csvout,rows,GateError

class FullPipelineTests(unittest.TestCase):
 def test_chemistry_position_features(self):
  s='U'*10+'C'+'G'*10; t='C'*10+'A'+'G'*10
  chem=[{'position_1based':str(i+1),'sugar':'2OME' if i%2 else '2F','link_to_next_3prime':('NONE_3PRIME' if i==len(s)-1 else ('PS' if i%3 else 'PO'))} for i in range(len(s))]
  f=transform(s,t,10,chem);self.assertGreater(len(f),80)
  with self.assertRaises(GateError):transform(s,t,10,chem[:-1])
 def test_vcf_panel_extraction(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'vep.vcf';p.write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|SYMBOL|Feature|HGVSc|HGVSp">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n12\t25398284\t.\tC\tT\t.\tPASS\tCSQ=T|KRAS|ENST00000256078|ENST00000256078:c.35G>A|p.Gly12Asp\n')
   out=Path(tmp)/'panel.csv';extract_panel(SimpleNamespace(annotated_vcf=str(p),transcript='ENST00000256078',assembly='GRCh38',pass_only=True,out=out))
   self.assertEqual(rows(out)[0]['hgvs_c'],'c.35G>A')
 def test_rank_abstains_without_measurement(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);guides=[{'guide_id':'G12D_1','variant_id':'G12D','wt_guide_target_base':'G','mutant_guide_target_base':'A'}]
   csvout(p/'guides.csv',guides,list(guides[0]))
   csvout(p/'bystanders.csv',[{'guide_id':'G12D_1','is_target':'NO'}],['guide_id','is_target'])
   priority(SimpleNamespace(guides=str(p/'guides.csv'),bystanders=str(p/'bystanders.csv'),offtarget_hits=None,experimental=None,out=str(p/'rank'),min_replicates=3,min_edit=.2,max_wt_bystander=.05,max_offtarget=.05,min_KRAS_GTP_reduction=.2))
   self.assertEqual(json.loads((p/'rank'/'ranking_status.json').read_text())['status'],'NO_THERAPEUTIC_RANKING_ABSTAIN')
   self.assertIn('WT_G_C_PERFECT',rows(p/'rank'/'engineering_triage.csv')[0]['WT_caveat'])
 def test_published_positions_transfer(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp);g={'guide_id':'g','variant_id':'G12D','guide_rna_5to3':'AACCC','target_rna_5to3':'GGGUU','guide_index_0based':'2','guide_length':'5'}
   csvout(d/'guides.csv',[g],list(g));chem=[{'transfer_id':'t1','guide_id':'g','position_1based':str(i+1),'kras_base':g['guide_rna_5to3'][i],'sugar':'RNA_2OH','link_to_next_3prime':'PO' if i<4 else 'NONE_3PRIME','template_source_key':'paper::ref'} for i in range(5)]
   csvout(d/'transfer.csv',chem,list(chem[0]));normalize(SimpleNamespace(guides=str(d/'guides.csv'),transferred=str(d/'transfer.csv'),out=str(d/'out')))
   self.assertEqual(len(rows(d/'out'/'constructs.csv')),1)

if __name__=='__main__':unittest.main()
