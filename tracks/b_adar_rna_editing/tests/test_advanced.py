import sys,tempfile,csv,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import advanced as adv
assert adv.revcomp('GAT')=='ATC'
with tempfile.TemporaryDirectory() as t:
 d=Path(t);input=d/'in';input.mkdir();out=d/'out'
 (input/'02_guides.csv').write_text('guide_id,variant_id,guide_rna_5to3,target_rna_5to3,guide_index_0based\ng1,G12D,AUCGAUCC,GGAUCGAA,3\n')
 (input/'03_bystander_sites.csv').write_text('guide_id,is_target\ng1,YES\ng1,NO\n')
 (d/'tx.fa').write_text('>tx1\nAAAAATCGATCCAAAA\n')
 class Args:pass
 a=Args();a.input=str(input);a.out=str(out);a.transcriptome=str(d/'tx.fa');a.max_mismatches=2;a.max_hits=100;a.md_manifest=None;a.measured_editing=None
 adv.run(a)
 assert len(adv.readcsv(out/'07_chemistry_positions.csv'))==16
 assert len(adv.readcsv(out/'12_candidate_evidence.csv'))==1
 assert adv.readcsv(out/'12_candidate_evidence.csv')[0]['therapeutic_priority']=='NOT_RANKABLE_WITHOUT_EXPERIMENTAL_EVIDENCE'
 assert json.loads((out/'STEPS_07_12_SUMMARY.json').read_text())['ML_status'].startswith('ABSTAIN')
print('PASS advanced: chemistry site manifest, transcriptome sequence scan, MD gate, ML abstention, honest evidence table.')
