import json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from full.evidence import prepared,topology
from full.parameterize import gen_tleap
from full.common import GateError
class EvidenceTests(unittest.TestCase):
 def test_no_automatic_validation_and_modification_hash_binding(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);pdb=p/'a.pdb';pdb.write_text('ATOM placeholder');lib=p/'a.lib';lib.write_text('LIB');fr=p/'a.frcmod';fr.write_text('PARAMS')
   evidence=p/'validated.json'
   prepared(SimpleNamespace(pdb=str(pdb),libs=[str(lib)],frcmods=[str(fr)],out=str(evidence)))
   v=json.loads(evidence.read_text());self.assertFalse(v['validated'])
   v.update({k:True for k in ('validated','stereochemistry_resolved','charge_and_torsion_checked','atom_mapping_checked')}); evidence.write_text(json.dumps(v))
   kw=dict(pdb=str(pdb),out=str(p/'out'),frcmod=[str(fr)],lib=[str(lib)],modified=True,parameter_validation_json=str(evidence),padding_A=12,salt_molar=.15,salt_pairs=None,execute=False)
   gen_tleap(SimpleNamespace(**kw))
   fr.write_text('CHANGED')
   with self.assertRaises(GateError):gen_tleap(SimpleNamespace(**kw))
 def test_md_checksum_scaffold_never_validated(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);f=p/'sys.prmtop';f.write_text('F');g=p/'sys.inpcrd';g.write_text('G')
   topology(SimpleNamespace(topology=str(f),coordinates=str(g),out=str(p/'e.json')))
   self.assertFalse(json.loads((p/'e.json').read_text())['validated'])
if __name__=='__main__':unittest.main()
