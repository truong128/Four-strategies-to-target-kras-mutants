"""Adversarial regression tests for scoped outputs and chemical linkage semantics."""
import tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import yaml
from full.workflow import run
from full.common import GateError
from full.learning import transform

class AdversarialV14(unittest.TestCase):
    def test_internal_bond_cannot_be_blank(self):
        seq='ACGUA'
        chem=[{'position_1based':str(i+1),'base':b,'sugar':'RNA_2OH','link_to_next_3prime':'PO' if i<4 else 'NONE_3PRIME'} for i,b in enumerate(seq)]
        self.assertIsInstance(transform(seq,'UGCAU',2,chem),list)
        chem[1]['link_to_next_3prime']=''
        with self.assertRaisesRegex(GateError,'internal linkage'):
            transform(seq,'UGCAU',2,chem)
    def test_terminal_bond_not_counted(self):
        seq='ACGUA'
        chem=[{'position_1based':str(i+1),'base':b,'sugar':'RNA_2OH','link_to_next_3prime':'PO' if i<4 else 'NONE_3PRIME'} for i,b in enumerate(seq)]
        features=transform(seq,'UGCAU',2,chem)
        self.assertNotEqual(features[-1],1.0)
    def test_outside_output_index_is_blocked(self):
        with tempfile.TemporaryDirectory() as td:
            cfg=yaml.safe_load((Path(__file__).resolve().parents[1]/'config/full_project.yaml').read_text())
            cfg['out']=str(Path(td)/'run')
            cfg['offtarget']['transcriptome_fasta']=str(Path(td)/'transcriptome.fa')
            cfg['offtarget']['index_prefix']=str(Path(td)/'old_run'/'index')
            cfg['discovery']['variants_csv']=str(Path(__file__).resolve().parents[1]/'config/synthetic_for_integration.csv')
            cfg['discovery']['allow_synthetic_panel']=True
            cfg['discovery']['transcript']='SYNTHETIC_TEST'
            path=Path(td)/'config.yaml';path.write_text(yaml.safe_dump(cfg))
            # The earlier discovery stage is mocked to isolate index safety.
            with patch('full.workflow.call'):
                with self.assertRaisesRegex(GateError,'OFFTARGET_INDEX_OUTSIDE_RUN'):
                    run(SimpleNamespace(config=str(path),execute_heavy=False))
if __name__=='__main__':unittest.main()
