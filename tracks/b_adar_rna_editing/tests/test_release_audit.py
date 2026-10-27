"""Release regression checks for planned vs produced outputs and allele evidence."""
import tempfile, unittest, json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from full.common import GateError
from full.workflow import run

class ReleaseAudit(unittest.TestCase):
    def config(self, directory):
        import yaml
        root=Path(__file__).resolve().parents[1]
        cfg=yaml.safe_load((root/'config/full_project.yaml').read_text())
        cfg['out']=str(Path(directory)/'result')
        return cfg

    def launch(self, directory, cfg, execute=False):
        import yaml
        path=Path(directory)/'setup.yaml'
        path.write_text(yaml.safe_dump(cfg))
        run(SimpleNamespace(config=str(path),execute_heavy=execute))

    def test_planned_VEP_is_not_mistaken_for_an_output(self):
        with tempfile.TemporaryDirectory() as td:
            cfg=self.config(td)
            cache=Path(td)/'vep_cache';cache.mkdir()
            # A dummy source VCF is sufficient for a non-executing VEP plan.
            vcf=Path(td)/'input.vcf';vcf.write_text('##fileformat=VCFv4.2\n')
            cfg['somatic']['vep_cache']=str(cache)
            cfg['somatic']['annotate_vcf']=str(vcf)
            self.launch(td,cfg)
            result=json.loads((Path(td)/'result'/'RUN_STATUS.json').read_text())
            self.assertEqual(result['status'],'INPUT_REQUIRED_NOT_EXECUTED')
            self.assertIn('PLANNED_NOT_RUN',[s['status'] for s in result['stages']])
            self.assertFalse((Path(td)/'result'/'somatic'/'KRAS_somatic_panel.csv').exists())

    def test_missing_annotated_vcf_never_falls_back_to_demo(self):
        with tempfile.TemporaryDirectory() as td:
            cfg=self.config(td)
            cfg['somatic']['annotated_vcf']=str(Path(td)/'missing.vcf')
            with self.assertRaisesRegex(GateError,'annotated_vcf does not exist'):
                self.launch(td,cfg)

    def test_two_variant_sources_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            cfg=self.config(td)
            cfg['discovery']['variants_csv']=str(Path(td)/'variants.csv')
            cfg['somatic']['annotate_vcf']=str(Path(td)/'raw.vcf')
            with self.assertRaisesRegex(GateError,'Ambiguous competing mutation sources'):
                self.launch(td,cfg)

if __name__=='__main__':unittest.main()
