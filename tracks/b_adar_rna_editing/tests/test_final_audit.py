"""Regression tests for the highest-impact E2E scientific and execution gates."""
import csv,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from full.common import GateError,csvout,rows
from full.somatic import gatk,extract_panel
from full.learning import transform,train
from full.structure import predict
from full.dynamics import simulate
from full.workflow import run as full_run
from scripts.adar_pipeline import edit_variant,translate

class FinalAuditTests(unittest.TestCase):
    def test_kras_g12d_genetic_code(self):
        assert edit_variant('GGT',['c.2G>A'])=='GAT'
        assert translate('GGT')=='G' and translate('GAT')=='D'

    def test_gatk_current_normal_flag_and_orientation_plan(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t); ref=d/'ref.fa';ref.write_text('>chr12\nACGT\n');(d/'ref.fa.fai').write_text('')
            for name in ('t','n'):
                (d/f'{name}.bam').write_bytes(b'')
                (d/f'{name}.bam.bai').write_bytes(b'')
            gatk(SimpleNamespace(reference=str(ref),tumor_bam=str(d/'t.bam'),normal_bam=str(d/'n.bam'),tumor_sample='T',normal_sample='N',out=str(d/'run'),intervals=None,germline_resource=None,panel_of_normals=None,orientation_priors=True,contamination_table=None,execute=False,assembly='GRCh38'))
            arg=json.loads((d/'run'/'01_mutect2.log.command.json').read_text())['argv']
            self.assertNotIn('-tumor',arg)
            self.assertIn('--f1r2-tar-gz',arg)
            self.assertEqual(arg[arg.index('-normal')+1],'N')
            prior=json.loads((d/'run'/'orientation_model.log.command.json').read_text())['argv']
            self.assertIn('LearnReadOrientationModel',prior)

    def test_mutect_vcf_filter_unfiltered_dot_is_not_pass(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);p=d/'a.vcf';p.write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|SYMBOL|Feature|HGVSc">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n12\t1\t.\tC\tT\t.\t.\tCSQ=T|KRAS|ENST00000256078|ENST00000256078:c.35G>A\n12\t2\t.\tC\tT\t.\tPASS\tCSQ=T|KRAS|ENST00000256078|ENST00000256078:c.35G>A\n')
            extract_panel(SimpleNamespace(annotated_vcf=str(p),transcript='ENST00000256078',assembly='GRCh38',pass_only=True,out=d/'p.csv'))
            self.assertEqual(len(rows(d/'p.csv')),1)

    def test_missing_MSA_must_not_execute_boltz(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);inputs=d/'inputs';inputs.mkdir();(inputs/'test.yaml').write_text('version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: MAAA\n  - rna:\n      id: B\n      sequence: ACGU\n')
            with self.assertRaises(GateError):
                predict(SimpleNamespace(input_dir=str(inputs),out=str(d/'out'),devices=1,samples=1,use_msa_server=False,execute=False))

    def test_bad_modification_base_mapping_fails(self):
        s='AACCC';t='GGGAU'
        chem=[{'position_1based':i+1,'base':'C' if i==0 else s[i],'sugar':'2OME','link_to_next_3prime':'PO' if i<4 else 'NONE_3PRIME'} for i in range(5)]
        with self.assertRaises(GateError):transform(s,t,2,chem)
        chem[0]['base']=s[0]
        with self.assertRaises(GateError):transform(s,t,2,chem+[chem[0]])

    def test_external_model_target_group_leakage_in_training(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            guide='AACCC';target='GGGAU'
            fields=['construct_id','guide_rna_5to3','target_rna_5to3','guide_index_0based','target_group','assay_context','editing_fraction']
            csvout(d/'l.csv',[dict(construct_id='c1',guide_rna_5to3=guide,target_rna_5to3=target,guide_index_0based='2',target_group=g,assay_context='assay',editing_fraction='0.3') for g in ('target1','target2')],fields)
            csvout(d/'chem.csv',[{'construct_id':'c1','position_1based':i+1,'base':guide[i],'sugar':'2OME','link_to_next_3prime':'PO' if i<4 else 'NONE_3PRIME'} for i in range(5)],['construct_id','position_1based','base','sugar','link_to_next_3prime'])
            with self.assertRaises(GateError):train(SimpleNamespace(labels=str(d/'l.csv'),chemistry=str(d/'chem.csv'),out=str(d/'ml'),min_samples=1,min_groups=2,cv_folds=2))

    def test_dynamics_plan_protocol_without_runtime(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);(d/'fake.prmtop').write_text('PLACEHOLDER');(d/'fake.inpcrd').write_text('PLACEHOLDER')
            simulate(SimpleNamespace(topology=str(d/'fake.prmtop'),coordinates=str(d/'fake.inpcrd'),out=str(d/'md'),modified=False,chemistry_validation=None,replicates=3,minimize_steps=3000,equil_steps=100,production_steps=300,temperature=310.15,pressure=1.0,timestep_fs=2.0,execute=False,restraint_k=1000,platform='auto',collision_per_ps=1,cutoff_nm=.9,seed=1,report_interval=100,checkpoint_interval=100))
            self.assertEqual(json.loads((d/'md'/'run_protocol.json').read_text())['equilibration_phases'],'restrained_NPT_then_unrestrained_NPT')
            self.assertFalse((d/'md'/'rep_01'/'production.dcd').exists())

    def test_nonempty_output_rejected_to_prevent_stale_artifacts(self):
        with tempfile.TemporaryDirectory() as t:
            import yaml
            d=Path(t); output=d/'run';output.mkdir();(output/'ml').mkdir()
            (output/'ml'/'chemistry_model.joblib').write_bytes(b'STALE')
            base=Path(__file__).resolve().parents[1]
            cfg=yaml.safe_load((base/'config/full_project.yaml').read_text())
            cfg['out']=str(output)
            cp=d/'conf.yaml';cp.write_text(yaml.safe_dump(cfg))
            with self.assertRaisesRegex(GateError, 'OUTPUT_DIRECTORY_NOT_EMPTY'):
                full_run(SimpleNamespace(config=str(cp),execute_heavy=False))

    def test_default_configuration_never_uses_demo_patient_data(self):
        with tempfile.TemporaryDirectory() as t:
            import yaml
            base=Path(__file__).resolve().parents[1]
            cfg=yaml.safe_load((base/'config/full_project.yaml').read_text())
            cfg['out']=str(Path(t)/'run')
            c=Path(t)/'conf.yaml';c.write_text(yaml.safe_dump(cfg))
            full_run(SimpleNamespace(config=str(c),execute_heavy=False))
            self.assertEqual(json.loads((Path(t)/'run'/'RUN_STATUS.json').read_text())['status'],'INPUT_REQUIRED_NOT_EXECUTED')
            cfg['discovery']['variants_csv']='config/variants.csv'
            c.write_text(yaml.safe_dump(cfg))
            with self.assertRaises(GateError):full_run(SimpleNamespace(config=str(c),execute_heavy=False))

if __name__=='__main__':unittest.main()
