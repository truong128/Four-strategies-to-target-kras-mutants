"""Full KRAS-ADAR executable orchestration. Input-heavy steps use real executors.
Missing scientific inputs yield explicit SKIPPED_INPUT_REQUIRED status, never a fabricated result.
"""
import argparse,datetime,json,os,shutil,subprocess,sys,traceback
from pathlib import Path
from full.common import GateError,jsonout,path_required
ROOT=Path(__file__).resolve().parents[1]

def resolve(value):
    if not value:return ''
    p=Path(value).expanduser()
    return str(p.resolve() if p.is_absolute() else (ROOT/p).resolve())

def call(name,args,report,needed=True):
    args=[str(a) for a in args]
    cmd=[sys.executable,*args]
    record={'stage':name,'command':cmd,'status':'RUNNING'}
    report.append(record)
    try:
        cp=subprocess.run(cmd,cwd=str(ROOT),check=False)
        if cp.returncode:raise GateError(f'Exit code {cp.returncode} from {name}')
        record['status']='COMPLETED'
    except Exception as e:
        record.update(status='FAILED',reason=str(e))
        raise

def skipped(name,report,why):
    report.append({'stage':name,'status':'SKIPPED_INPUT_REQUIRED','reason':why})

def run(args):
    try:import yaml
    except ImportError as e:raise GateError('pip install PyYAML') from e
    config=path_required(args.config,'--config')
    cfg=yaml.safe_load(config.read_text())
    if not isinstance(cfg,dict):raise GateError('Configuration must be YAML mapping')
    out=Path(resolve(cfg.get('out','results/full_e2e')))
    # A new research run must never silently consume outputs from an earlier run.
    if out.exists() and any(out.iterdir()):
        raise GateError('OUTPUT_DIRECTORY_NOT_EMPTY: choose a fresh config.out directory for each run. Existing models, VCFs, and off-target hits must never be reused implicitly.')
    out.mkdir(parents=True,exist_ok=True)
    report=[]
    try:
        d=cfg.get('discovery',{})
        sm=cfg.get('somatic',{})
        if d.get('variants_csv'):
            variant_path=resolve(d['variants_csv'])
            example_panel=(Path(variant_path).name in ('variants.csv','synthetic_for_integration.csv')
                           and Path(variant_path).parent == ROOT/'config')
            example_tx=str(d.get('transcript','')).upper().startswith('SYNTHETIC')
            if (example_panel or example_tx) and not d.get('allow_synthetic_panel',False):
                raise GateError('DEMO INPUT BLOCKED: set discovery.variants_csv to your verified tumor variants. For synthetic testing only, set allow_synthetic_panel: true.')
        if d.get('variants_csv') and (sm.get('annotated_vcf') or sm.get('annotate_vcf') or (sm.get('vep_cache') and (sm.get('tumor_bam') or sm.get('tumor_fastq',{}).get('read1')))):
            raise GateError('Ambiguous competing mutation sources: choose discovery.variants_csv OR somatic.annotated_vcf, not both.')
        # Optional FASTQ->BAM and matched tumor/normal calling, fully executable on user hardware.
        for sample in ('tumor','normal'):
            fastq=sm.get(sample+'_fastq',{})
            if fastq.get('read1') and fastq.get('read2'):
                cmd=['-m','full.alignment','--reference',resolve(sm['reference_fasta']),'--read1',resolve(fastq['read1']),'--read2',resolve(fastq['read2']),'--sample',fastq['sample_id'],'--out',str(out/'alignment'/sample)]
                if sm.get('known_sites'):cmd+=['--known-sites',resolve(sm['known_sites'])]
                if args.execute_heavy:cmd.append('--execute')
                call('00_align_'+sample,cmd,report)
                if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
        if args.execute_heavy:
            for sample in ('tumor','normal'):
                if not sm.get(sample+'_bam') and sm.get(sample+'_fastq',{}).get('read1'):
                    f=out/'alignment'/sample
                    sm[sample+'_bam']=str((f/'recalibrated.bam') if sm.get('known_sites') else (f/'duplicate_marked.bam'))
        if sm.get('tumor_bam') and sm.get('normal_bam'):
            cmd=['-m','full.somatic','mutect2','--reference',resolve(sm['reference_fasta']),'--tumor-bam',resolve(sm['tumor_bam']),'--normal-bam',resolve(sm['normal_bam']),'--tumor-sample',sm['tumor_sample'],'--normal-sample',sm['normal_sample'],'--out',str(out/'somatic')]
            if sm.get('orientation_priors',False):cmd+=['--orientation-priors']
            if sm.get('contamination_table'):cmd+=['--contamination-table',resolve(sm['contamination_table'])]
            if sm.get('intervals'):cmd+=['--intervals',resolve(sm['intervals'])]
            if sm.get('germline_resource'):cmd+=['--germline-resource',resolve(sm['germline_resource'])]
            if sm.get('panel_of_normals'):cmd+=['--panel-of-normals',resolve(sm['panel_of_normals'])]
            if args.execute_heavy:cmd.append('--execute')
            call('00_mutect2_somatic_calling',cmd,report)
            if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
        if sm.get('vep_cache') and not sm.get('annotated_vcf'):
            vcf=sm.get('annotate_vcf') or str(out/'somatic'/'tumor_normal.filtered.vcf.gz')
            if Path(resolve(vcf)).is_file() or (args.execute_heavy and sm.get('tumor_bam') and sm.get('normal_bam')):
                cmd=['-m','full.somatic','vep','--vcf',resolve(vcf),'--cache',resolve(sm['vep_cache']),'--assembly',sm.get('assembly','GRCh38'),'--out',str(out/'somatic'/'vep_annotated.vcf')]
                if args.execute_heavy:cmd.append('--execute')
                call('00_variant_effect_prediction',cmd,report)
                if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
                elif not sm.get('annotated_vcf'):sm['annotated_vcf']=str(out/'somatic'/'vep_annotated.vcf')
        panel_vcf=sm.get('annotated_vcf')
        if panel_vcf and not Path(resolve(panel_vcf)).is_file():
            raise GateError('Configured annotated_vcf does not exist: '+str(panel_vcf))
        if not d.get('variants_csv') and panel_vcf:
            outfile=out/'somatic'/'KRAS_somatic_panel.csv'
            call('00_extract_KRAS_HGVS_panel',['-m','full.somatic','panel','--annotated-vcf',resolve(panel_vcf),'--out',str(outfile),'--transcript',d.get('transcript','ENST00000256078'),'--pass-only'],report)
            d['variants_csv']=str(outfile)
        if not d.get('variants_csv'):
            if sm.get('vep_cache') and not args.execute_heavy:
                skipped('00_extract_KRAS_HGVS_panel',report,'VEP output is planned but not executed; cannot extract HGVS from a nonexistent VCF')
            warning='Verified tumor variant panel or completed VEP annotation is required; example mutation panel is not automatically used.'
            skipped('01_variant_annotation_and_guide_design',report,warning)
            jsonout(out/'RUN_STATUS.json',{'status':'INPUT_REQUIRED_NOT_EXECUTED','time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':report,'required_next':'Set discovery.variants_csv to verified variants or configure and run tumor-normal somatic/vep workflow; make sure genome build, transcript, and allele phase are reviewed.','important':'No candidates were produced and no clinical claim is made.'})
            (out/'INPUT_REQUIREMENTS.md').write_text('# Input required before KRAS ADAR discovery\n\n- Verified somatic KRAS variants in HGVS c. notation and a version-matched KRAS CDS.\n- Or complete matched tumor/normal reads, GRCh38, VEP cache and associated quality controls.\n- Do not interpret config/variants.csv as patient observations.\n- Downstream structures, chemistry parameters, edit-efficiency labels, and matched controls must be independently validated.\n')
            if args.execute_heavy:raise GateError(warning)
            print('INPUT_REQUIRED_NOT_EXECUTED: see '+str(out/'RUN_STATUS.json'))
            return
        command=['scripts/adar_pipeline.py','--variants',resolve(d['variants_csv']),'--out',str(out/'discovery'),'--transcript',d.get('transcript','ENST00000256078'),'--lengths',*[str(x) for x in d.get('lengths',[33,40,59])]]
        if d.get('cds_fasta'):command+=['--cds-fasta',resolve(d['cds_fasta'])]
        call('01_variant_annotation_and_guide_design',command,report)
        call('02_bystanders_structure_approximation',['scripts/advanced.py','--input',str(out/'discovery'),'--out',str(out/'discovery')],report)
        r=cfg.get('restore2',{})
        reference=Path(resolve(r.get('publisher_dir','reference/publisher')))
        if r.get('supplement_xlsx') and r.get('source_xlsx'):
            call('03_import_publisher_workbooks',['scripts/ingest_publisher_workbooks.py','--supplement',resolve(r['supplement_xlsx']),'--source-data',resolve(r['source_xlsx']),'--out',str(reference)],report)
        if (reference/'publisher_constructs.csv').is_file() and (reference/'publisher_residues.csv').is_file():
            dest=out/'chemistry_transfer';dest.mkdir(exist_ok=True)
            call('04_source_verified_chemistry_transfer',['scripts/transfer_publisher_chemistry.py','--reference-dir',str(reference),'--guides',str(out/'discovery'/'02_guides.csv'),'--out',str(dest)],report)
            call('05_chemistry_normalization',['-m','full.chemistry','--guides',str(out/'discovery'/'02_guides.csv'),'--transferred',str(dest/'transferred_per_residue.csv'),'--out',str(out/'chemistry')],report)
        else:skipped('03_05_published_chemistry',report,'publisher workbook reference library missing')
        of=cfg.get('offtarget',{})
        if of.get('transcriptome_fasta'):
            index=resolve(of.get('index_prefix') or str(out/'index'/'transcriptome'))
            if not Path(index).is_relative_to(out.resolve()):
                raise GateError('OFFTARGET_INDEX_OUTSIDE_RUN: index must live inside the fresh run directory to prevent stale cross-run reuse')
            call('06_transcriptome_build_index',['-m','full.offtarget','index','--transcriptome',resolve(of['transcriptome_fasta']),'--out',index,'--threads',str(of.get('threads',8))]+(['--execute'] if args.execute_heavy else []),report)
            scan=['-m','full.offtarget','screen','--guides',str(out/'discovery'/'02_guides.csv'),'--index',index,'--out',str(out/'offtarget'),'--threads',str(of.get('threads',8))]
            if of.get('gene_map'):scan+=['--gene-map',resolve(of['gene_map']),'--on-target-gene',str(of.get('on_target_gene','KRAS'))]
            call('07_transcriptome_local_scan',scan+(['--execute'] if args.execute_heavy else []),report)
            if not args.execute_heavy:
                report[-1]['status']='PLANNED_NOT_RUN'
                report[-2]['status']='PLANNED_NOT_RUN'
        else:skipped('06_07_transcriptome_screen',report,'offtarget.transcriptome_fasta not set')
        mod=cfg.get('structure',{})
        if mod.get('adar_fasta'):
            prep=['-m','full.structure','prepare','--guides',str(out/'discovery'/'02_guides.csv'),'--adar-fasta',resolve(mod['adar_fasta']),'--msa',mod.get('msa','auto'),'--out',str(out/'boltz_input')]
            if mod.get('allow_single_sequence'):prep.append('--allow-single-sequence')
            if mod.get('variants'):prep+=['--variants',*[str(x) for x in mod['variants']]]
            if mod.get('max_guides'):prep+=['--max-guides',str(mod['max_guides'])]
            call('08_generate_ADAR_RNA_complex_inputs',prep,report)
            pred=['-m','full.structure','predict','--input-dir',str(out/'boltz_input'),'--out',str(out/'boltz_models'),'--devices',str(mod.get('devices',1)),'--samples',str(mod.get('samples',5))]
            if mod.get('use_msa_server'):pred.append('--use-msa-server')
            if mod.get('no_kernels',True):pred.append('--no-kernels')
            if args.execute_heavy:pred.append('--execute')
            call('09_predict_protein_RNA_canonical_complexes',pred,report)
            if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
            else:
                call('09b_ADAR_RNA_geometry',['-m','full.structure','geometry','--index',str(out/'boltz_input'/'structure_index.csv'),'--models',str(out/'boltz_models'),'--out',str(out/'boltz_models'/'geometry.csv')],report)
        else:skipped('08_09_ADAR_RNA_complex_modeling',report,'structure.adar_fasta required to generate protein-RNA Boltz-2 jobs')
        ff=cfg.get('forcefield',{})
        if ff.get('build_specs_csv'):
            import csv
            with open(resolve(ff['build_specs_csv'])) as file:specs=list(csv.DictReader(file))
            for spec in specs:
                ident=spec['construct_id'];folder=out/'topologies'/ident
                cmd=['-m','full.parameterize','--pdb',resolve(spec['pdb']),'--out',str(folder)]
                if spec.get('modified','').lower()=='yes':cmd+=['--modified','--parameter-validation-json',resolve(spec['validation_json'])]
                if spec.get('libs'):cmd+=['--lib',*[resolve(x) for x in spec['libs'].split(';') if x]]
                if spec.get('frcmods'):cmd+=['--frcmod',*[resolve(x) for x in spec['frcmods'].split(';') if x]]
                cmd+=['--padding-A',str(ff.get('padding_A',12)),'--salt-molar',str(ff.get('salt_molar',0.15))]
                if spec.get('salt_pairs'):cmd+=['--salt-pairs',str(spec['salt_pairs'])]
                if args.execute_heavy:cmd.append('--execute')
                call('10_build_validated_Amber_'+ident,cmd,report)
                if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
        else:skipped('10_modified_forcefield_build',report,'forcefield.build_specs_csv needed: validated structure, residue libs/frcmods and chirality evidence')
        md=cfg.get('md',{})
        jobs=[]
        if md.get('jobs_csv'):
            import csv
            with open(resolve(md['jobs_csv'])) as file:jobs=list(csv.DictReader(file))
        elif args.execute_heavy and ff.get('build_specs_csv') and md.get('auto_jobs_from_forcefield',True):
            for spec in specs:
                folder=out/'topologies'/spec['construct_id']
                jobs.append({'construct_id':spec['construct_id'],'topology_prmtop':str(folder/'system.prmtop'),'coordinates_inpcrd':str(folder/'system.inpcrd'),'modified':spec.get('modified',''),'validation_json':spec.get('md_validation_json','')})
        if jobs:
            for job in jobs:
                ident=job['construct_id'];cmd=['-m','full.dynamics','simulate','--topology',resolve(job['topology_prmtop']),'--coordinates',resolve(job['coordinates_inpcrd']),'--out',str(out/'md'/ident),'--replicates',str(md.get('replicates',3)),'--production-steps',str(md.get('production_steps',5000000)),'--equil-steps',str(md.get('equil_steps',250000)),'--platform',md.get('platform','auto')]
                if job.get('modified','').lower()=='yes':
                    if not job.get('validation_json'):raise GateError('Missing post-parameterization MD validation JSON for '+ident+'. It must contain topology and coordinate SHA256 hashes.')
                    cmd+=['--modified','--chemistry-validation',resolve(job['validation_json'])]
                if args.execute_heavy:cmd.append('--execute')
                call('11_physiological_MD_'+ident,cmd,report)
                if not args.execute_heavy:report[-1]['status']='PLANNED_NOT_RUN'
                if args.execute_heavy:
                    traces=[str(out/'md'/ident/f'rep_{i:02d}'/'production.dcd') for i in range(1,int(md.get('replicates',3))+1)]
                    # Output PDB is produced after equilibration, retaining atom order/topology.
                    if all(Path(t).is_file() for t in traces):
                        call('11b_MD_conformational_analysis_'+ident,['-m','full.dynamics','analyze','--trajectories',*traces,'--topology-pdb',str(out/'md'/ident/'rep_01'/'equilibrated.pdb'),'--out',str(out/'md'/ident/'analysis')],report)
                    else:skipped('11b_MD_conformational_analysis_'+ident,report,'Production trajectories missing')
        else:skipped('11_modified_MD',report,'md.jobs_csv with verified solvated topology/coordinates required')
        ml=cfg.get('machine_learning',{})
        if ml.get('curated_labels_csv') and ml.get('curated_chemistry_positions_csv'):
            call('12_cross_target_ML_training',['-m','full.learning','train','--labels',resolve(ml['curated_labels_csv']),'--chemistry',resolve(ml['curated_chemistry_positions_csv']),'--out',str(out/'ml')],report)
            trained=out/'ml'/'chemistry_model.joblib'
            report_path=out/'ml'/'ml_report.json'
            model_valid=(trained.is_file() and report_path.is_file()
                         and json.loads(report_path.read_text()).get('status')=='EXPLORATORY_CROSS_TARGET_MODEL_NOT_EXTERNALLY_VALIDATED_ON_KRAS')
            if model_valid and ml.get('external_labels_csv') and ml.get('external_chemistry_positions_csv'):
                call('12b_external_independent_target_validation',['-m','full.learning','validate-external','--model',str(trained),'--provenance',str(out/'ml'/'ml_provenance.json'),'--labels',resolve(ml['external_labels_csv']),'--chemistry',resolve(ml['external_chemistry_positions_csv']),'--out',str(out/'ml'/'external')],report)
            elif model_valid:skipped('12b_external_independent_target_validation',report,'External labeled holdout with nonoverlapping target groups is missing')
            if model_valid and (out/'chemistry'/'constructs.csv').is_file():
                call('13_predict_KRAS_exploratory_OOD',['-m','full.learning','infer','--model',str(trained),'--chemistry',str(out/'chemistry'/'positions.csv'),'--candidates',str(out/'chemistry'/'constructs.csv'),'--out',str(out/'ml'/'KRAS_out_of_domain_predictions.csv')],report)
            else:skipped('13_KRAS_ML_predictions',report,'No model beating baseline or no verified KRAS chemistry designs')
        else:skipped('12_13_chemistry_ML',report,'Require curated editing yields linked to exact constructs + position chemistry')
        exp=cfg.get('experimental_rnaseq',{})
        if all(exp.get(k) for k in ('treated_bam','control_bam','reference_fasta','sites_csv')):
            call('13b_RNA_edit_site_counts_NOT_allele_specific_rate',['-m','full.experimental','--treated-bam',resolve(exp['treated_bam']),'--control-bam',resolve(exp['control_bam']),'--reference-fasta',resolve(exp['reference_fasta']),'--sites',resolve(exp['sites_csv']),'--out',str(out/'experimental'/'editing_site_raw_counts.csv')],report)
        else:skipped('13b_experimental_RNAseq_evidence',report,'Requires treated/control RNA BAM, reference FASTA and nominated edit sites; does not measure mutant-specific rate by itself')
        last=['-m','full.prioritize','--guides',str(out/'discovery'/'02_guides.csv'),'--bystanders',str(out/'discovery'/'03_bystander_sites.csv'),'--out',str(out/'ranking')]
        if cfg.get('ranking',{}).get('experimental_csv'):last+=['--experimental',resolve(cfg['ranking']['experimental_csv'])]
        oth=out/'offtarget'/'offtarget_alignments.csv'
        if of.get('transcriptome_fasta') and args.execute_heavy and oth.exists():last+=['--offtarget-hits',str(oth)]
        call('14_engineering_triage_and_conditional_ranking',last,report)
    except Exception:
        jsonout(out/'RUN_STATUS.json',{'status':'FAILED','time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':report,'last_exception':traceback.format_exc()[-2500:]})
        raise
    jsonout(out/'RUN_STATUS.json',{'status':'FINISHED_WITH_CONDITIONAL_MODULES','time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stages':report,'warning':'FINISHED means scripts completed, NOT validated ADAR editing or therapy; skipped stages and experimental validation remain mandatory.'})
    print(f'Completed orchestration, see {out/"RUN_STATUS.json"}')

def doctor(args):
    env=[]
    for executable in ['python3','boltz','RNAcofold','RNAfold','bowtie2-build','bowtie2','gatk','vep','tleap','samtools']:
        env.append({'executable':executable,'path':shutil.which(executable) or '', 'installed':bool(shutil.which(executable))})
    for module in ['yaml','gemmi','openmm','mdtraj','pysam','sklearn','numpy','joblib']:
        try:__import__(module);ok=True
        except ImportError:ok=False
        env.append({'python_module':module,'installed':ok})
    jsonout(args.out,{'inventory':env,'notes':['External bioinformatics packages are optional only until their stage is requested','Pretrained Alu ADAR models do not validate modified-guide editing percentages','Chromosome/annotation genome builds and exact ADAR1 isoform must be configured consistently']})
    for entry in env:print(entry)

def preflight(args):
    """Check every input, executable and resource that a configured run will touch, before it starts."""
    import yaml,importlib.util
    config=path_required(args.config,'--config');cfg=yaml.safe_load(config.read_text())
    results=[]
    def add(level,name,detail):results.append({'level':level,'check':name,'detail':detail})
    def need_file(section,key,label,required=True):
        value=(cfg.get(section) or {}).get(key)
        if not value:
            if required:add('FAIL',label,f'{section}.{key} is not set')
            return None
        path=Path(resolve(value))
        if path.is_file():add('PASS',label,str(path));return path
        add('FAIL',label,f'missing file {path}');return None
    out=Path(resolve(cfg.get('out','results/full_e2e')))
    if out.exists() and any(out.iterdir()):add('FAIL','fresh output directory',f'{out} is not empty; choose a new out')
    else:add('PASS','fresh output directory',str(out))
    free=shutil.disk_usage(out.parent if out.parent.exists() else ROOT).free/1e9
    add('PASS' if free>=args.min_free_gb else 'FAIL','free disk',f'{free:.0f} GB free, need {args.min_free_gb:.0f} GB')
    d=cfg.get('discovery',{});sm=cfg.get('somatic',{})
    variants=need_file('discovery','variants_csv','variant table',required=not (sm.get('annotated_vcf') or sm.get('vep_cache')))
    cds=need_file('discovery','cds_fasta','KRAS CDS FASTA',required=False)
    if not cds:add('WARN','KRAS CDS FASTA','not set; discovery will call the Ensembl REST service')
    if variants and cds:
        spec=importlib.util.spec_from_file_location('adar_pipeline',ROOT/'scripts'/'adar_pipeline.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        try:
            sequence=mod.loadfasta(cds)
            if len(sequence)%3 or sequence[:3]!='ATG' or sequence[33:36] not in ('GGT','GGC','GGA','GGG'):raise ValueError(f'CDS sanity check failed, codon12={sequence[33:36]}')
            add('PASS','KRAS CDS sanity',f'{len(sequence)} nt, codon 12 {sequence[33:36]}')
            bad=[]
            for row in mod.variants(variants):
                try:mod.edit_variant(sequence,[x.strip() for x in row['hgvs_c'].split(';') if x.strip()])
                except Exception as e:bad.append(f"{row['variant_id']}: {e}")
            add('FAIL' if bad else 'PASS','variant reference bases','; '.join(bad) if bad else 'every HGVS reference base matches the CDS')
        except Exception as e:add('FAIL','KRAS CDS sanity',str(e))
        if str(d.get('transcript','')).upper()=='ENST00000256078':add('WARN','transcript','ENST00000256078 is KRAS-4A; MANE Select is ENST00000311936 (KRAS-4B)')
    of=cfg.get('offtarget',{})
    if of.get('transcriptome_fasta'):
        need_file('offtarget','transcriptome_fasta','transcriptome FASTA')
        if of.get('gene_map'):need_file('offtarget','gene_map','transcript to gene map')
        else:add('WARN','transcript to gene map','not set; KRAS on-target hits cannot be separated from off-target hits')
        for exe in ('bowtie2','bowtie2-build'):add('PASS' if shutil.which(exe) else 'FAIL',exe,shutil.which(exe) or 'not on PATH')
    mod_cfg=cfg.get('structure',{})
    if mod_cfg.get('adar_fasta'):
        need_file('structure','adar_fasta','ADAR FASTA')
        add('PASS' if shutil.which('boltz') else 'FAIL','boltz',shutil.which('boltz') or 'not on PATH')
        weights=Path.home()/'.boltz'/'boltz2_conf.ckpt'
        add('PASS' if weights.is_file() else 'WARN','Boltz-2 weights',str(weights) if weights.is_file() else 'not cached; first prediction downloads about 6 GB')
        smi=shutil.which('nvidia-smi')
        add('PASS' if smi else 'WARN','GPU',smi or 'nvidia-smi not found')
        if mod_cfg.get('use_msa_server'):add('WARN','MSA server','ADAR protein sequence is sent to the public ColabFold MSA server during prediction')
    for exe in ('RNAcofold','RNAfold'):add('PASS' if shutil.which(exe) else 'FAIL',exe,shutil.which(exe) or 'not on PATH')
    if cfg.get('forcefield',{}).get('build_specs_csv'):add('PASS' if shutil.which('tleap') else 'FAIL','tleap',shutil.which('tleap') or 'not on PATH')
    if cfg.get('md',{}).get('jobs_csv') or cfg.get('forcefield',{}).get('build_specs_csv'):
        try:
            import openmm;names=[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())]
            add('PASS' if 'CUDA' in names else 'WARN','OpenMM platforms',','.join(names))
        except ImportError:add('FAIL','OpenMM','not importable')
    if sm.get('tumor_bam') or sm.get('tumor_fastq',{}).get('read1') or sm.get('vep_cache'):
        for exe in ('gatk','samtools','bcftools','vep'):add('PASS' if shutil.which(exe) else 'FAIL',exe,shutil.which(exe) or 'not on PATH')
    for r in results:print(f"{r['level']:<5} {r['check']:<28} {r['detail']}")
    failed=[r for r in results if r['level']=='FAIL'];warned=[r for r in results if r['level']=='WARN']
    print(f"\n{len(results)-len(failed)-len(warned)} passed, {len(warned)} warnings, {len(failed)} failed")
    if failed:raise GateError('Preflight failed; fix the FAIL lines before running.')


if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
    a=sp.add_parser('run');a.add_argument('--config',default='config/full_project.yaml');a.add_argument('--execute-heavy',action='store_true')
    a=sp.add_parser('preflight');a.add_argument('--config',default='config/full_project.yaml');a.add_argument('--min-free-gb',type=float,default=20.0)
    a=sp.add_parser('doctor');a.add_argument('--out',default='results/full_doctor.json')
    a=p.parse_args();{'run':run,'doctor':doctor,'preflight':preflight}[a.cmd](a)
