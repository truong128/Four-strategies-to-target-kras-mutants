"""Executable Boltz-2 canonical protein/RNA predictions and structural post-analysis.
Canonical RNA models are STARTING POSES; do not label them chemically modified or catalytic.
"""
import argparse, math, json, re, shutil
from pathlib import Path
from full.common import rows, jsonout, csvout, path_required, GateError, run, which, bounded_id, digest

def read_fasta(path):
    chars=[]
    for l in Path(path).read_text().splitlines():
        if l.startswith('>'):continue
        chars.append(l.strip())
    return ''.join(chars).upper()

def input_yaml(args):
    try:import yaml
    except ImportError as e:raise GateError('pip install PyYAML') from e
    protein=read_fasta(path_required(args.adar_fasta,'ADAR FASTA'))
    if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',protein):raise GateError('ADAR FASTA must be unmodified single protein chain sequence')
    if args.msa!='auto' and args.msa!='empty':path_required(args.msa,'--msa')
    if args.msa=='empty' and not args.allow_single_sequence:raise GateError('MSA empty has low reliability; explicitly set --allow-single-sequence')
    guidedata=rows(path_required(args.guides,'guide table'))
    if args.variants:
        wanted=set(args.variants)
        guidedata=[r for r in guidedata if r['variant_id'] in wanted]
    if args.max_guides:guidedata=guidedata[:args.max_guides]
    if not guidedata:raise GateError('No guides left after --variants/--max-guides filtering')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    index=[]
    for r in guidedata:
        gid=bounded_id(r['guide_id']);target=r['target_rna_5to3'].upper().replace('T','U');guide=r['guide_rna_5to3'].upper().replace('T','U')
        if not re.fullmatch('[ACGU]+',target+guide):raise GateError(f'Noncanonical RNA in {gid}')
        yaml_path=out/(gid+'.yaml')
        # Boltz 2: protein chain A, RNA target B, RNA guide C. No CCD PS backbone claim.
        model={'version':1,'sequences':[{'protein':{'id':'A','sequence':protein,'msa':args.msa if args.msa=='empty' else ('empty' if args.msa=='auto' else str(Path(args.msa).resolve()))}},{'rna':{'id':'B','sequence':target}},{'rna':{'id':'C','sequence':guide}}]}
        if args.msa=='auto':del model['sequences'][0]['protein']['msa']
        yaml_path.write_text(yaml.safe_dump(model,sort_keys=False))
        index.append({'guide_id':gid,'input_yaml':str(yaml_path.resolve()),'protein_chain':'A','target_rna_chain':'B','guide_rna_chain':'C','target_A_residue_1based':int(r['target_index_0based'])+1,'guide_orphan_C_1based':int(r['guide_index_0based'])+1,'chemical_model':'UNMODIFIED_RNA_SURROGATE','status':'STRUCTURE_PREDICTION_INPUT'})
    csvout(out/'structure_index.csv',index,list(index[0]) if index else ['guide_id'])
    jsonout(out/'structure_disclaimer.json',{'scope':'Canonical RNA-protein complex pose screening','NOT':'Exact RESTORE2 2OMe/2F/PS atomic structure, RNA deamination chemistry, metal-center validation, or binding affinity','protein_source_fasta_sha256':digest(args.adar_fasta),'msa_mode':args.msa})
    return out

def predict(args):
    folder=Path(args.input_dir)
    if not folder.is_dir():raise GateError('Boltz input YAML folder not found')
    import yaml
    yfiles=sorted(folder.glob('*.yaml'))
    if not yfiles:raise GateError('No Boltz YAML inputs found')
    # Auto MSA requires server access; custom .a3m or explicitly empty MSA is an alternative.
    missing=[]
    for f in yfiles:
        spec=yaml.safe_load(f.read_text())
        for entity in spec.get('sequences',[]):
            if 'protein' in entity and 'msa' not in entity['protein'] and not args.use_msa_server:
                missing.append(f.name)
    if missing:raise GateError('Boltz protein MSA is missing. Enable --use-msa-server or provide --msa .a3m (or explicitly empty). YAML: '+','.join(missing[:4]))
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    # Boltz rejects any non-FASTA/YAML file in its input folder, so stage only the YAMLs.
    staged=out/'boltz_yaml'
    staged.mkdir(parents=True,exist_ok=True)
    for f in yfiles:shutil.copy2(f,staged/f.name)
    cmd=[which('boltz') if args.execute else 'boltz','predict',str(staged.resolve()),'--out_dir',str(out.resolve()),'--devices',str(args.devices),'--diffusion_samples',str(args.samples)]
    if args.use_msa_server:cmd.append('--use_msa_server')
    if args.no_kernels:cmd.append('--no_kernels')
    run(cmd,out/'boltz.log',args.execute)

def geometry(args):
    """Report every Boltz diffusion sample separately, never choose a pose silently."""
    try:import gemmi
    except ImportError as e:raise GateError('pip install gemmi for mmCIF/PDB geometry evaluation') from e
    idx=rows(path_required(args.index,'--index'))
    modeldir=Path(args.models)
    if not modeldir.is_dir():raise GateError('Prediction output directory missing')
    if args.map_file and len(idx)!=1:
        raise GateError('--map-file can only be used when structure_index contains exactly one guide')
    result=[]
    for g in idx:
        gid=g['guide_id']
        if args.map_file:matches=[path_required(args.map_file,'--map-file')]
        else:
            matches=sorted(set(modeldir.rglob('*'+gid+'*model_*.cif')) | set(modeldir.rglob('*'+gid+'*model_*.pdb')))
        if not matches:
            result.append({'guide_id':gid,'model_path':'','status':'NO_PREDICTION_FOUND','n_candidates':0})
            continue
        for path in matches:
            try:
                st=gemmi.read_structure(str(path));m=st[0];chains={c.name:c for c in m}
                if not {'A','B','C'}.issubset(chains):raise ValueError('expected protein A, target RNA B, guide RNA C chain IDs')
                targetpos=int(g['target_A_residue_1based'])
                target=chains['B']
                if not 1<=targetpos<=len(target):raise ValueError('Target nucleotide index outside modeled RNA chain')
                targetres=target[targetpos-1]
                if targetres.name.strip().upper() not in {'A','ADE','RA'}:raise ValueError('target RNA nucleotide is not A')
                reaction_atoms=[at for at in targetres if at.name.strip() in ('N6','C6','N1')]
                if not reaction_atoms:raise ValueError('N6/C6/N1 atoms absent at target A')
                protein_atoms=[at for res in chains['A'] for at in res if at.element.name not in ('H','D')]
                if not protein_atoms:raise ValueError('Protein chain contains no heavy atoms')
                closest=min(at.pos.dist(pr.pos) for at in reaction_atoms for pr in protein_atoms)
                zinc=[at for chain in m for res in chain for at in res if at.element.name=='Zn']
                zdist=min((at.pos.dist(z.pos) for at in reaction_atoms for z in zinc),default=None)
                result.append({'guide_id':gid,'model_path':str(path),'status':'GEOMETRY_ONLY_UNVALIDATED_NO_CATALYTIC_INFERENCE','n_candidates':len(matches),'target_A_to_nearest_protein_Ang':round(closest,3),'targetA_to_Zn_Ang':round(zdist,3) if zdist is not None else '', 'zinc_present':bool(zinc)})
            except Exception as e:
                result.append({'guide_id':gid,'model_path':str(path),'status':'GEOMETRY_PARSE_ERROR:'+str(e)[:100],'n_candidates':len(matches)})
    columns=['guide_id','model_path','status','n_candidates','target_A_to_nearest_protein_Ang','targetA_to_Zn_Ang','zinc_present']
    csvout(args.out,result,columns)
    jsonout(str(args.out)+'.scope.json',{'meaning':'Each predicted model is listed independently. Geometry is not a probability of editing or affinity.','limits':['Boltz canonical RNA input omits backbone chemistry and ion/metal coordination verification','ADAR1 isoforms require template/domain and biological-context checks','Compare multiple replicates or conformational ensembles, never select visually appealing single predictions']})

if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
    s=sp.add_parser('prepare');s.add_argument('--guides',required=True);s.add_argument('--adar-fasta',required=True);s.add_argument('--msa',default='auto');s.add_argument('--allow-single-sequence',action='store_true');s.add_argument('--variants',nargs='*');s.add_argument('--max-guides',type=int);s.add_argument('--out',required=True)
    s=sp.add_parser('predict');s.add_argument('--input-dir',required=True);s.add_argument('--out',required=True);s.add_argument('--devices',type=int,default=1);s.add_argument('--samples',type=int,default=5);s.add_argument('--use-msa-server',action='store_true');s.add_argument('--no-kernels',action='store_true');s.add_argument('--execute',action='store_true')
    s=sp.add_parser('geometry');s.add_argument('--index',required=True);s.add_argument('--models',required=True);s.add_argument('--map-file');s.add_argument('--out',required=True)
    a=p.parse_args();{'prepare':input_yaml,'predict':predict,'geometry':geometry}[a.cmd](a)
