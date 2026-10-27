"""Amber tleap executors for validated modified nucleotide residue libraries.
This program NEVER hallucinates partial charges or PS chirality. Builds run-ready
solvated systems from supplied .off/.lib/.frcmod files and user PDB structures.
"""
import argparse, csv, json, math, subprocess
from pathlib import Path
from full.common import rows, jsonout, csvout, path_required, GateError, run, which, bounded_id, digest

def gen_tleap(args):
    input_structure=path_required(args.pdb,'--pdb');out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    ff_files=[path_required(f,'--frcmod') for f in args.frcmod]
    lib_files=[path_required(f,'--lib') for f in args.lib]
    if args.modified and not (ff_files and lib_files):raise GateError('Modified residues require both validated Amber libraries and frcmod files')
    if args.modified and not args.parameter_validation_json:raise GateError('Modified system requires explicit parameter validation evidence JSON')
    if args.parameter_validation_json:
        p=path_required(args.parameter_validation_json,'parameter validation JSON');v=json.loads(p.read_text())
        required=('validated','stereochemistry_resolved','charge_and_torsion_checked','atom_mapping_checked')
        if any(v.get(x) is not True for x in required):raise GateError('Parameter validation evidence does not pass: '+str(required))
        if args.modified:
            if v.get('pdb_sha256')!=digest(input_structure):raise GateError('Validated PDB hash missing/mismatched; evidence not bound to input structure')
            hashes=v.get('lib_sha256',{})
            if not isinstance(hashes,dict) or not hashes:raise GateError('Validated residue library hash map missing')
            for file in lib_files:
                if hashes.get(file.name)!=digest(file):raise GateError('Residue library hash mismatch: '+file.name)
            hashes=v.get('frcmod_sha256',{})
            if not isinstance(hashes,dict) or not hashes:raise GateError('Validated frcmod hash map missing')
            for file in ff_files:
                if hashes.get(file.name)!=digest(file):raise GateError('Parameter file hash mismatch: '+file.name)
    if not 5<=args.padding_A<=30 or not 0<=args.salt_molar<1.5:raise GateError('Unreasonable padding or salt concentration')
    base=['source leaprc.protein.ff14SB','source leaprc.RNA.OL3','source leaprc.DNA.OL15','source leaprc.water.tip3p']
    for file in lib_files:base+=['loadoff '+json.dumps(str(file))]
    for file in ff_files:base+=['loadamberparams '+json.dumps(str(file))]
    base += ['sys = loadpdb '+json.dumps(str(input_structure)),'check sys',f'solvateoct sys TIP3PBOX {args.padding_A}', 'addIons2 sys Na+ 0','addIons2 sys Cl- 0']
    # First-stage water count gives a physically motivated salt-pair estimate:
    # concentration mol/L approximately N_pairs/N_water*55.5M.
    stage1=base+[f'savepdb sys {json.dumps(str((out/"stage1_water.pdb").resolve()))}','quit']
    (out/'stage1.leap').write_text('\n'.join(stage1)+'\n')
    if args.execute:
        run([which('tleap'),'-f',str((out/'stage1.leap').resolve())],out/'stage1_tleap.log',True)
        waterfile=path_required(out/'stage1_water.pdb','First-stage solvated PDB')
        water_keys=set()
        with waterfile.open() as fp:
            for line in fp:
                if line.startswith(('ATOM  ','HETATM')) and line[17:20].strip() in ('WAT','HOH','TIP3'):
                    water_keys.add((line[21:22],line[22:26],line[26:27]))
        nwater=len(water_keys)
        if nwater<10:raise GateError('Tleap first stage did not produce a water box: inspect stage1 logs')
        pairs=args.salt_pairs if args.salt_pairs is not None else round(nwater*args.salt_molar/55.5)
    else:
        nwater=None
        pairs=args.salt_pairs
    build=base[:]
    if pairs is not None and pairs>0:
        build += [f'addIonsRand sys Na+ {pairs}',f'addIonsRand sys Cl- {pairs}']
    elif pairs is None:
        build += ['# AUTO additional salt pairs from stage1 water number at execute time']
    build += [f'saveamberparm sys {json.dumps(str((out/"system.prmtop").resolve()))} {json.dumps(str((out/"system.inpcrd").resolve()))}',f'savepdb sys {json.dumps(str((out/"solvated.pdb").resolve()))}','quit']
    (out/'build.leap').write_text('\n'.join(build)+'\n')
    jsonout(out/'parameterization_scope.json',{'input_pdb':str(input_structure),'modified':args.modified,'salt_target_M':args.salt_molar,'salt_pairs':pairs,'n_waters_after_neutralization':nwater,'approx_bulk_salt_M':pairs/nwater*55.5 if (nwater and pairs is not None) else 'NOT_YET_CALCULATED','temperature_K':310.15,'note':'Ion strength is approximate. Validate using actual solvent box and ions. No constant-pH simulation. Force field quality is user validated independently; tleap success is not proof.'})
    if args.execute:
        run([which('tleap'),'-f',str((out/'build.leap').resolve())],out/'tleap.log',True)
        if not (out/'system.prmtop').exists() or not (out/'system.inpcrd').exists():raise GateError('tleap completed without output topology/coords')
        log=(out/'tleap.log').read_text(errors='replace')
        if 'FATAL' in log.upper() or 'could not find atom' in log.lower():raise GateError('tleap output indicates unresolved errors')

def main():
    p=argparse.ArgumentParser();p.add_argument('--pdb',required=True);p.add_argument('--out',required=True)
    p.add_argument('--lib',nargs='*',default=[]);p.add_argument('--frcmod',nargs='*',default=[])
    p.add_argument('--modified',action='store_true');p.add_argument('--parameter-validation-json');p.add_argument('--padding-A',type=float,default=12.0);p.add_argument('--salt-pairs',type=int);p.add_argument('--salt-molar',type=float,default=0.15);p.add_argument('--execute',action='store_true')
    a=p.parse_args();gen_tleap(a)
if __name__=='__main__':main()
