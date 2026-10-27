"""Production MD executor: minimization -> thermal equilibration -> NPT -> production.
Supports validated Amber prmtop/inpcrd files including custom modifications.
No quantum chemistry / A->I chemical reaction simulation is claimed.
"""
import argparse, json, math, os
from pathlib import Path
from full.common import GateError, path_required, jsonout, digest

def simulate(args):
    top=path_required(args.topology,'--topology');pos=path_required(args.coordinates,'--coordinates')
    if args.modified:
        evidence=path_required(args.chemistry_validation,'--chemistry-validation')
        checks=json.loads(evidence.read_text())
        for key in ('validated','stereochemistry_resolved','charge_and_torsion_checked','atom_mapping_checked'):
            if checks.get(key) is not True:raise GateError(f'Cannot run modified model without {key}=true evidence')
        for key,path in (('topology_sha256',top),('coordinates_sha256',pos)):
            if checks.get(key)!=digest(path):raise GateError(f'Modified model evidence missing/mismatched {key}')
    if args.replicates<1 or args.replicates>50:raise GateError('Replicates must be 1..50')
    if not 0<=args.restraint_k<=100000:raise GateError('Invalid positional restraint constant')
    if min(args.minimize_steps,args.equil_steps,args.production_steps)<0:raise GateError('Negative simulation steps')
    if args.production_steps<1:raise GateError('production_steps must be positive')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    metadata={'topology':str(top),'coordinates':str(pos),'topology_sha256':digest(top),'coordinates_sha256':digest(pos),'n_replicates':args.replicates,'temperature_K':args.temperature,'pressure_atm':args.pressure,'time_step_fs':args.timestep_fs,'equil_steps':args.equil_steps,'equilibration_phases':'restrained_NPT_then_unrestrained_NPT','restraint_k_kJ_mol_nm2':args.restraint_k,'production_steps':args.production_steps,'salt_molar':'DEFINED_BY_INPUT_TOPOLOGY_NOT_BY_MD','chemistry':'CUSTOM_VALIDATED' if args.modified else 'CANONICAL','execution':bool(args.execute)}
    jsonout(out/'run_protocol.json',metadata)
    if not args.execute:return
    try:
        from openmm import unit, Platform, LangevinMiddleIntegrator, MonteCarloBarostat, CustomExternalForce
        from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, PME, HBonds, DCDReporter, StateDataReporter, CheckpointReporter, PDBFile
    except ImportError as exc:raise GateError('Install OpenMM to execute: micromamba install -c conda-forge openmm') from exc
    ambertop=AmberPrmtopFile(str(top));amberpos=AmberInpcrdFile(str(pos))
    periodic=amberpos.boxVectors is not None
    if not periodic:raise GateError('Periodic solvent/ions required for NPT; coordinates contain no periodic box')
    for rep in range(1,args.replicates+1):
        d=out/f'rep_{rep:02d}';d.mkdir(parents=True,exist_ok=True)
        system=ambertop.createSystem(nonbondedMethod=PME,nonbondedCutoff=args.cutoff_nm*unit.nanometer,constraints=HBonds,rigidWater=True)
        system.addForce(MonteCarloBarostat(args.pressure*unit.atmosphere,args.temperature*unit.kelvin,25))
        # Positional solute-heavy-atom restraints for the first NPT equilibration stage.
        # They are released before unrestrained NPT equilibration and production.
        restraint=CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
        restraint.addGlobalParameter('k',args.restraint_k)
        for component in ('x0','y0','z0'):restraint.addPerParticleParameter(component)
        xyz=amberpos.positions.value_in_unit(unit.nanometer)
        restricted=0
        for ai,atom in enumerate(ambertop.topology.atoms()):
            # Restrain protein and nucleic-acid heavy atoms only, not solvent or ions.
            if atom.residue.name.upper() in ('WAT','HOH','NA','CL','NA+','CL-','K','K+'):
                continue
            if atom.element is None or atom.element.symbol in ('H','D'):
                continue
            restraint.addParticle(ai,[xyz[ai][0],xyz[ai][1],xyz[ai][2]])
            restricted+=1
        if restricted<5:raise GateError('Too few solute heavy atoms to restrain; review topology')
        system.addForce(restraint)
        integ=LangevinMiddleIntegrator(args.temperature*unit.kelvin,args.collision_per_ps/unit.picosecond,args.timestep_fs*unit.femtoseconds)
        integ.setRandomNumberSeed(args.seed+rep*104729)
        if args.platform=='auto':sim=Simulation(ambertop.topology,system,integ)
        else:
            plat=Platform.getPlatformByName(args.platform)
            sim=Simulation(ambertop.topology,system,integ,plat)
        sim.context.setPositions(amberpos.positions)
        sim.context.setPeriodicBoxVectors(*amberpos.boxVectors)
        sim.context.setVelocitiesToTemperature(args.temperature*unit.kelvin,args.seed+rep*104729)
        if args.minimize_steps:sim.minimizeEnergy(maxIterations=args.minimize_steps)
        if args.equil_steps:
            sim.reporters.append(StateDataReporter(str(d/'equil.csv'),max(1,min(args.report_interval,args.equil_steps)),step=True,potentialEnergy=True,temperature=True,volume=True,separator=','))
            restrained_steps=max(1,args.equil_steps//2)
            sim.step(restrained_steps)
            sim.context.setParameter('k',0.0)
            sim.step(args.equil_steps-restrained_steps)
            sim.reporters.clear()
        else:
            # Production must never retain the initial positional restraint.
            sim.context.setParameter('k',0.0)
        # Write starting coordinates after equilibration; file is not a pre-equilibrated structure until equil steps executed.
        state=sim.context.getState(getPositions=True)
        with (d/'equilibrated.pdb').open('w') as f:
            PDBFile.writeFile(sim.topology,state.getPositions(),f,keepIds=True)
        sim.reporters.append(DCDReporter(str(d/'production.dcd'),max(1,args.report_interval)))
        sim.reporters.append(StateDataReporter(str(d/'production.csv'),max(1,args.report_interval),step=True,potentialEnergy=True,kineticEnergy=True,temperature=True,volume=True,separator=','))
        sim.reporters.append(CheckpointReporter(str(d/'production.chk'),max(1,args.checkpoint_interval)))
        sim.step(args.production_steps)
        with (d/'final.pdb').open('w') as f:
            PDBFile.writeFile(sim.topology,sim.context.getState(getPositions=True).getPositions(),f,keepIds=True)
        sim.saveState(str(d/'final_state.xml'))
        jsonout(d/'provenance.json',{'replicate':rep,'status':'PRODUCTION_COMPLETED_NOT_CONVERGENCE_VALIDATED','production_ns':args.production_steps*args.timestep_fs*1e-6,'equilibration_ns':args.equil_steps*args.timestep_fs*1e-6,'seed':args.seed+rep*104729})

def analyze(args):
    try:import mdtraj as md
    except ImportError as e:raise GateError('Install mdtraj for MD analysis') from e
    import numpy as np
    pdb=path_required(args.topology_pdb,'--topology-pdb')
    trajectories=[Path(x) for x in args.trajectories]
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    res=[]
    for trpath in trajectories:
        t=md.load(str(path_required(str(trpath),'--trajectories')),top=str(pdb))
        solute=t.topology.select('not water and not element Na and not element Cl')
        if len(solute)<5:raise GateError('Solute selection unexpectedly small')
        if t.topology.n_bonds and t.unitcell_vectors is not None:
            t.image_molecules(inplace=True)
        sol=t.atom_slice(solute)
        sol.superpose(sol,0)
        rmsd=md.rmsd(sol,sol,0)
        rg=md.compute_rg(sol)
        res.append({'trajectory':str(trpath),'frames':t.n_frames,'rmsd_median_nm':float(np.median(rmsd)),'rmsd_p95_nm':float(np.percentile(rmsd,95)),'radius_gyration_mean_nm':float(np.mean(rg)),'warning':'RMSD alone cannot demonstrate convergence or editing; compare independent replicas and local RNA geometry'})
    from full.common import csvout
    csvout(out/'md_summary.csv',res,list(res[0]))
    jsonout(out/'analysis_scope.json',{'notes':['Analyze RNA hydrogen bonds/base-pair step parameters separately with cpptraj nastruct, after defining RNA masks for the specific structure','Do not compare raw RMSD of different systems without consistent alignment and atoms','Base flipping barrier and ADAR catalysis require enhanced sampling/QM-MM, not this classical MD executor']})

if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
    a=sp.add_parser('simulate');a.add_argument('--topology',required=True);a.add_argument('--coordinates',required=True);a.add_argument('--out',required=True)
    a.add_argument('--modified',action='store_true');a.add_argument('--chemistry-validation');a.add_argument('--execute',action='store_true')
    a.add_argument('--replicates',type=int,default=3);a.add_argument('--restraint-k',type=float,default=1000.0);a.add_argument('--minimize-steps',type=int,default=3000);a.add_argument('--equil-steps',type=int,default=250000);a.add_argument('--production-steps',type=int,default=5000000)
    a.add_argument('--report-interval',type=int,default=5000);a.add_argument('--checkpoint-interval',type=int,default=50000);a.add_argument('--seed',type=int,default=1976);a.add_argument('--temperature',type=float,default=310.15);a.add_argument('--pressure',type=float,default=1.0);a.add_argument('--collision-per-ps',type=float,default=1.0);a.add_argument('--cutoff-nm',type=float,default=0.9);a.add_argument('--timestep-fs',type=float,default=2.0);a.add_argument('--platform',default='auto',choices=['auto','CUDA','OpenCL','CPU','Reference'])
    a=sp.add_parser('analyze');a.add_argument('--trajectories',required=True,nargs='+');a.add_argument('--topology-pdb',required=True);a.add_argument('--out',required=True)
    a=p.parse_args();{'simulate':simulate,'analyze':analyze}[a.cmd](a)
