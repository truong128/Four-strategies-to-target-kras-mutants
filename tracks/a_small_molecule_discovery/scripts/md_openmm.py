#!/usr/bin/env python
"""Optional explicit-solvent OpenMM protein-ligand MD driver."""
import argparse
from openmm import unit,LangevinMiddleIntegrator,MonteCarloBarostat
from openmm.app import PDBFile,Modeller,ForceField,Simulation,DCDReporter,StateDataReporter,PME,HBonds
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule
ap=argparse.ArgumentParser();ap.add_argument("--protein",required=True);ap.add_argument("--ligand",required=True);ap.add_argument("--ns",type=float,default=5);ap.add_argument("--output-prefix",required=True)
a=ap.parse_args();p=PDBFile(a.protein);lig=Molecule.from_file(a.ligand);lt=lig.to_topology().to_openmm();lp=lig.conformers[0].to_openmm()
m=Modeller(p.topology,p.positions);m.add(lt,lp);ff=ForceField("amber14-all.xml","amber14/tip3p.xml");gen=SMIRNOFFTemplateGenerator(molecules=lig, forcefield="openff-2.2.1");ff.registerTemplateGenerator(gen.generator)
m.addHydrogens(ff);m.addSolvent(ff,padding=1*unit.nanometer,ionicStrength=.15*unit.molar)
sys=ff.createSystem(m.topology,nonbondedMethod=PME,nonbondedCutoff=1*unit.nanometer,constraints=HBonds);sys.addForce(MonteCarloBarostat(1*unit.bar,300*unit.kelvin,25))
integ=LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,.002*unit.picoseconds);sim=Simulation(m.topology,sys,integ);sim.context.setPositions(m.positions)
sim.minimizeEnergy(maxIterations=5000);sim.context.setVelocitiesToTemperature(300*unit.kelvin);sim.step(50000);steps=int(a.ns*500000)
sim.reporters.append(DCDReporter(a.output_prefix+".dcd",5000));sim.reporters.append(StateDataReporter(a.output_prefix+".log",5000,step=True,potentialEnergy=True,temperature=True,density=True,progress=True,totalSteps=steps,remainingTime=True,speed=True));sim.step(steps)
st=sim.context.getState(getPositions=True)
with open(a.output_prefix+"_final.pdb","w") as f:PDBFile.writeFile(m.topology,st.getPositions(),f)
