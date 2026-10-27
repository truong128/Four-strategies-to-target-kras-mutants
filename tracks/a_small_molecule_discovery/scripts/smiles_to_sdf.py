#!/usr/bin/env python
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
ap=argparse.ArgumentParser()
ap.add_argument("--smiles",required=True); ap.add_argument("--id",required=True)
ap.add_argument("--output",required=True); a=ap.parse_args()
m=Chem.MolFromSmiles(a.smiles)
if m is None: raise SystemExit("Invalid SMILES")
m=Chem.AddHs(m); p=AllChem.ETKDGv3(); p.randomSeed=42
if AllChem.EmbedMolecule(m,p)!=0: raise SystemExit("3D embedding failed")
try: AllChem.MMFFOptimizeMolecule(m,maxIters=250)
except Exception: AllChem.UFFOptimizeMolecule(m,maxIters=250)
m.SetProp("_Name",a.id); w=Chem.SDWriter(a.output); w.write(m); w.close()
