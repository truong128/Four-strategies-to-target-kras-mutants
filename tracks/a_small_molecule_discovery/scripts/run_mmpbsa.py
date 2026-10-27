#!/usr/bin/env python
"""Strict AmberTools MMPBSA.py wrapper. Requires Amber-compatible topology/trajectory."""
import argparse,subprocess,shutil
ap=argparse.ArgumentParser()
ap.add_argument("--exe",default="MMPBSA.py");ap.add_argument("--input",default="mmpbsa.in")
for x in ["complex","receptor","ligand","trajectory","output"]:ap.add_argument("--"+x,required=True)
a=ap.parse_args()
if shutil.which(a.exe) is None:raise SystemExit("MMPBSA.py not found")
raise SystemExit(subprocess.run([a.exe,"-O","-i",a.input,"-cp",a.complex,"-rp",a.receptor,"-lp",a.ligand,"-y",a.trajectory,"-o",a.output]).returncode)
