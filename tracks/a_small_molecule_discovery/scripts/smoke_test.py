#!/usr/bin/env python
from pathlib import Path
import subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from common import fingerprint_and_desc,scaffold_smiles,enrichment_factor
for s in ["CCO","CCN","c1ccccc1","CC(=O)O","c1ccncc1"]:
    x=fingerprint_and_desc(s,2,256,True);assert x is not None and len(x)==266
assert scaffold_smiles("c1ccccc1O")
assert enrichment_factor([0,1,0,1],[.1,.9,.2,.8],.5)>1
from benchmark_16_models import CLASSICAL, NEURAL, build_classifiers
assert len(CLASSICAL) == 16 and len(NEURAL) == 4
ml = {"dl_epochs": 1, "dl_batch_size": 64, "dl_patience": 1, "dl_learning_rate": 1e-3, "dl_max_rows": 300}
models = build_classifiers(1, ml)
assert set(models) == set(CLASSICAL) | set(NEURAL)
import numpy as np
rng = np.random.default_rng(0)
X = rng.normal(size=(240, 48)).astype(np.float32)
y = np.zeros(240, dtype=int); y[:30] = 1
for name in NEURAL:
    p = models[name].fit(X, y).predict_proba(X[:8])
    assert p.shape == (8, 2)
import joblib, tempfile, os
blob = tempfile.NamedTemporaryFile(suffix=".joblib", delete=False)
blob.close()
joblib.dump(models["realmlp"], blob.name)
loaded = joblib.load(blob.name)
assert loaded.predict_proba(X[:4]).shape == (4, 2)
os.remove(blob.name)
for p in (ROOT/"scripts").glob("*.py"):
    q=subprocess.run([sys.executable,"-m","py_compile",str(p)],capture_output=True,text=True)
    if q.returncode:raise RuntimeError(f"{p.name}: {q.stderr}")
print("SMOKE TEST PASSED")
print("Python scripts compiled successfully; external docking/ADMET/Boltz/Amber executables were not invoked.")
