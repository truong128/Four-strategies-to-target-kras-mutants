from __future__ import annotations
import bz2, gzip, math, os, yaml, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

DESCRIPTOR_ORDER = ["MW","cLogP","TPSA","HBD","HBA","RotB","RingCount",
                    "HeavyAtoms","FormalCharge","FractionCSP3"]

def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)

def open_text(path, mode="rt"):
    p=str(path)
    if p.endswith(".gz"): return gzip.open(p, mode)
    if p.endswith(".bz2"): return bz2.open(p, mode)
    return open(p, mode)

def library_paths(path):
    suffixes=(".smi",".smiles",".cxsmiles",".cxsmiles.bz2",".cxsmiles.gz",".smi.gz",".smi.bz2",".csv",".tsv",".csv.gz")
    if os.path.isdir(path):
        out=[]
        for root, dirs, files in os.walk(path):
            for name in files:
                if name.startswith(".") or name.endswith(".part"): continue
                if name.endswith(suffixes): out.append(os.path.join(root, name))
        return sorted(out)
    return [str(path)]

def iter_smiles_ids(path):
    for fp in library_paths(path):
        with open_text(fp) as fh:
            for line in fh:
                line=line.strip()
                if not line or line.startswith("#"): continue
                if "\t" in line:
                    fields=line.split("\t")
                    smi=fields[0].split()[0]
                    ident=(fields[1] if len(fields)>1 else fields[0]).split(",")[0].strip()
                    if smi=="smiles": continue
                    if smi and ident: yield smi, ident
                    continue
                parts=line.split()
                if len(parts)>=2: yield parts[0], parts[1]

def passes_library_filter(desc, cfg):
    return (desc["MW"]<=cfg["max_mw"] and desc["cLogP"]<=cfg["max_clogp"] and
            desc["TPSA"]<=cfg["max_tpsa"] and desc["RotB"]<=cfg["max_rotatable_bonds"] and
            abs(desc["FormalCharge"])<=cfg["allowed_formal_charge_abs"])

def standardize_smiles(smiles):
    m=Chem.MolFromSmiles(str(smiles))
    if m is None: return None,None
    try: Chem.SanitizeMol(m)
    except Exception: return None,None
    return Chem.MolToSmiles(m,canonical=True,isomericSmiles=True),m

def basic_descriptors(m):
    return dict(
        MW=Descriptors.MolWt(m), cLogP=Crippen.MolLogP(m),
        TPSA=rdMolDescriptors.CalcTPSA(m), HBD=Lipinski.NumHDonors(m),
        HBA=Lipinski.NumHAcceptors(m), RotB=Lipinski.NumRotatableBonds(m),
        RingCount=Lipinski.RingCount(m), HeavyAtoms=m.GetNumHeavyAtoms(),
        FormalCharge=Chem.GetFormalCharge(m), FractionCSP3=rdMolDescriptors.CalcFractionCSP3(m))

def fingerprint_array(m,radius=2,nbits=2048,include_descriptors=True,desc=None):
    fp=AllChem.GetMorganFingerprintAsBitVect(m,int(radius),nBits=int(nbits))
    arr=np.zeros((int(nbits),),dtype=np.float32)
    from rdkit.DataStructs import ConvertToNumpyArray
    ConvertToNumpyArray(fp,arr)
    if include_descriptors:
        d=desc if desc is not None else basic_descriptors(m)
        arr=np.concatenate([arr,np.asarray([d[k] for k in DESCRIPTOR_ORDER],dtype=np.float32)])
    return arr

def fingerprint_and_desc(smiles,radius=2,nbits=2048,include_descriptors=True):
    m=Chem.MolFromSmiles(str(smiles))
    if m is None: return None
    return fingerprint_array(m,radius,nbits,include_descriptors)

def scaffold_smiles(smiles):
    m=Chem.MolFromSmiles(str(smiles))
    if m is None: return ""
    s=MurckoScaffold.MurckoScaffoldSmiles(mol=m,includeChirality=False)
    return s if s else Chem.MolToSmiles(m,canonical=True)

def read_smiles_table(path):
    rows=[]
    with open_text(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            p=line.split()
            if len(p)>=2: rows.append((p[0],p[1]))
    return pd.DataFrame(rows,columns=["smiles","id"])

def read_any_table(path):
    p=str(path)
    if ".csv" in p: return pd.read_csv(p)
    if ".tsv" in p: return pd.read_csv(p,sep="\t")
    return read_smiles_table(p)

def enrichment_factor(y,score,frac=0.01):
    y=np.asarray(y).astype(int); s=np.asarray(score); n=len(y)
    if n==0 or y.sum()==0: return float("nan")
    k=max(1,int(math.ceil(frac*n))); idx=np.argsort(-s)[:k]
    return float(y[idx].mean()/y.mean())

def recall_at_fraction(y,score,frac=0.01):
    y=np.asarray(y).astype(int); s=np.asarray(score)
    if y.sum()==0: return float("nan")
    k=max(1,int(math.ceil(frac*len(y)))); idx=np.argsort(-s)[:k]
    return float(y[idx].sum()/y.sum())

def model_score(model,X):
    if hasattr(model,"predict_proba"):
        return np.asarray(model.predict_proba(X))[:,1]
    if hasattr(model,"decision_function"):
        z=np.asarray(model.decision_function(X))
        return 1/(1+np.exp(-z))
    return np.asarray(model.predict(X),dtype=float)
