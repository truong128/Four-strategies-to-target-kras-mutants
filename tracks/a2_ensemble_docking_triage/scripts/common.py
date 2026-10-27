from pathlib import Path
import yaml

def root_dir():
    return Path(__file__).resolve().parents[1]

def load_config():
    root = root_dir()
    with open(root/'config'/'config.yaml') as f:
        return yaml.safe_load(f)

def ensure_dirs():
    root = root_dir()
    for p in ['work/raw_pdb','work/aligned_pdb','work/receptors','work/ligands/sdf','work/ligands/pdbqt','work/docking','results/figures','results/tables','logs']:
        (root/p).mkdir(parents=True, exist_ok=True)
    return root
