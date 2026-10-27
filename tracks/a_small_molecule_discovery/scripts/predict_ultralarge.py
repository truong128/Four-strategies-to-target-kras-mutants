#!/usr/bin/env python
"""Score a streamed library in parallel and keep only the top hits.

Each worker reads its own compressed files, so bzip2 decompression and
fingerprinting both use many cores. Only the top screen_top_k molecules
that pass the conformal gate are written.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
import argparse, csv, gzip, heapq, multiprocessing as mp, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from rdkit import Chem
from common import (load_yaml, library_paths, iter_smiles_ids, open_text, fingerprint_array,
                    standardize_smiles, basic_descriptors, passes_library_filter, model_score)

def _pvals(p, A0, A1):
    i0=np.searchsorted(A0, p, side="left"); i1=np.searchsorted(A1, 1-p, side="left")
    return (len(A0)-i0+1)/(len(A0)+1), (len(A1)-i1+1)/(len(A1)+1)

def _init(bundle, config):
    import joblib
    global B, M, LIB, F, PROD, ISO, A0, A1, EPS, K, CHUNK
    B=joblib.load(bundle); cfg=load_yaml(config); M=cfg["ml"]; LIB=cfg["library"]
    F=B["feature_config"]; PROD=B["production"]; ISO=B["isotonic"]
    A0=np.asarray(B["alpha0"]); A1=np.asarray(B["alpha1"]); EPS=float(B["epsilon"])
    K=int(M.get("screen_top_k", cfg["redock"]["max_hits"]))
    CHUNK=int(M.get("inference_chunk_size", 4000))
    try:
        import torch
        torch.set_num_threads(1)
    except Exception:
        pass

def _consider(heap, row):
    # heap entries are (probability, tie, smiles, id, p0, p1). Smallest probability is dropped.
    prob=row[0]
    if len(heap)<K:
        heapq.heappush(heap, row)
    elif prob>heap[0][0]:
        heapq.heapreplace(heap, row)

def _flush(batch, heap, tie):
    if not batch: return tie
    X=[]; meta=[]
    for smi, mid, x in batch:
        X.append(x); meta.append((smi, mid))
    batch.clear()
    X=np.asarray(X, dtype=np.float32)
    raw=np.mean([model_score(m, X) for m in PROD], axis=0) if isinstance(PROD, list) else model_score(PROD, X)
    pr=np.asarray(ISO.transform(raw), dtype=float)
    p0,p1=_pvals(pr, A0, A1)
    keep=(p1>EPS)&(p0<=EPS)&(pr>=float(M.get("keep_min_probability", 0)))
    for (smi, mid), q, a0, a1, z in zip(meta, pr, p0, p1, keep):
        if not z: continue
        tie+=1
        _consider(heap, (float(q), tie, smi, mid, float(a0), float(a1)))
    return tie

def _iter_prefiltered(fp, lib):
    """Yield (smiles, id) and skip rows whose precomputed Enamine MW, sLogP, TPSA
    or rotatable-bond columns already fail the library filter. MW, sLogP and TPSA match
    RDKit to rounding, so a 0.01 margin is safe. Enamine rotatable bonds can exceed
    RDKit's count by 1, so that column gets a margin of 1."""
    with open_text(fp) as fh:
        header=fh.readline().rstrip("\n").split("\t")
        col={name:i for i, name in enumerate(header)}
        need=("MW", "sLogP", "TPSA", "RotBonds")
        if header[:2]!=["smiles", "id"] or any(k not in col for k in need):
            for rec in iter_smiles_ids(fp):
                yield rec
            return
        iw, il, it, ir = (col[k] for k in need)
        mw_max=lib["max_mw"]+0.01; lp_max=lib["max_clogp"]+0.01
        tp_max=lib["max_tpsa"]+0.01; rb_max=lib["max_rotatable_bonds"]+1
        for line in fh:
            f=line.rstrip("\n").split("\t")
            try:
                if (float(f[iw])>mw_max or float(f[il])>lp_max or
                        float(f[it])>tp_max or float(f[ir])>rb_max):
                    continue
            except (ValueError, IndexError):
                continue
            yield f[0], f[1]

def _read_files(files, queue, batch_n, lib=None):
    batch=[]
    for fp in files:
        for rec in (_iter_prefiltered(fp, lib) if lib else iter_smiles_ids(fp)):
            batch.append(rec)
            if len(batch)>=batch_n:
                queue.put(batch); batch=[]
    if batch: queue.put(batch)

def _score_queue(queue, out, bundle, config):
    _init(bundle, config)
    heap=[]; tie=0; seen=0; feat=[]
    while True:
        batch=queue.get()
        if batch is None: break
        seen+=len(batch)
        for smi, mid in batch:
            m=Chem.MolFromSmiles(smi)
            if m is None: continue
            d=basic_descriptors(m)
            if not passes_library_filter(d, LIB): continue
            x=fingerprint_array(m, F["radius"], F["nbits"], F["include_descriptors"], d)
            feat.append((smi, mid, x))
            if len(feat)>=CHUNK:
                tie=_flush(feat, heap, tie)
        tie=_flush(feat, heap, tie)
    out.put((heap, seen))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--bundle", required=True)
    ap.add_argument("--config", required=True); ap.add_argument("--output", required=True)
    a=ap.parse_args(); cfg=load_yaml(a.config); M=cfg["ml"]
    files=library_paths(a.input)
    cores=os.cpu_count() or 1
    workers=max(1, min(int(M.get("inference_workers") or 12), cores-1))
    readers=max(1, min(int(M.get("reader_workers") or 4), len(files), cores-workers))
    groups=[[] for _ in range(readers)]
    for i, fp in enumerate(files): groups[i%readers].append(fp)
    ctx=mp.get_context("spawn")
    queue=ctx.Queue(maxsize=48); out=ctx.Queue()
    rprocs=[ctx.Process(target=_read_files, args=(g, queue, 2000, cfg["library"])) for g in groups if g]
    wprocs=[ctx.Process(target=_score_queue, args=(queue, out, a.bundle, a.config)) for _ in range(workers)]
    for p in rprocs+wprocs: p.start()
    for p in rprocs: p.join()
    for _ in wprocs: queue.put(None)
    heaps=[]; seen=0
    for _ in wprocs:
        h, n=out.get(); heaps.append(h); seen+=n
    for p in wprocs: p.join()
    K=int(M.get("screen_top_k", cfg["redock"]["max_hits"]))
    merged=[]
    for h in heaps: merged.extend(h)
    best=sorted(merged, key=lambda r: r[0], reverse=True)[:K]
    op=gzip.open if str(a.output).endswith(".gz") else open
    with op(a.output, "wt", newline="") as o:
        w = csv.writer(o)
        w.writerow(["smiles", "id", "ml_probability", "p0", "p1", "deltaP"])
        for q, tie, smi, mid, a0, a1 in best:
            cs, _ = standardize_smiles(smi)
            if cs is None: continue
            w.writerow([cs, mid, f"{q:.8g}", f"{a0:.8g}", f"{a1:.8g}", f"{(a1 - a0):.8g}"])
    print(f"scored_lines {seen} hits_kept {len(best)}")

if __name__=="__main__":
    main()
