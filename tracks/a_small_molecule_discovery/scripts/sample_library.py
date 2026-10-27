#!/usr/bin/env python
"""Draw a random training sample from a very large library.

Directory inputs with library.sample_expected_lines set use independent
Bernoulli sampling per line (geometric skipping). Several compressed files
stream at once through lbzip2 (multi-core bzip2), so a 328 GB Enamine REAL
read uses most of the machine. Everything else uses an exact reservoir.
"""
import argparse, math, multiprocessing as mp, os, random, shutil, subprocess, sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (load_yaml, library_paths, iter_smiles_ids, standardize_smiles,
                    basic_descriptors, passes_library_filter)

CHUNK_BYTES = 1 << 26


def find_exe(names, configured=None):
    if configured and os.path.exists(configured):
        return configured
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def parse_line(raw):
    line = raw.decode("utf-8", "replace").strip()
    if not line or line.startswith("#"):
        return None
    if "\t" in line:
        f = line.split("\t")
        smi = f[0].split()[0] if f[0].split() else ""
        ident = (f[1] if len(f) > 1 else f[0]).split(",")[0].strip()
        if smi == "smiles" or not smi or not ident:
            return None
        return smi, ident
    parts = line.split()
    if len(parts) >= 2 and parts[0] != "smiles":
        return parts[0], parts[1]
    return None


def open_stream(fp, threads, lbzip2):
    if fp.endswith(".bz2"):
        cmd = [lbzip2, "-dc", "-n", str(threads), fp] if lbzip2 else ["bzip2", "-dc", fp]
    elif fp.endswith(".gz"):
        pigz = shutil.which("pigz")
        cmd = [pigz, "-dc", "-p", str(threads), fp] if pigz else ["gzip", "-dc", fp]
    else:
        cmd = ["cat", fp]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=CHUNK_BYTES)


def sample_file(job):
    fp, p, seed, threads, lbzip2 = job
    rng = random.Random(seed)
    lp = math.log1p(-p)

    def skip():
        return int(math.log(1.0 - rng.random()) / lp)

    picks = []; base = 0; nxt = skip(); tail = b""
    proc = open_stream(fp, threads, lbzip2)
    while True:
        buf = proc.stdout.read(CHUNK_BYTES)
        if not buf:
            break
        buf = tail + buf
        cut = buf.rfind(b"\n")
        if cut < 0:
            tail = buf
            continue
        tail = buf[cut + 1:]
        body = buf[:cut]
        n = body.count(b"\n") + 1
        if nxt - base < n:
            lines = body.split(b"\n")
            while nxt - base < n:
                picks.append(lines[nxt - base])
                nxt += 1 + skip()
        base += n
    if tail:
        if nxt == base:
            picks.append(tail)
        base += 1
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"decompression of {fp} failed with exit code {rc}")
    recs = [r for r in map(parse_line, picks) if r is not None]
    return fp, base, recs


def _clean(args):
    chunk, lib = args
    out = []
    for smi, mid in chunk:
        cs, m = standardize_smiles(smi)
        if m is None:
            continue
        if not passes_library_filter(basic_descriptors(m), lib):
            continue
        out.append((cs, mid))
    return out


def clean_parallel(records, lib, n, workers):
    size = 4000
    chunks = [(records[i:i + size], lib) for i in range(0, len(records), size)]
    kept = []
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
        for part in pool.map(_clean, chunks):
            kept.extend(part)
            if len(kept) >= n:
                break
    return kept[:n]


def run_parallel(path, c, n, draw, rng):
    files = library_paths(path)
    cores = os.cpu_count() or 1
    expected = float(c["sample_expected_lines"])
    p = min(0.5, draw / expected)
    lbzip2 = find_exe(["lbzip2"], c.get("lbzip2_exe"))
    concurrent = max(1, min(int(c.get("sample_file_parallel", 4)), len(files)))
    threads = max(1, int(c.get("sample_decompress_threads", max(1, (cores - concurrent) // concurrent))))
    if lbzip2 is None:
        print("WARNING lbzip2 not found; falling back to single-core bzip2", flush=True)
    print(f"parallel sampling: files {len(files)} concurrent {concurrent} threads_per_file {threads} "
          f"p {p:.3e} lbzip2 {lbzip2}", flush=True)
    seed0 = int(c.get("reservoir_seed", 42))
    order = sorted(range(len(files)), key=lambda i: -os.path.getsize(files[i]))
    jobs = [(files[i], p, seed0 + i, threads, lbzip2) for i in order]
    records = []; total = 0
    with ProcessPoolExecutor(max_workers=concurrent, mp_context=mp.get_context("spawn")) as pool:
        for fp, lines, recs in pool.map(sample_file, jobs):
            total += lines; records.extend(recs)
            print(f"done {os.path.basename(fp)} lines {lines} picked {len(recs)}", flush=True)
    rng.shuffle(records)
    return records, total


def run_reservoir(path, draw, rng):
    R = []; i = 0
    for rec in iter_smiles_ids(path):
        i += 1
        if len(R) < draw:
            R.append(rec)
        else:
            j = rng.randrange(i)
            if j < draw:
                R[j] = rec
    return R, i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--output", required=True)
    ap.add_argument("--config", required=True); a = ap.parse_args()
    c = load_yaml(a.config)["library"]; n = int(c["training_sample_n"])
    draw = n * int(c.get("sample_oversample", 4))
    rng = random.Random(int(c.get("reservoir_seed", 42)))
    use_parallel = os.path.isdir(a.input) and float(c.get("sample_expected_lines") or 0) > 0
    if use_parallel:
        records, total = run_parallel(a.input, c, n, draw, rng)
    else:
        records, total = run_reservoir(a.input, draw, rng)
    workers = max(1, (os.cpu_count() or 1) - 1)
    kept = clean_parallel(records, c, n, workers)
    if len(kept) < n:
        print(f"WARNING only {len(kept)} molecules passed the filters (target {n})", flush=True)
    with open(a.output, "w") as fo:
        for cs, mid in kept:
            fo.write(f"{cs}\t{mid}\n")
    print(f"sampled {len(kept)} of {total} lines", flush=True)


if __name__ == "__main__":
    main()
