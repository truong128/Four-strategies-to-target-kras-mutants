#!/usr/bin/env python3
"""Engineering triage of designed guides for a first wet-lab pilot.

This is a transparent, rule-based ordering of candidates for experimental testing. It is NOT a prediction of
editing efficiency, allele selectivity or therapeutic benefit.
"""
import argparse, csv, itertools, json
from collections import defaultdict
from pathlib import Path

BASES = 'TCAG'
AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {''.join(c): a for c, a in zip(itertools.product(BASES, repeat=3), AA)}


def read_fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()


def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True)
    ap.add_argument('--cds-fasta', required=True)
    ap.add_argument('--variant', default='KRAS_G12D')
    ap.add_argument('--out', required=True)
    ap.add_argument('--top', type=int, default=6)
    a = ap.parse_args()
    run = Path(a.run)
    cds = read_fasta(a.cds_fasta)
    guides = [g for g in rows(run / 'discovery/02_guides.csv') if g['variant_id'] == a.variant]
    by = defaultdict(list)
    for r in rows(run / 'discovery/03_bystander_sites.csv'):
        if r['is_target'] != 'YES':
            by[r['guide_id']].append(r)
    off = {r['guide_id']: r for r in rows(run / 'offtarget/offtarget_gene_summary.csv')}
    chem = defaultdict(int)
    for r in rows(run / 'chemistry/constructs.csv'):
        chem[r['guide_id']] += 1
    cof = {r['guide_id']: r for r in rows(run / 'discovery/08_cofold.csv')}
    geo = {}
    gp = run / 'boltz_models/geometry.csv'
    if gp.exists():
        for r in rows(gp):
            geo.setdefault(r['guide_id'], []).append(r)

    out = []
    for g in guides:
        gid = g['guide_id']
        nonsyn, syn, noncoding = [], [], 0
        for b in by.get(gid, []):
            pos = int(b['cds_position'])
            if not 1 <= pos <= len(cds):
                noncoding += 1
                continue
            ci = (pos - 1) // 3
            codon = cds[ci * 3:ci * 3 + 3]
            if len(codon) < 3 or cds[pos - 1] != 'A':
                continue
            k = (pos - 1) % 3
            new = codon[:k] + 'G' + codon[k + 1:]
            label = f'c.{pos}A>G p.{CODON[codon]}{ci + 1}{CODON[new]}'
            (syn if CODON[codon] == CODON[new] else nonsyn).append(label)
        o = off.get(gid, {})
        n_off = int(o.get('off_target_genes', 0) or 0)
        out.append({
            'guide_id': gid, 'length': int(g['guide_length']), 'gc_percent': g['gc_percent'],
            'offtarget_non_KRAS_genes': n_off, 'offtarget_genes': o.get('genes', ''),
            'bystander_A_total': len(nonsyn) + len(syn),
            'bystander_nonsynonymous': len(nonsyn), 'bystander_synonymous': len(syn),
            'bystander_nonsyn_detail': '; '.join(nonsyn), 'bystander_syn_detail': '; '.join(syn),
            'published_chemistry_transferable_designs': chem.get(gid, 0),
            'cofold_mfe_kcal_mol': cof.get(gid, {}).get('mfe_kcal_mol', ''),
            'boltz_models': len(geo.get(gid, [])),
            'target_A_to_nearest_protein_A': min([float(x['target_A_to_nearest_protein_Ang']) for x in geo.get(gid, []) if x.get('target_A_to_nearest_protein_Ang')] or [float('nan')]),
        })
    out.sort(key=lambda r: (r['offtarget_non_KRAS_genes'], r['bystander_nonsynonymous'], r['bystander_A_total'],
                            -r['published_chemistry_transferable_designs'], r['length']))
    for i, r in enumerate(out, 1):
        r['triage_rank'] = i
    fields = ['triage_rank'] + [k for k in out[0] if k != 'triage_rank']
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out)
    for r in out[:a.top]:
        print({k: r[k] for k in ['triage_rank', 'guide_id', 'length', 'offtarget_non_KRAS_genes', 'bystander_nonsynonymous', 'bystander_synonymous', 'published_chemistry_transferable_designs']})


if __name__ == '__main__':
    main()
