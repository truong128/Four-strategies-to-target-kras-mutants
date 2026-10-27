#!/usr/bin/env python3
"""Pull public ADAR/KRAS evidence from the cBioPortal API (CCLE 2025, TCGA-PAAD GDC, CPTAC-PDAC).

DepMap's own download endpoint now requires a browser captcha; cBioPortal serves the same
CCLE/DepMap release (ccle_broad_2025) without authentication.
"""
import argparse, csv, json, sys, time, urllib.error, urllib.request
from pathlib import Path

API = 'https://www.cbioportal.org/api'
GENES = {'ADAR': 103, 'ADARB1': 104, 'ADARB2': 105, 'KRAS': 3845, 'IFIH1': 64135, 'STAT1': 6772,
         'ISG15': 9636, 'IFIT1': 3434, 'MX1': 4599, 'OAS1': 4938, 'DDX58': 23586, 'ZBP1': 81030,
         'DHX9': 1660, 'GAPDH': 2597}


def call(path, body=None):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(API + path, data=data, headers={'Content-Type': 'application/json', 'Accept': 'application/json', 'User-Agent': 'curl/8.5.0'})
    for attempt in range(8):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code != 429 or attempt == 7:
                raise
            time.sleep(15 * (attempt + 1))


def expression(profile, study, sample_list=None):
    sl = sample_list or f'{study}_all'
    out = call(f'/molecular-profiles/{profile}/molecular-data/fetch?projection=SUMMARY',
               {'entrezGeneIds': list(GENES.values()), 'sampleListId': sl})
    rev = {v: k for k, v in GENES.items()}
    table = {}
    for x in out:
        table.setdefault(x['sampleId'], {})[rev[x['entrezGeneId']]] = x['value']
    return table


def kras_mutations(profile, study, sample_list=None):
    sl = sample_list or f'{study}_all'
    out = call(f'/molecular-profiles/{profile}/mutations/fetch?projection=SUMMARY',
               {'entrezGeneIds': [3845], 'sampleListId': sl})
    res = {}
    for m in out:
        res.setdefault(m['sampleId'], []).append(m.get('proteinChange', ''))
        DETAILS.append({'study': study, 'sample': m['sampleId'], 'proteinChange': m.get('proteinChange', ''),
                        'startPosition': m.get('startPosition'), 'refAllele': m.get('referenceAllele'), 'varAllele': m.get('variantAllele'),
                        't_alt_count': m.get('tumorAltCount'), 't_ref_count': m.get('tumorRefCount')})
    return res


DETAILS = []


def write(path, rows, fields):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    study = 'ccle_broad_2025'
    expr = expression(f'{study}_rna_seq_mrna', study)
    mut = kras_mutations(f'{study}_mutations', study)
    rows = []
    for s, d in expr.items():
        lineage = s.split('_', 1)[1] if '_' in s else ''
        row = {'sample': s, 'lineage': lineage, 'KRAS_mutations': ';'.join(sorted(set(mut.get(s, [])))) or 'none_called'}
        row.update({g: d.get(g, '') for g in GENES})
        rows.append(row)
    write(out / 'ccle_expression_kras.csv', rows, ['sample', 'lineage', 'KRAS_mutations'] + list(GENES))

    for tag, prof, st, ml in [('tcga_paad_gdc', 'paad_tcga_gdc_mrna_seq_tpm', 'paad_tcga_gdc', 'paad_tcga_gdc_mutations'),
                              ('cptac_pdac_mrna', 'paad_cptac_2021_mrna', 'paad_cptac_2021', 'paad_cptac_2021_mutations')]:
        try:
            e = expression(prof, st)
            m = kras_mutations(ml, st)
        except Exception as ex:
            print('FAILED', tag, ex, file=sys.stderr)
            continue
        rows = []
        for s, d in e.items():
            row = {'sample': s, 'KRAS_mutations': ';'.join(sorted(set(m.get(s, [])))) or 'none_called'}
            row.update({g: d.get(g, '') for g in GENES})
            rows.append(row)
        write(out / f'{tag}_expression_kras.csv', rows, ['sample', 'KRAS_mutations'] + list(GENES))
        print(tag, len(rows), 'samples')

    prot = expression('paad_cptac_2021_protein_quantification', 'paad_cptac_2021')
    cell = {}
    for ctype in ('SAMPLE', 'PATIENT'):
        for x in call(f'/studies/paad_cptac_2021/clinical-data?clinicalDataType={ctype}&attributeId=NEOPLASTIC_CELLULARITY&projection=SUMMARY'):
            cell[x['patientId']] = x['value']
    rows = [{'sample': s, 'neoplastic_cellularity': cell.get(s, ''), **{g: d.get(g, '') for g in GENES}} for s, d in prot.items()]
    write(out / 'cptac_pdac_protein.csv', rows, ['sample', 'neoplastic_cellularity'] + list(GENES))
    write(out / 'kras_mutation_details.csv', DETAILS, ['study', 'sample', 'proteinChange', 'startPosition', 'refAllele', 'varAllele', 't_alt_count', 't_ref_count'])
    print('protein', len(rows), 'samples; KRAS mutation records', len(DETAILS))


if __name__ == '__main__':
    main()
