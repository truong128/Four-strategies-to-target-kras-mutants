#!/usr/bin/env python3
"""Fail-fast input validity, sequence coordinate and sample-provenance checks."""
import argparse,csv,json,re
from pathlib import Path

def run(args):
 cds=''.join(x.strip() for x in Path(args.cds).read_text().splitlines() if not x.startswith('>')).upper().replace('U','T')
 if not re.fullmatch('[ACGT]+',cds) or len(cds)%3 or cds[:3]!='ATG':raise ValueError('Invalid CDS: must be exact coding DNA with ATG and triplet length')
 if cds[33:36] not in ['GGT','GGC','GGA','GGG']:raise ValueError(f'KRAS codon12 wrong: {cds[33:36]}')
 with open(args.variants,newline='') as f:vs=list(csv.DictReader(f))
 errors=[]
 for v in vs:
  for m in v['hgvs_c'].split(';'):
   m=m.strip();match=re.fullmatch(r'c\.(\d+)([ACGT])>([ACGT])',m)
   if not match:errors.append(f'{v["variant_id"]}:unsupported {m}');continue
   pos,old,new=int(match[1]),match[2],match[3]
   if not 1<=pos<=len(cds) or cds[pos-1]!=old:errors.append(f'{v["variant_id"]}:reference mismatch {m}')
 if errors:raise ValueError('; '.join(errors))
 print(json.dumps({'CDS_length':len(cds),'variants_checked':len(vs),'KRAS_codon12':cds[33:36],'reference_checks':'PASS','warning':'Not a confirmation of tumor-specific somatic genotype or variant phasing'},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cds',required=True);p.add_argument('--variants',default='config/variants.csv');run(p.parse_args())
