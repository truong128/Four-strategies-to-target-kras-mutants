"""Pinned public reference download utilities with recorded provenance.
Use consistent GENCODE GRCh38 release and Ensembl transcript versions.
"""
import argparse,gzip,hashlib,re,shutil,urllib.request
from pathlib import Path
from full.common import GateError,jsonout,digest

def download(url,destination):
    dest=Path(destination);dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=Path(str(dest)+'.partial')
    with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'KRAS-ADAR-code-study/1.0'}),timeout=60) as r,tmp.open('wb') as o:
        shutil.copyfileobj(r,o)
    if tmp.stat().st_size<100:raise GateError('Download too small')
    tmp.replace(dest);return dest

def gencode(args):
    rel=int(args.release)
    if rel<20 or rel>100:raise GateError('Human GENCODE release must be a numbered GRCh38 release, not GRCh37')
    name=f'gencode.v{rel}.transcripts.fa.gz'
    url=f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{rel}/{name}'
    dest=Path(args.out)/name
    if args.execute:
        if not dest.exists() or args.replace:download(url,dest)
        outfile=dest.with_suffix('')
        if not outfile.exists() or args.replace:
            with gzip.open(dest,'rb') as a,outfile.open('wb') as b:shutil.copyfileobj(a,b)
    else:outfile=dest.with_suffix('')
    jsonout(Path(args.out)/'transcriptome_reference_manifest.json',{'release':rel,'assembly':'GRCh38','source':url,'file':str(outfile),'downloaded':bool(args.execute),'sha256_gzip':digest(dest) if dest.exists() else None,'sha256_fasta':digest(outfile) if outfile.exists() else None,'note':'Transcriptome is sense-oriented transcript sequences; download pinned GTF/genome from the matching release for mapping analyses'})

def protein(args):
    accession=args.accession
    if not re.fullmatch(r'[A-Z0-9-]{6,20}',accession):raise GateError('Supply exact verified UniProt protein/isoform accession')
    url=f'https://rest.uniprot.org/uniprotkb/{accession}.fasta'
    dest=Path(args.out)/f'{accession}.fasta'
    if args.execute:download(url,dest)
    jsonout(Path(args.out)/f'{accession}_provenance.json',{'source':url,'accession':accession,'file':str(dest),'downloaded':bool(args.execute),'sha256':digest(dest) if dest.exists() else None,'mandatory':'Verify the isoform is ADAR1 p110/p150 (not ADAR2), localization, transcript expression and domain identity before structure prediction'})

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='cmd',required=True)
    x=sub.add_parser('gencode');x.add_argument('--release',type=int,required=True);x.add_argument('--out',required=True);x.add_argument('--execute',action='store_true');x.add_argument('--replace',action='store_true')
    x=sub.add_parser('protein');x.add_argument('--accession',required=True);x.add_argument('--out',required=True);x.add_argument('--execute',action='store_true')
    args=p.parse_args();{'gencode':gencode,'protein':protein}[args.cmd](args)
