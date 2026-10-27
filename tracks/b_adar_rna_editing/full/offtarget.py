"""Bowtie2 transcript-sense antisense-guide screening. Complementarity != editing.
Includes local/seed sites and indels. Enforce reverse-strand orientation of guide to sense RNAs.
"""
import argparse,collections,csv,re,subprocess
from pathlib import Path
from full.common import GateError, rows, csvout, jsonout, path_required, which, run

def load_gene_map(path):
    if not path:return {}
    with open(path_required(path,'--gene-map'),newline='',encoding='utf-8') as f:
        return {r['transcript_id']:r.get('gene_name') or r.get('gene_id','') for r in csv.DictReader(f,delimiter='\t')}

def build(args):
    fa=path_required(args.transcriptome,'--transcriptome');out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    run([which('bowtie2-build') if args.execute else 'bowtie2-build','--threads',str(max(1,args.threads)),fa,out],str(out)+'.index.log',args.execute)

def query(args):
    guide_rows=rows(path_required(args.guides,'--guides'));out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    fa=out/'guides.fa'
    with fa.open('w') as f:
        for g in guide_rows:f.write(f'>{g["guide_id"]}\n{g["guide_rna_5to3"].replace("U","T")}\n')
    idx=Path(args.index)
    cmd=[which('bowtie2') if args.execute else 'bowtie2','-x',str(idx),'-f','-U',str(fa),'--very-sensitive-local','--nofw','-k',str(args.max_alignments),'-p',str(max(1,args.threads)),'-S',str(out/'aligned.sam')]
    run(cmd,out/'bowtie2.log',args.execute)
    if not args.execute:return
    genes=load_gene_map(args.gene_map);hits=[]
    with (out/'aligned.sam').open() as f:
        for line in f:
            if line.startswith('@'):continue
            v=line.rstrip().split('\t')
            if len(v)<11:continue
            flag=int(v[1]);
            if flag&4 or not flag&16:continue
            name=v[0];ref=v[2];start=int(v[3]);cigar=v[5]
            optional={x.split(':',2)[0]:x.split(':',2)[-1] for x in v[11:] if x.count(':')>=2}
            nm=optional.get('NM','')
            # local alignment may clip guide; calculate matched query coverage from CIGAR
            ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
            aligned=sum(int(n) for n,op in ops if op in ('M','=','X','I'))
            coverage=aligned/len(v[9]) if v[9] else 0
            if coverage>=args.min_coverage:
                gene=genes.get(ref.split('.')[0],'') if genes else ''
                hits.append({'guide_id':name,'transcript_id':ref,'gene_name':gene,'is_on_target_gene':('YES' if gene.upper()==args.on_target_gene.upper() else 'NO') if genes else 'UNKNOWN','position_1based':start,'strand':'REVERSE_GUIDE_VS_SENSE_RNA','NM':nm,'query_coverage_fraction':round(coverage,3),'cigar':cigar,'mapq':v[4],'evidence':'COMPLEMENTARITY_NOT_EDITING'})
    csvout(out/'offtarget_alignments.csv',hits,['guide_id','transcript_id','gene_name','is_on_target_gene','position_1based','strand','NM','query_coverage_fraction','cigar','mapq','evidence'])
    per_gene=collections.defaultdict(set)
    for h in hits:
        if h['is_on_target_gene']!='YES':per_gene[h['guide_id']].add(h['gene_name'] or h['transcript_id'])
    summary=[{'guide_id':r['guide_id'],'off_target_genes':len(per_gene.get(r['guide_id'],())),'genes':';'.join(sorted(per_gene.get(r['guide_id'],())))[:500]} for r in guide_rows]
    csvout(out/'offtarget_gene_summary.csv',summary,['guide_id','off_target_genes','genes'])
    jsonout(out/'scope.json',{'number_of_alignments':len(hits),'on_target_gene':args.on_target_gene if genes else 'NOT_ANNOTATED','saturation_possible':True,'max_alignments_per_guide':args.max_alignments,'limits':['Local alignments can find partial matches; does not predict guide occupancy or ADAR editing','Transcriptome FASTA must be strand-normalized and contain relevant isoforms','Transcript expression, ADAR1 isoform, genomic SNPs and RNA editing require separate evidence','A mapping report cannot establish treatment safety']})
if __name__=='__main__':
    p=argparse.ArgumentParser();sp=p.add_subparsers(dest='cmd',required=True)
    a=sp.add_parser('index');a.add_argument('--transcriptome',required=True);a.add_argument('--out',required=True);a.add_argument('--threads',type=int,default=8);a.add_argument('--execute',action='store_true')
    a=sp.add_parser('screen');a.add_argument('--guides',required=True);a.add_argument('--index',required=True);a.add_argument('--out',required=True);a.add_argument('--min-coverage',type=float,default=0.7);a.add_argument('--max-alignments',type=int,default=1000);a.add_argument('--threads',type=int,default=8);a.add_argument('--gene-map');a.add_argument('--on-target-gene',default='KRAS');a.add_argument('--execute',action='store_true')
    a=p.parse_args();{'index':build,'screen':query}[a.cmd](a)
