"""RNA-seq pilot read counts at nominated RNA-editing sites.
A-to-G RNA difference is confounded by genomic allele G; in KRAS G12D
heterozygotes, G reads alone cannot prove editing. Reports raw fractions only.
"""
import argparse,collections,math,json
from pathlib import Path
from full.common import GateError,rows,csvout,jsonout,path_required

def count(bam,fa,contig,position,quality):
    import pysam
    with pysam.AlignmentFile(str(bam),'rb',reference_filename=str(fa)) as f:
        d=collections.Counter()
        for pile in f.pileup(contig,position-1,position,truncate=True,stepper='samtools',min_base_quality=quality,min_mapping_quality=20):
            if pile.reference_pos!=position-1:continue
            for read in pile.pileups:
                if read.is_del or read.is_refskip:continue
                al=read.alignment
                if al.is_duplicate or al.is_secondary or al.is_supplementary or al.is_qcfail:continue
                if read.query_position is None:continue
                base=al.query_sequence[read.query_position].upper()
                if base in 'ACGT':d[base]+=1
        return d

def run(args):
    try:import pysam
    except ImportError as e:raise GateError('pip install pysam') from e
    treatment=path_required(args.treated_bam,'treated BAM');control=path_required(args.control_bam,'control BAM');ref=path_required(args.reference_fasta,'reference FASTA')
    sites=rows(args.sites)
    out=[]
    for s in sites:
        chrom=s['contig'];pos=int(s['position_1based']);strand=s['transcript_strand']
        if strand not in '+-':raise GateError('transcript_strand must be + or -')
        if s['genomic_reference_base'].upper() not in 'ACGT':raise GateError('genomic_reference_base required')
        ctl=count(control,ref,chrom,pos,args.min_baseq);trt=count(treatment,ref,chrom,pos,args.min_baseq)
        a,g=('A','G') if strand=='+' else ('T','C')
        ct=ctl[a]+ctl[g];tt=trt[a]+trt[g]
        out.append({'site_id':s['site_id'],'contig':chrom,'position_1based':pos,'strand':strand,'control_A_or_T':ctl[a],'control_G_or_C':ctl[g],'treated_A_or_T':trt[a],'treated_G_or_C':trt[g],'control_G_like_fraction':ctl[g]/ct if ct else '', 'treated_G_like_fraction':trt[g]/tt if tt else '', 'min_covered':min(ct,tt),'status':'ALLELE_MIXTURE_NOT_EDITING_RATE_FOR_HETEROZYGOUS_KRAS_G12D','interpretation':'Requires matched DNA and phased reads / orthogonal mutant protein or ADAR perturbation evidence'})
    csvout(args.out,out,['site_id','contig','position_1based','strand','control_A_or_T','control_G_or_C','treated_A_or_T','treated_G_or_C','control_G_like_fraction','treated_G_like_fraction','min_covered','status','interpretation'])
    jsonout(str(args.out)+'.scope.json',{'warning':'At c.35G>A, edited mutant RNA (I read as G) and unedited WT RNA both produce G. Raw G counts cannot separate them. RNA allele proportions may also change due to expression and selection. Do not label G/(A+G) editing efficiency.','controls':['matched DNA VAF + RNA allelic expression','phasing via nearby heterozygous sites/long reads where possible','ADAR knockdown or catalytic inactive control','protein-level Gly12/Asp12 and KRAS-GTP/pERK assays','full RNAseq differential editing including mapping bias and SNP filters']})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--treated-bam',required=True);p.add_argument('--control-bam',required=True);p.add_argument('--reference-fasta',required=True);p.add_argument('--sites',required=True);p.add_argument('--out',required=True);p.add_argument('--min-baseq',type=int,default=25);run(p.parse_args())
