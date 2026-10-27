"""Matched tumor/normal FASTQ -> QCed, RG-tagged, duplicate-marked BAM.
User controls actual sample IDs, reference, known sites and genome build.
"""
import argparse,subprocess,json,shutil
from pathlib import Path
from full.common import GateError,path_required,jsonout,which,run

def align(args):
    ref=path_required(args.reference,'--reference');r1=path_required(args.read1,'--read1');r2=path_required(args.read2,'--read2')
    if not Path(str(ref)+'.bwt.2bit.64').exists() and args.execute:
        raise GateError('bwa-mem2 indexed reference not detected: run bwa-mem2 index REF')
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    rg=f'@RG\\tID:{args.sample}.1\\tSM:{args.sample}\\tPL:ILLUMINA\\tLB:{args.sample}.library'
    pre=out/'coordinate_sorted.bam';md=out/'duplicate_marked.bam'
    commands={
        'align':['bwa-mem2','mem','-t',str(args.threads),'-R',rg,str(ref),str(r1),str(r2)],
        'sort':['samtools','sort','-@',str(args.threads),'-o',str(pre),'-'],
        'duplicates':['gatk','MarkDuplicates','-I',str(pre),'-O',str(md),'-M',str(out/'dup_metrics.txt'),'--CREATE_INDEX','true'],
        'flagstat':['samtools','flagstat',str(md)]
    }
    jsonout(out/'alignment_commands.json',commands)
    if not args.execute:return
    for exe in ('bwa-mem2','samtools','gatk'):which(exe)
    with (out/'alignment.log').open('w') as log:
        a=subprocess.Popen(commands['align'],stdout=subprocess.PIPE,stderr=log)
        b=subprocess.Popen(commands['sort'],stdin=a.stdout,stdout=log,stderr=log)
        a.stdout.close()
        bcode=b.wait();acode=a.wait()
        if bcode or acode:raise GateError(f'Alignment / sort failed: bwa={acode}, sort={bcode}')
    run(commands['duplicates'],out/'duplicates.log',True)
    with (out/'flagstat.txt').open('w') as fh:
        subprocess.run(commands['flagstat'],stdout=fh,check=True)
    if args.known_sites:
        known=path_required(args.known_sites,'known-sites VCF')
        recal=out/'recal.table'
        run(['gatk','BaseRecalibrator','-R',str(ref),'-I',str(md),'--known-sites',str(known),'-O',str(recal)],out/'bqsr.log',True)
        run(['gatk','ApplyBQSR','-R',str(ref),'-I',str(md),'--bqsr-recal-file',str(recal),'-O',str(out/'recalibrated.bam')],out/'apply_bqsr.log',True)
        run(['samtools','index',str(out/'recalibrated.bam')],out/'recalibrated_index.log',True)
    jsonout(out/'alignment_status.json',{'status':'EXECUTED_NEEDS_MANUAL_QC','sample':args.sample,'reference':str(ref),'dedup_bam':str(md),'bqsr_bam':str(out/'recalibrated.bam') if args.known_sites else None,'warning':'Assume matched patient samples, same genome build, appropriate library orientation and sequencing quality.'})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--read1',required=True);p.add_argument('--read2',required=True);p.add_argument('--sample',required=True);p.add_argument('--out',required=True);p.add_argument('--known-sites');p.add_argument('--threads',type=int,default=8);p.add_argument('--execute',action='store_true')
    align(p.parse_args())
