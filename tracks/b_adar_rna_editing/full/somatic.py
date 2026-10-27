"""Executable matched tumor-normal calling (GATK Mutect2), VEP, VCF->KRAS HGVS panel.
Requires the user to provide correctly matched references, tumor/normal BAMs, and VEP cache.
No inference of cis/trans phase from unphased records.
"""
import argparse, csv, gzip, json, re
from pathlib import Path
from full.common import GateError, path_required, csvout, jsonout, run, which

def gatk(args):
    ref=path_required(args.reference,'--reference')
    tumor=path_required(args.tumor_bam,'--tumor-bam')
    normal=path_required(args.normal_bam,'--normal-bam')
    o=Path(args.out);o.mkdir(parents=True,exist_ok=True)
    gatk_exe=which('gatk') if args.execute else 'gatk'
    for bam in (tumor,normal):
        if not any(Path(str(bam)+suffix).exists() for suffix in ('.bai','.csi')) and not Path(str(bam).replace('.bam','.bai')).exists():
            raise GateError(f'BAM index missing for {bam}')
    if not Path(str(ref)+'.fai').exists():raise GateError('Reference FASTA .fai missing')
    vcf=o/'tumor_normal.unfiltered.vcf.gz'
    cmd=[gatk_exe,'Mutect2','-R',ref,'-I',tumor,'-I',normal,'-normal',args.normal_sample,'-O',vcf]
    # Modern GATK infers tumor sample from SM tags; -tumor is obsolete.
    if args.orientation_priors:cmd+=['--f1r2-tar-gz',o/'f1r2.tar.gz']
    if args.intervals:cmd+=['-L',path_required(args.intervals,'intervals')]
    if args.germline_resource:cmd+=['--germline-resource',path_required(args.germline_resource,'germline_resource')]
    if args.panel_of_normals:cmd+=['--panel-of-normals',path_required(args.panel_of_normals,'panel_of_normals')]
    run(cmd,o/'01_mutect2.log',args.execute)
    cmd2=[gatk_exe,'FilterMutectCalls','-R',ref,'-V',vcf,'-O',o/'tumor_normal.filtered.vcf.gz']
    if args.orientation_priors:
        orientation=[gatk_exe,'LearnReadOrientationModel','-I',o/'f1r2.tar.gz','-O',o/'orientation_priors.tar.gz']
        run(orientation,o/'orientation_model.log',args.execute)
        cmd2+=['--orientation-bias-artifact-priors',o/'orientation_priors.tar.gz']
    if args.contamination_table:
        cmd2+=['--contamination-table',path_required(args.contamination_table,'contamination table')]
    run(cmd2,o/'02_filter.log',args.execute)
    jsonout(o/'calling_scope.json',{'status':'EXECUTED' if args.execute else 'PLANNED','genome_assembly':args.assembly,'tumor_sample':args.tumor_sample,'normal_sample':args.normal_sample,'limitations':['Matched BAMs must have been aligned and QCed with the same genome build','Somatic variants should be supported by orthogonal QC, artifact/orientation filters, contamination and CN/LOH analyses','No automatic allele phase inference; the KRAS panel will not group different variants in cis without separate evidence','No purity/CNV/LOH correction or FFPE artifact clearance without appropriate samples and controls','Orientation model optional; contamination table is user supplied, not automatically estimated']})

def vep(args):
    inp=path_required(args.vcf,'--vcf');cache=Path(args.cache).expanduser().resolve()
    if not cache.is_dir():raise GateError('VEP cache directory missing')
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    exe=which('vep') if args.execute else 'vep'
    cmd=[exe,'--input_file',inp,'--output_file',out,'--offline','--cache','--dir_cache',cache,'--assembly',args.assembly,'--species','homo_sapiens','--vcf','--hgvs','--symbol','--canonical','--mane','--force_overwrite']
    if args.fasta:cmd+=['--fasta',path_required(args.fasta,'--fasta')]
    run(cmd,str(out)+'.log',args.execute)

def vcf_records(file):
    with (gzip.open(file,'rt') if str(file).endswith('.gz') else open(file)) as fh:
        csq=[]
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ,'):
                m=re.search(r'Format: ([^"<>]+)',line)
                if m:csq=m.group(1).strip().split('|')
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)<8:continue
            yield fields,csq

def extract_panel(args):
    panel=[]; excluded=[]
    for cells,headers in vcf_records(path_required(args.annotated_vcf,'--annotated-vcf')):
        chrom,pos,vid,ref,alts,qual,flt,info=cells[:8]
        if args.pass_only and flt!='PASS':continue
        if not headers:raise GateError('VEP CSQ Format header missing. Run VEP --vcf --hgvs first.')
        payload=next((x[4:] for x in info.split(';') if x.startswith('CSQ=')),None)
        if payload is None:continue
        for allele in payload.split(','):
            data=dict(zip(headers,allele.split('|')))
            if data.get('Allele') and data['Allele'] not in alts.split(','):continue
            if data.get('SYMBOL')!='KRAS' or not data.get('HGVSc'):continue
            tx=data.get('Feature','')
            if args.transcript and tx.split('.')[0]!=args.transcript.split('.')[0]:continue
            hgvs=data['HGVSc'].split(':')[-1]
            if re.fullmatch(r'c\.\d+[ACGT]>[ACGT]',hgvs):
                panel.append({'variant_id':f'KRAS_{chrom}_{pos}_{ref}_{alts}_{len(panel)+1}','hgvs_c':hgvs,'transcript':tx,'hgvs_p':data.get('HGVSp',''),'genomic_variant':f'{chrom}:{pos}:{ref}>{alts}','genome_assembly':args.assembly,'evidence':'VEP_ANNOTATION_NOT_PHASED','filter':flt,'variant_allele_frequency':'', 'phase':'UNDETERMINED'})
            else:excluded.append({'chrom':chrom,'pos':pos,'hgvs':hgvs,'reason':'Not single-nucleotide c.HGVS substitution'})
    o=Path(args.out);csvout(o,panel,['variant_id','hgvs_c','transcript','hgvs_p','genomic_variant','genome_assembly','evidence','filter','variant_allele_frequency','phase'])
    csvout(o.with_name(o.stem+'_excluded.csv'),excluded,['chrom','pos','hgvs','reason'])
    if not panel:raise GateError('No eligible KRAS SNV annotations. Verify transcript, VEP CSQ fields, assembly and filters.')

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('mutect2');a.add_argument('--orientation-priors',action='store_true');a.add_argument('--contamination-table');a.add_argument('--reference',required=True);a.add_argument('--tumor-bam',required=True);a.add_argument('--normal-bam',required=True);a.add_argument('--tumor-sample',required=True);a.add_argument('--normal-sample',required=True);a.add_argument('--intervals');a.add_argument('--panel-of-normals');a.add_argument('--germline-resource');a.add_argument('--out',default='results/somatic');a.add_argument('--assembly',default='GRCh38');a.add_argument('--execute',action='store_true')
    a=sub.add_parser('vep');a.add_argument('--vcf',required=True);a.add_argument('--cache',required=True);a.add_argument('--fasta');a.add_argument('--assembly',default='GRCh38');a.add_argument('--out',required=True);a.add_argument('--execute',action='store_true')
    a=sub.add_parser('panel');a.add_argument('--annotated-vcf',required=True);a.add_argument('--transcript',default='ENST00000256078');a.add_argument('--assembly',default='GRCh38');a.add_argument('--out',required=True);a.add_argument('--pass-only',action='store_true')
    a=p.parse_args();{'mutect2':gatk,'vep':vep,'panel':extract_panel}[a.cmd](a)
