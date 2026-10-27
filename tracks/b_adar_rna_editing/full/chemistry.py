"""Canonicalizes independently source-checked RESTORE2 per-residue transfer tables.
Does not infer PS Rp/Sp chirality from a stereo-random published oligonucleotide.
"""
import argparse,collections
from pathlib import Path
from full.common import rows,csvout,jsonout,GateError

def normalize(args):
    guides={g['guide_id']:g for g in rows(args.guides)}
    raw=rows(args.transferred);groups=collections.defaultdict(list)
    for r in raw:groups[r['transfer_id']].append(r)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    records=[];positions=[];problems=[]
    for tid,ss in groups.items():
        try:
            ss=sorted(ss,key=lambda x:int(x['position_1based']))
            gid=ss[0]['guide_id'];g=guides[gid];L=int(g['guide_length'])
            if len(ss)!=L or [int(t['position_1based']) for t in ss]!=list(range(1,L+1)):
                raise GateError('Noncontiguous chemistry positions')
            s=''.join(t['kras_base'] for t in ss)
            if s.replace('T','U')!=g['guide_rna_5to3'].replace('T','U'):
                raise GateError('KRAS guide sequence mismatch')
            for r in ss:
                link=r['link_to_next_3prime'].strip().upper()
                sugar=r['sugar'].strip().upper()
                if 'SPECIAL' in sugar or 'TERMINAL' in link or 'OTHER' in sugar:raise GateError('Special chemistry requires validated custom parser/parameters')
                if int(r['position_1based'])==L and link not in ('NONE_3PRIME','NONE'):
                    raise GateError('3prime terminal nucleotide must not have a linkage to a nonexistent next residue')
                if not link or (link!='PO' and not link.startswith('PS') and link not in ('NONE_3PRIME','NONE')):
                    raise GateError('Unrecognized backbone linkage')
            cid=tid;ps=sum(r['link_to_next_3prime'].upper().startswith('PS') for r in ss)
            records.append({'construct_id':cid,'guide_id':gid,'variant_id':g['variant_id'],'guide_rna_5to3':g['guide_rna_5to3'],'target_rna_5to3':g['target_rna_5to3'],'guide_index_0based':g['guide_index_0based'],'source':ss[0]['template_source_key'],'chemistry_published_on_original_target':'YES','chemistry_transfer_to_KRAS_unvalidated':'YES','PS_stereochemistry_unresolved':'YES' if ps else 'NO','PS_linkages':ps,'MD_ready':'NO_REQUIRES_CUSTOM_VALIDATED_PARAMETERS'})
            for r in ss:
                positions.append({'construct_id':cid,'position_1based':r['position_1based'],'base':r['kras_base'],'sugar':r['sugar'],'link_to_next_3prime':r['link_to_next_3prime'],'stereo':'UNSPECIFIED_RP_SP' if r['link_to_next_3prime'].strip().upper().startswith('PS') else 'NOT_APPLICABLE','origin':r['template_source_key']})
        except Exception as e:problems.append({'transfer_id':tid,'reason':str(e)})
    csvout(out/'constructs.csv',records,['construct_id','guide_id','variant_id','guide_rna_5to3','target_rna_5to3','guide_index_0based','source','chemistry_published_on_original_target','chemistry_transfer_to_KRAS_unvalidated','PS_stereochemistry_unresolved','PS_linkages','MD_ready'])
    csvout(out/'positions.csv',positions,['construct_id','position_1based','base','sugar','link_to_next_3prime','stereo','origin'])
    csvout(out/'excluded.csv',problems,['transfer_id','reason'])
    jsonout(out/'summary.json',{'constructs':len(records),'excluded':len(problems),'chemistry_origin':'PUBLISHER_SEQUENCE_PATTERN_TRANSFER','validation':'NOT_VALIDATED_IN_KRAS','atomic_MD':'REQUIRES_3D_TOPOLOGY_AND_RESIDUE_PARAMETERS'})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--transferred',required=True);p.add_argument('--guides',required=True);p.add_argument('--out',required=True);normalize(p.parse_args())
