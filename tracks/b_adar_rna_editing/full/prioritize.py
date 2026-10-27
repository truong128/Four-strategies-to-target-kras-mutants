"""Engineering triage plus *conditional* measured-evidence ranking.
No made-up treatment efficacy scores; rank only with measured KRAS edit and functional data.
"""
import argparse,collections,json,math
from pathlib import Path
from full.common import GateError, rows, csvout, jsonout

def priority(args):
    guides=rows(args.guides);bystanders=rows(args.bystanders)
    bct=collections.Counter(r['guide_id'] for r in bystanders if r['is_target']=='NO')
    experimental={r['guide_id']:r for r in rows(args.experimental)} if args.experimental else {}
    off=collections.Counter();off_genes=collections.defaultdict(set);on_target=collections.Counter()
    if args.offtarget_hits:
        for r in rows(args.offtarget_hits):
            if r.get('is_on_target_gene')=='YES':on_target[r['guide_id']]+=1;continue
            off[r['guide_id']]+=1
            off_genes[r['guide_id']].add(r.get('gene_name') or r['transcript_id'])
    triage=[];evaluable=[]
    for g in guides:
        gid=g['guide_id'];a=bct[gid]
        hit=off.get(gid,0) if args.offtarget_hits else 'NOT_SCREENED'
        hit_genes=len(off_genes.get(gid,())) if args.offtarget_hits else 'NOT_SCREENED'
        on_hits=on_target.get(gid,0) if args.offtarget_hits else 'NOT_SCREENED'
        # The G12D mismatch issue is explicit: WT G:C pair can have greater stability than mutant A:C.
        wt_base=g['wt_guide_target_base'];mut_base=g['mutant_guide_target_base']
        if mut_base=='A' and wt_base=='G':allele='WT_G_C_PERFECT_PAIR_MUTANT_A_C_MISMATCH__BINDING_SPECIFICITY_NOT_ASSUMED'
        else:allele='ALLELE_DISCRIMINATION_NOT_VALIDATED'
        triage.append({'guide_id':gid,'variant_id':g['variant_id'],'bystander_A_count_not_editing_rate':a,'complementary_transcript_hits_not_editing':hit,'complementary_off_target_genes':hit_genes,'on_target_gene_transcript_hits':on_hits,'WT_caveat':allele,'therapeutic_evidence':'NONE' if gid not in experimental else 'USER_PROVIDED_CHECK_INPUT'})
        if gid in experimental:
            v=experimental[gid]
            req=['mutant_edit_fraction','wt_bystander_fraction','max_offtarget_edit_fraction','KRAS_GTP_reduction_fraction','cell_viability_effect_fraction','assay_batch','biological_replicates']
            if any(not v.get(k) for k in req):continue
            try:
                z=[float(v[k]) for k in req[:5]]
                n=int(v['biological_replicates'])
            except ValueError:continue
            if n<args.min_replicates or any(not 0<=t<=1 for t in z):continue
            if z[0]<args.min_edit or z[1]>args.max_wt_bystander or z[2]>args.max_offtarget or z[3]<args.min_KRAS_GTP_reduction:continue
            # Ranking uses measured normalized on-target and KRAS biochemical effect, with safety gate;
            # do not turn structural MD/model prediction into a clinical safety score.
            e=(z[0]+z[3])/2
            evaluable.append({'guide_id':gid,'measured_technical_score':round(e,4),'mutant_edit_fraction':z[0],'KRAS_GTP_reduction_fraction':z[3],'cell_viability_effect_fraction':z[4],'biological_replicates':n,'status':'EXPERIMENTALLY_GATED_IN_VITRO_RESEARCH_PRIORITY_ONLY'})
    if evaluable:evaluable.sort(key=lambda x:-x['measured_technical_score'])
    o=Path(args.out);o.mkdir(parents=True,exist_ok=True)
    csvout(o/'engineering_triage.csv',triage,list(triage[0]) if triage else ['guide_id'])
    csvout(o/'conditional_measured_ranking.csv',evaluable,['guide_id','measured_technical_score','mutant_edit_fraction','KRAS_GTP_reduction_fraction','cell_viability_effect_fraction','biological_replicates','status'])
    jsonout(o/'ranking_status.json',{'total_guides':len(guides),'with_complete_measured_evidence_passing_gates':len(evaluable),'status':'IN_VITRO_MEASURED_TECHNICAL_RANKING' if evaluable else 'NO_THERAPEUTIC_RANKING_ABSTAIN','important':['Assay results must be quality-controlled and biological replicates independent','Not a model of delivery, in vivo toxicity, tumor response, PK or clinical efficacy','Bystander/off-target thresholds are research design defaults, not clinical safe thresholds']})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--guides',required=True);p.add_argument('--bystanders',required=True);p.add_argument('--offtarget-hits');p.add_argument('--experimental');p.add_argument('--out',required=True)
    p.add_argument('--min-replicates',type=int,default=3);p.add_argument('--min-edit',type=float,default=0.2);p.add_argument('--max-wt-bystander',type=float,default=0.05);p.add_argument('--max-offtarget',type=float,default=0.05);p.add_argument('--min-KRAS-GTP-reduction',type=float,default=0.2)
    priority(p.parse_args())
