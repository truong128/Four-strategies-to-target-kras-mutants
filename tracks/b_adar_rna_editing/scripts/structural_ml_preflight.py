#!/usr/bin/env python3
"""Independent, fail-closed structural/ML scientific readiness audit.

This program deliberately does not equate an existing 3D coordinate set with
validated phosphorothioate stereochemistry or an editing-efficiency model.
"""
import argparse, csv, json, shutil, hashlib
from pathlib import Path
from collections import Counter

def rows(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root', default='.')
    ap.add_argument('--guides', default='results/run/02_guides.csv')
    ap.add_argument('--transferred',default='results/run/restore2_published/transferred_per_residue.csv')
    ap.add_argument('--md-manifest', default='')
    ap.add_argument('--measured-labels', default='')
    ap.add_argument('--out', default='results/structural_ml_preflight')
    a=ap.parse_args()
    root=Path(a.root).resolve(); out=root/a.out;out.mkdir(parents=True,exist_ok=True)
    if not (root/a.guides).exists():
        raise SystemExit(f'ERROR: {root/a.guides} is missing. Run ./run_all.sh e2e first or supply --guides results/synthetic/02_guides.csv for offline smoke data.')
    guides=rows(root/a.guides); transferred=rows(root/a.transferred) if (root/a.transferred).exists() else []
    md={r['construct_id']:r for r in rows(root/a.md_manifest)} if a.md_manifest and (root/a.md_manifest).exists() else {}
    by={}
    for r in transferred:by.setdefault(r['transfer_id'],[]).append(r)
    guide_ids={g['guide_id'] for g in guides}; issues=[]; tasks=[]
    for cid,parts in by.items():
        parts.sort(key=lambda r:int(r['position_1based']))
        target=parts[0]['guide_id'];L=len(next(g['guide_rna_5to3'] for g in guides if g['guide_id']==target))
        n=len(parts); positions=[int(r['position_1based']) for r in parts]
        chemistry_valid=(n==L and positions==list(range(1,L+1)) and all(r.get('sugar') and r.get('link_to_next_3prime') for r in parts))
        ps=sum('PS' in r['link_to_next_3prime'].upper() for r in parts)
        task=dict(construct_id=cid,guide_id=target,chemical_position_manifest='PASS' if chemistry_valid else 'FAIL',ps_linkages=ps,stereo_control='REQUIRES_EXPLICIT_RP_SP_OR_DEFINED_STEREOENSEMBLE' if ps else 'NO_PS',parameterization='REQUIRES_VALIDATED_PARAMETERS',complex_template='REQUIRES_ADAR_ISOFORM_TEMPLATE_AND_ALIGNMENT',md='BLOCKED',editing_prediction='NOT_A_VALIDATED_PREDICTION')
        if cid in md:
            m=md[cid]
            criteria=['chemistry_verified','stereochemistry_verified','parameterization_validated','template_validated','equilibrated_explicit_solvent_ions']
            files=['structure_pdb','topology_prmtop','coordinates_inpcrd']
            ok=all(m.get(k,'').upper()=='YES' for k in criteria) and all(m.get(k) and Path(m[k]).exists() for k in files)
            task['md']='READY_FOR_INDEPENDENT_MD_EXECUTION' if ok and chemistry_valid else 'BLOCKED_INCOMPLETE_MD_EVIDENCE'
        tasks.append(task)
    labelstatus='NO_LABELED_DATA'; nlabels=0; groups=0; chem_cols=[]
    if a.measured_labels and (root/a.measured_labels).exists():
        labeled=rows(root/a.measured_labels);nlabels=len(labeled)
        if labeled:
            names=set(labeled[0]);chem_cols=[x for x in ('positionwise_chemistry','ribose_modification_pattern','backbone_pattern','chemical_construct_id') if x in names]
            groups=len({r.get('study_or_target_group','') for r in labeled if r.get('study_or_target_group')})
            if not {'guide_rna_5to3','target_rna_5to3','editing_fraction','study_or_target_group'}.issubset(names):labelstatus='BLOCKED_MISSING_MANDATORY_LABEL_COLUMNS'
            elif not chem_cols:labelstatus='BLOCKED_LABELS_NO_CHEMISTRY_FEATURES'
            elif groups<3:labelstatus='BLOCKED_INSUFFICIENT_INDEPENDENT_TARGET_GROUPS'
            else: labelstatus='DATA_PRESENT_NOT_INDEPENDENTLY_VALIDATED_REQUIRES_GROUP_SPLIT_AND_KRAS_EXTERNAL_SET'
    tools={x:bool(shutil.which(x)) for x in ('RNAfold','RNAcofold','tleap','antechamber','parmchk2','cpptraj','pmemd.cuda','sander','samtools','bcftools')}
    with (out/'structural_tasks.csv').open('w',newline='') as f:
        fields=['construct_id','guide_id','chemical_position_manifest','ps_linkages','stereo_control','parameterization','complex_template','md','editing_prediction'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(tasks)
    report={'guide_count':len(guides),'transferred_constructs':len(by),'constructs_with_complete_position_map':sum(x['chemical_position_manifest']=='PASS' for x in tasks),'md_ready_independent_evidence':sum(x['md']=='READY_FOR_INDEPENDENT_MD_EXECUTION' for x in tasks),'label_count':nlabels,'independent_target_groups':groups,'chemistry_columns_found':chem_cols,'editing_ml_status':labelstatus,'external_commands_available':tools,'scientific_blockers':['A noncovalent ADAR-RNA coordinate pose is not a validated catalytic transition-state model','Modified residue 3D topology/charges and Rp/Sp stereochemistry need independent chemical validation','RESTORE 2.0 guide-editing percentages cannot be extrapolated to KRAS without external KRAS validation','Mutant-versus-WT allele discrimination and off-target editing require experimental testing','MD conformational stability cannot establish cellular editing, delivery, or tumor killing'],'source_hashes':{str(x.relative_to(root)):sha(x) for x in [root/'reference/publisher/publisher_constructs.csv',root/'reference/publisher/publisher_residues.csv'] if x.exists()}}
    (out/'SCIENTIFIC_READINESS.json').write_text(json.dumps(report,indent=2))
    report_md=['# KRAS-ADAR scientific completion audit','',f'Guides: {len(guides)}; publisher chemistry-transferred constructs: {len(by)}; MD independently qualified: {report["md_ready_independent_evidence"]}.',f'Editing ML: {labelstatus}.','','## Modules that cannot be declared completed','']+['- '+x for x in report['scientific_blockers']]+['','## Minimum evidence to unblock','- Curated ADAR isoform-specific protein-RNA template and documented target sequence alignment.','- Atom-resolved modified guide chemistry, verified atom mapping, validated charges, torsions, and each PS stereochemical representation.','- Solvated/equilibrated Amber topology and coordinates with a documented quality-control report.','- Measured guide-specific editing data with chemistry features, target-level grouped validation, and a held-out KRAS experimental dataset.','- WT/mutant discrimination and transcriptome-wide RNA-editing measurement.']
    (out/'READINESS_REPORT.md').write_text('\n'.join(report_md)+'\n')
    print(json.dumps({'guide_count':len(guides),'chemistry_constructs':len(by),'qualified_md':report['md_ready_independent_evidence'],'ml':labelstatus,'report':str(out/'READINESS_REPORT.md')},indent=2))
if __name__=='__main__':main()
