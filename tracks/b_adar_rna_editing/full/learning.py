"""Chemistry-aware ML training + grouped cross-target validation.
Trains only from independently curated assay labels; reports OOD KRAS status.
"""
import argparse, collections, json, math
from pathlib import Path
from full.common import GateError, rows, csvout, jsonout, path_required

SUGARS=('RNA_2OH','2OME','2F','DNA_2H','LNA','OTHER')
BONDS=('PO','PS','OTHER')

def normalized_chem(row):
    v=row.get('sugar','').upper().replace("'",'').replace('-','').replace('_','').replace(' ','')
    if v in ('RNA2OH','2OH','RNA','OH'):return 'RNA_2OH'
    if v in ('2OME','OM E','2OMETHYL'):return '2OME'
    if v in ('2F','2FLUORO'):return '2F'
    if v in ('DNA','DNA2H','2H'):return 'DNA_2H'
    if v.startswith('LNA'):return 'LNA'
    return 'OTHER'

def transform(seq,target,site,chemical_rows):
    seq=seq.upper().replace('T','U');target=target.upper().replace('T','U');length=len(seq)
    if length!=len(target):raise GateError('guide and target window lengths differ')
    if not 0<=site<length:raise GateError('Orphan position invalid')
    # Relative positions in GUIDE 5->3; no spurious chemical inference from sequence.
    ret=[length,site/length,(seq.count('G')+seq.count('C'))/length,
         (target.count('A')/length),target[(length-1)-site]=='A']
    ret=[float(x) for x in ret]
    for b in 'ACGU':ret.append(seq.count(b)/length)
    for b in 'ACGU':ret.append(target.count(b)/length)
    for p in (0,1,2,3,-1,-2,-3,-4):
        t=(length-1-site)+p
        for b in 'ACGU':ret.append(float(0<=t<length and target[t]==b))
    chem={int(r['position_1based'])-1:r for r in chemical_rows}
    if len(chem)!=len(chemical_rows):raise GateError('Duplicate chemistry position rows')
    if len(chem)!=length or set(chem)!=set(range(length)):raise GateError('Chemistry manifest must include every guide nucleotide exactly once')
    for i in range(length):
        base=chem[i].get('base','').upper().replace('T','U')
        if base and base not in ('A','C','G','U'):raise GateError(f'Noncanonical base at chemistry position {i+1}')
        if base and base!=seq[i]:raise GateError(f'Chemistry base/guide sequence mismatch at {i+1}')
    mods=[normalized_chem(chem[i]) for i in range(length)]
    if 'OTHER' in mods:raise GateError('Unrecognized chemical modification; cannot infer standardized ML feature')
    terminal=chem[length-1].get('link_to_next_3prime','').upper()
    if terminal not in ('NONE_3PRIME','NONE','', 'NA'):
        raise GateError('Terminal nucleotide cannot have an internucleotide 3-prime linkage')
    bond=[]
    for i in range(length):
        if i == length-1:
            bond.append(None)  # no physical 3-prime internucleotide bond at terminal residue
            continue
        raw=chem[i].get('link_to_next_3prime','').strip().upper()
        if raw == 'PO':bond.append('PO')
        elif raw in ('PS','PS_RP','PS_SP','PS_UNSPECIFIED'):bond.append('PS')
        else:raise GateError(f'Unknown or missing internal linkage at position {i+1}: {raw!r}')
    for v in SUGARS:ret.append(sum(m==v for m in mods)/length)
    for v in BONDS:ret.append(sum(b==v for b in bond)/max(1,length-1))
    # Distances around orphan C; directional on guide.
    for delta in (-5,-4,-3,-2,-1,0,1,2,3,4,5):
        i=site+delta
        for v in SUGARS:ret.append(float(0<=i<length and mods[i]==v))
        for v in BONDS:ret.append(float(0<=i<length and bond[i]==v))
    return ret

def grouped_library(filename, key_name):
    chem=collections.defaultdict(list)
    for x in rows(filename):chem[x[key_name]].append(x)
    return chem

def dataset(records,manifest,require_labels=True):
    examples=[];excluded=[]
    for r in records:
        try:
            cid=r['construct_id'];seq=r['guide_rna_5to3'];target=r['target_rna_5to3'];site=int(r['guide_index_0based'])
            if cid not in manifest:raise GateError('Missing per-position chemistry')
            x=transform(seq,target,site,manifest[cid])
            label=float(r['editing_fraction']) if require_labels else None
            if require_labels and not 0<=label<=1:raise GateError('Editing fraction outside 0..1')
            examples.append((r,x,label))
        except Exception as e:
            excluded.append({'construct_id':r.get('construct_id',''), 'reason':str(e)})
    return examples,excluded

def train(args):
    try:
        import numpy as np
        import joblib
        from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
        from sklearn.model_selection import GroupKFold
        from sklearn.metrics import mean_absolute_error, mean_squared_error
    except ImportError as e:raise GateError('Install numpy scikit-learn joblib for ML training') from e
    raw=rows(path_required(args.labels,'--labels'));chem=grouped_library(path_required(args.chemistry,'--chemistry'),'construct_id')
    required={'construct_id','guide_rna_5to3','target_rna_5to3','guide_index_0based','target_group','assay_context','editing_fraction'}
    if any(not required.issubset(r) for r in raw):raise GateError('Required columns: '+','.join(sorted(required)))
    data,excluded=dataset(raw,chem)
    o=Path(args.out);o.mkdir(parents=True,exist_ok=True)
    csvout(o/'excluded.csv',excluded,['construct_id','reason'])
    groups=[r['target_group'] for r,x,y in data]
    for r,x,y in data:
        if not r.get('target_group','').strip():raise GateError('target_group required for group-held-out validation')
    sequence_groups=collections.defaultdict(set)
    for r,x,y in data:
        sequence_groups[r['target_rna_5to3'].upper().replace('T','U')].add(r['target_group'])
    if any(len(g)>1 for g in sequence_groups.values()):
        raise GateError('DATA LEAKAGE: an identical RNA target sequence is assigned multiple target_group labels')
    if len(data)<args.min_samples or len(set(groups))<args.min_groups:
        jsonout(o/'ml_report.json',{'status':'NOT_TRAINABLE_INSUFFICIENT_DATA','samples':len(data),'independent_targets':len(set(groups)),'min_samples':args.min_samples,'min_groups':args.min_groups})
        return
    X=np.array([x for r,x,y in data],dtype=float);Y=np.array([y for r,x,y in data]);groups=np.array(groups)
    cv=GroupKFold(n_splits=min(args.cv_folds,len(set(groups))))
    models={'RandomForest':RandomForestRegressor(n_estimators=250,min_samples_leaf=3,random_state=7,n_jobs=-1),
            'ExtraTrees':ExtraTreesRegressor(n_estimators=300,min_samples_leaf=3,random_state=7,n_jobs=-1)}
    summaries=[];best=None;preds_all=[]
    for model_name,model in models.items():
        preds=np.zeros(len(Y));baseline=np.zeros(len(Y))
        for trainidx,testidx in cv.split(X,Y,groups):
            model.fit(X[trainidx],Y[trainidx]);preds[testidx]=np.clip(model.predict(X[testidx]),0,1)
            baseline[testidx]=Y[trainidx].mean()
        mae=float(mean_absolute_error(Y,preds));rmae=float(mean_absolute_error(Y,baseline))
        sums={'model':model_name,'group_CV_MAE':mae,'group_CV_RMSE':float(np.sqrt(mean_squared_error(Y,preds))),'fold_matched_baseline_MAE':rmae,'better_than_baseline':mae<rmae,'target_groups':len(set(groups)),'samples':len(Y)}
        summaries.append(sums)
        preds_all.extend({'model':model_name,'construct_id':data[i][0]['construct_id'],'target_group':groups[i],'observed_editing_fraction':Y[i],'cv_predicted_fraction':preds[i]} for i in range(len(Y)))
        if mae<rmae and (best is None or mae<best[0]):best=(mae,model_name)
    csvout(o/'cross_target_cv_predictions.csv',preds_all,list(preds_all[0]))
    if not best:
        jsonout(o/'ml_report.json',{'status':'ABSTAIN_NO_MODEL_BEATS_GROUPED_BASELINE','results':summaries})
        return
    chosen=models[best[1]];chosen.fit(X,Y)
    joblib.dump(chosen,o/'chemistry_model.joblib')
    from full.common import digest
    jsonout(o/'ml_provenance.json',{'labels_sha256':digest(args.labels),'chemistry_sha256':digest(args.chemistry),'feature_schema':'site_relative_sugar_backbone_v1','target_groups':sorted(set(groups)),'training_data_assay_specific':'YES','KRAS_independent_external_validation':'NOT_PERFORMED','training_mean_editing_fraction':float(Y.mean())})
    jsonout(o/'ml_report.json',{'status':'EXPLORATORY_CROSS_TARGET_MODEL_NOT_EXTERNALLY_VALIDATED_ON_KRAS','chosen':best[1],'results':summaries,'feature_count':X.shape[1],'chemistry_model_file':'chemistry_model.joblib','external_KRAS_validation':'REQUIRED','warning':'Predictions are not measured editing yields, not efficacy claims; group CV only on available literature targets'})

def infer(args):
    try:import joblib,numpy as np
    except ImportError as e:raise GateError('Install joblib numpy') from e
    model=joblib.load(path_required(args.model,'--model'));chem=grouped_library(path_required(args.chemistry,'--chemistry'),'construct_id')
    data,bad=dataset(rows(args.candidates),chem,False)
    out=[]
    for r,x,_ in data:
        pred=float(np.clip(model.predict(np.asarray([x]))[0],0,1))
        out.append({'construct_id':r['construct_id'],'guide_id':r.get('guide_id',''),'predicted_editing_fraction_exploratory':round(pred,4),'status':'OUT_OF_DOMAIN_NOT_VALIDATED_ON_KRAS','not_therapeutic_efficiency':'YES'})
    csvout(args.out,out,['construct_id','guide_id','predicted_editing_fraction_exploratory','status','not_therapeutic_efficiency'])
    csvout(str(args.out)+'.excluded.csv',bad,['construct_id','reason'])


def external_validate(args):
    """Evaluate frozen fitted model on truly independent target groups, without refitting."""
    try:import joblib,numpy as np
    except ImportError as e:raise GateError('Install joblib and numpy') from e
    from sklearn.metrics import mean_absolute_error, mean_squared_error
    meta=json.loads(path_required(args.provenance,'--provenance').read_text())
    train_groups=set(meta.get('target_groups',[]))
    labels=rows(path_required(args.labels,'--labels'))
    if not labels:raise GateError('External validation labels are empty')
    holdout_groups={r.get('target_group','') for r in labels}
    if '' in holdout_groups:raise GateError('External labels must declare target_group')
    overlap=holdout_groups&train_groups
    if overlap:raise GateError('Target group leakage into external test data: '+','.join(sorted(overlap)))
    chem=grouped_library(path_required(args.chemistry,'--chemistry'),'construct_id')
    examples,excluded=dataset(labels,chem,True)
    if not examples:raise GateError('No externally validatable examples')
    model=joblib.load(path_required(args.model,'--model'))
    x=np.asarray([v for r,v,y in examples]);y=np.asarray([y for r,v,y in examples]);pred=np.clip(model.predict(x),0,1)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    csvout(out/'external_predictions.csv',[{'construct_id':examples[i][0]['construct_id'],'target_group':examples[i][0]['target_group'],'measured':float(y[i]),'predicted':float(pred[i])} for i in range(len(y))],['construct_id','target_group','measured','predicted'])
    csvout(out/'excluded.csv',excluded,['construct_id','reason'])
    baseline=np.repeat(float(meta['training_mean_editing_fraction']),len(y))
    jsonout(out/'independent_validation.json',{'status':'EXTERNAL_EVALUATION_ONLY_NOT_KRAS_CLINICAL_VALIDATION','n_samples':len(y),'n_independent_targets':len(holdout_groups),'MAE':float(mean_absolute_error(y,pred)),'RMSE':float(np.sqrt(mean_squared_error(y,pred))),'train_mean_baseline_MAE':float(mean_absolute_error(y,baseline)),'target_group_overlap':[],'warning':'External accuracy in a novel target or cell context does not prove KRAS G12D therapeutic efficiency'})

if __name__=='__main__':
    p=argparse.ArgumentParser();s=p.add_subparsers(dest='cmd',required=True)
    a=s.add_parser('train');a.add_argument('--labels',required=True);a.add_argument('--chemistry',required=True);a.add_argument('--out',required=True);a.add_argument('--min-samples',type=int,default=40);a.add_argument('--min-groups',type=int,default=4);a.add_argument('--cv-folds',type=int,default=5)
    a=s.add_parser('infer');a.add_argument('--model',required=True);a.add_argument('--chemistry',required=True);a.add_argument('--candidates',required=True);a.add_argument('--out',required=True)
    a=s.add_parser('validate-external');a.add_argument('--model',required=True);a.add_argument('--provenance',required=True);a.add_argument('--labels',required=True);a.add_argument('--chemistry',required=True);a.add_argument('--out',required=True)
    args=p.parse_args();{'train':train,'infer':infer,'validate-external':external_validate}[args.cmd](args)
