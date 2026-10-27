import sys
sys.path.insert(0,'scripts')
from restore2_reference import validate,guide_mapping
seq='ACCU';items=[]
for i,b in enumerate(seq,1):items.append(dict(construct_id='published_x',target_id='test',guide_5to3=seq,guide_position_1based=str(i),base=b,sugar='2OME',link_to_next_3prime='NONE' if i==len(seq) else 'PO',source_table='test',source_row=str(i),evidence_level='PUBLISHED_EXACT',is_orphan_C='YES' if i==2 else 'NO'))
assert len(validate(items))==1
trans,rej=guide_mapping(items,[dict(guide_id='test_guide',guide_rna_5to3='GCCU',guide_index_0based='1')]);assert len(trans)==4 and len(rej)==0
try:validate(items[:-1]);raise AssertionError('Should reject missing positions')
except ValueError:pass
print('RESTORE2 reference validation/transfer smoke: PASS')
