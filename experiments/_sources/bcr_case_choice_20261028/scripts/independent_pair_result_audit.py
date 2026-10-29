"""Independent checks of completed pair outputs; no model calls or rule changes."""
from pathlib import Path
from collections import defaultdict
import csv,hashlib,json,math,statistics

R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
def read(name):
    with (R/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def flag(x):
    assert x in ['True','False'],x
    return x=='True'
def near(a,b,tol=2e-7):
    assert math.isclose(float(a),float(b),rel_tol=tol,abs_tol=tol),(a,b)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
raw=read('paired_all_residue_predictions.csv');number=read('paired_numbering.csv')
number_map={(n['state_id'],n['chain'],n['author'],n['aa']):n for n in number}
pred={(r['state_id'],r['variant'],r['chain'],r['author']):r for r in raw};assert len(pred)==len(raw)
for r in raw:
    assert flag(r['hchi_valid'])==(flag(r['model_valid']) and flag(r['geometry_valid']) and r['aa'] not in ['A','G'])
contrasts=json.loads((R/'pair_screen_contrasts.json').read_text());conmap={c['id']:c for c in contrasts}
delta=read('paired_complete_position_deltas.csv');regions=read('paired_predefined_region_deltas.csv')
by_group=defaultdict(list);changed=[]
for r in delta:
    con=conmap[r['contrast']];v=r['variant'];ch=r['chain'];au=r['author']
    b=pred[con['before'],v,ch,au];a=pred[con['after'],v,ch,au]
    assert r['aa_before']==b['aa'] and r['aa_after']==a['aa']
    same=b['aa']==a['aa'];assert flag(r['same_aa'])==same
    eligible=all(flag(x[y]) for x in [a,b] for y in ['model_valid','geometry_valid'])
    chi=eligible and flag(a['hchi_valid']) and flag(b['hchi_valid']) and same
    assert flag(r['rmsf_eligible'])==eligible and flag(r['hchi_eligible'])==chi
    for m in ['rmsf','hchi']:
        near(r[m+'_before'],b[m]);near(r[m+'_after'],a[m]);near(r['delta_'+m],float(a[m])-float(b[m]))
    if not same:
        assert ch=='H' and au=='103' and {a['aa'],b['aa']}=={'P','A'}
        assert not chi
        changed.append({'contrast':r['contrast'],'variant':v,'chain':ch,'author':au,'before':b['aa'],'after':a['aa'],'hchi_excluded':not chi})
    by_group[r['contrast'],v].append(r)

# Verify coverage by set intersection, not row count or order.
coverage=[]
for c in contrasts:
    for v in ['full','intrachain']:
        b={(r['chain'],r['author']) for r in raw if r['state_id']==c['before'] and r['variant']==v}
        a={(r['chain'],r['author']) for r in raw if r['state_id']==c['after'] and r['variant']==v}
        got={(r['chain'],r['author']) for r in by_group[c['id'],v]}
        assert got==a&b and len(got)==len(by_group[c['id'],v])
        coverage.append({'contrast':c['id'],'variant':v,'shared_positions':len(got),'before_only':sorted(b-a),'after_only':sorted(a-b)})

def select_scope(row,scope):
    ch=row['chain'];author=int(row['author_number_before'])
    if scope=='all_shared':return True
    if scope=='all_unchanged':return flag(row['same_aa'])
    if scope=='complete_variable_domains':return bool(row['imgt_number_before']) and bool(row['imgt_number_after'])
    if scope.startswith('chain_'):return ch==scope[-1]
    if scope in [x+'_'+y for x in ['H','L'] for y in ['CDR1','CDR2','CDR3']]:
        return ch==scope[0] and row['region_before']==scope[2:] and row['region_after']==scope[2:]
    if scope=='predefined_network_three_sites':return (ch,author) in [('H',31),('H',97),('A',67)]
    if scope=='predefined_network_complete_neighbors_plusminus2':return any(ch==c and abs(author-n)<=2 for c,n in [('H',31),('H',97),('A',67)])
    raise ValueError(scope)
for s in regions:
    m=s['metric'];selected=[r for r in by_group[s['contrast'],s['variant']] if select_scope(r,s['scope']) and flag(r[m+'_eligible'])]
    assert int(s['n'])==len(selected)>0
    for target,key in [('mean_before',m+'_before'),('mean_after',m+'_after'),('mean_delta','delta_'+m)]:near(s[target],statistics.mean(float(r[key]) for r in selected))
    near(s['mean_absolute_delta'],statistics.mean(abs(float(r['delta_'+m])) for r in selected))
    assert int(s['n_increase'])==sum(float(r['delta_'+m])>0 for r in selected)
    assert int(s['n_decrease'])==sum(float(r['delta_'+m])<0 for r in selected)

network_rows=read('paired_predefined_network_all_outputs.csv')
expected_network=[r for r in delta if (r['chain'],int(r['author_number_before'])) in [('H',31),('H',97),('A',67)]]
assert len(network_rows)==len(expected_network)
netmap={(r['contrast'],r['variant'],r['chain'],r['author']):r for r in network_rows}
for r in expected_network:
    n=netmap[r['contrast'],r['variant'],r['chain'],r['author']]
    assert n==r

observed_network=[]
for contrast in ['P103A_minus_WT_copy1','P103A_minus_WT_copy2','P103A_minus_WT_copy1_shared','P103A_minus_WT_copy2_shared']:
    for r in by_group[contrast,'full']:
        if (r['chain'],int(r['author_number_before'])) in [('H',31),('H',97),('A',67)]:
            assert r['aa_before']==r['aa_after']=={('H','31'):'D',('H','97'):'E',('A','67'):'R'}[r['chain'],r['author']]
            observed_network.append({'contrast':contrast,'site':r['chain']+r['author']+r['aa_before'],'delta_rmsf':float(r['delta_rmsf']),'delta_hchi':float(r['delta_hchi'])})
by_net={(r['contrast'],r['site']):r for r in observed_network}
network_signs=[]
for site in ['H31D','H97E','A67R']:
    a=by_net['P103A_minus_WT_copy1',site];b=by_net['P103A_minus_WT_copy2',site]
    network_signs.append({'site':site,'RMSF_sign_same_between_copies':a['delta_rmsf']*b['delta_rmsf']>0,
                          'Hchi_sign_same_between_copies':a['delta_hchi']*b['delta_hchi']>0})

h3=[s for s in regions if s['variant']=='full' and s['scope']=='H_CDR3' and s['contrast'] in ['P103A_minus_WT_copy1','P103A_minus_WT_copy2','P103A_minus_WT_copy1_shared','P103A_minus_WT_copy2_shared','WT_copy2_minus_copy1','P103A_copy2_minus_copy1']]
absolute=[r for r in read('6A77_complete_region_metrics.csv') if r['scope']=='heavy_CDR3']
report={'status':'passed','method':'Independent csv-module recomputation from raw before/after predictions and numbering; no original analysis function called.',
        'raw_prediction_rows_checked':len(raw),'delta_rows_checked':len(delta),'region_summaries_recomputed':len(regions),'network_rows_crosschecked':len(network_rows),
        'all_shared_position_intersections_exact':True,'all_delta_values_and_eligibility_match':True,'all_changed_H103_Hchi_excluded':all(x['hchi_excluded'] for x in changed),
        'changed_rows':changed,'coverage_by_contrast':coverage,'predefined_network_full_model':observed_network,'network_copy_direction_check':network_signs,
        'complete_H_CDR3_deltas_and_copy_controls':h3,'WT_Fab_6A77_absolute_HCDR3_metrics_not_mutation_truth':absolute,
        'interpretation_zh':[
            'H31D与H97E的RMSF和Hχ变化方向均在两份晶体拷贝间反转；共同覆盖裁剪没有解决这种反转。不能挑第一拷贝或第二拷贝来证明一种固定的盐桥稳定与柔性变化对应关系。',
            '抗原A67R的RMSF与Hχ在两拷贝都下降，但该一致局部信号没有WT/突变匹配动态真值，文献R477的编号依据也尚未解决；它不足以支持整条亲和机制论证。',
            '重链CDR3的平均RMSF突变差近于零且两拷贝符号不同，均值差小于同一状态的拷贝均值差。逐位平均绝对突变差约0.017–0.021Å，与拷贝差约0.012–0.014Å同一量级，不能说区域完全无局部响应，也不能说其均值稳定刚化。',
            '重链CDR3同身份侧链Hχ平均增加约0.021–0.032 nats，两拷贝同向；这不等于结合熵、溶剂熵或亲和力改善，也没有相应MD真值验证。',
            'WT Fab 6A77的H-CDR3绝对预测确实较好（RMSF MAE0.116 vs链内0.290Å，Hχ MAE0.160 vs0.255nats），但它是Fab结合态的单状态参考，不能定量验证scFv WT/P103A差值。',
            '当前不适合写成已经解释单点突变亲和提高的成功案例。可保留为WT抗体识别环动态预测案例与待检验突变线索；若要机制主案例，需要可比构建体的独立动态参考和稳健局部响应。'],
        'source_sha256':{str(p):sha(p) for p in [R/'paired_all_residue_predictions.csv',R/'paired_complete_position_deltas.csv',R/'paired_predefined_region_deltas.csv',R/'paired_predefined_network_all_outputs.csv',R/'6A77_complete_region_metrics.csv',R/'pair_screen_contract.json',Path(__file__)]}}
(R/'audits/paired_result_independent_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','raw_prediction_rows_checked','delta_rows_checked','region_summaries_recomputed','network_copy_direction_check','interpretation_zh']},ensure_ascii=False,indent=2))
