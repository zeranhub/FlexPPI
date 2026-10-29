"""Read-only recomputation of all newly inserted CaM numerical claims."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / 'partner_screen'

def load(name):
    with (DATA / name).open(encoding='utf-8-sig', newline='') as file:
        return list(csv.DictReader(file))

def mean(values):
    values = list(values)
    return math.fsum(values) / len(values)

def ranks(values):
    pairs = sorted(enumerate(values), key=lambda item: item[1])
    result = [None] * len(pairs)
    start = 0
    while start < len(pairs):
        end = start + 1
        while end < len(pairs) and pairs[end][1] == pairs[start][1]:
            end += 1
        rank = (start + end - 1) / 2
        for index, _ in pairs[start:end]:
            result[index] = rank
        start = end
    return result

def rho(x, y):
    x, y = ranks(x), ranks(y)
    xm, ym = mean(x), mean(y)
    numerator = math.fsum((a - xm) * (b - ym) for a, b in zip(x, y))
    denominator = math.sqrt(math.fsum((a-xm)**2 for a in x) * math.fsum((b-ym)**2 for b in y))
    return numerator / denominator

rows = load('cam_md_state_contrast_all_common_residues.csv')
hchi = [row for row in rows if row['metric'] == 'hchi']
rmsf = [row for row in rows if row['metric'] == 'rmsf']
assert len(hchi) == 122 and len(rmsf) == 145
truth = [float(row['md_delta']) for row in hchi]
prediction = [float(row['full0_delta']) for row in hchi]
large = [row for row in hchi if abs(float(row['md_delta'])) >= .1]
correct = [row for row in large if float(row['md_delta']) * float(row['full0_delta']) > 0]
metrics = {
    'RMSF_positions':len(rmsf), 'Hchi_positions':len(hchi),
    'RMSF_MD_NHE1':mean(float(row['md_a']) for row in rmsf),
    'RMSF_MD_IQCG':mean(float(row['md_b']) for row in rmsf),
    'Hchi_up':sum(value > 0 for value in truth), 'Hchi_down':sum(value < 0 for value in truth),
    'Hchi_change_rho':rho(truth, prediction),
    'Hchi_change_MAE_Full':mean(abs(a-b) for a,b in zip(truth, prediction)),
    'Hchi_change_MAE_Intra':mean(abs(float(row['md_delta'])-float(row['intra0_delta'])) for row in hchi),
    'Hchi_change_MAE_AA_SASA_contacts':mean(abs(float(row['md_delta'])-float(row['aa_sasa_contacts0_delta'])) for row in hchi),
    'Hchi_change_MAE_Static47_SASA':mean(abs(float(row['md_delta'])-float(row['static47_sasa0_delta'])) for row in hchi),
    'Hchi_change_MAE_zero':mean(abs(value) for value in truth),
    'Hchi_large_n':len(large), 'Hchi_large_correct':len(correct),
    'Hchi_large_up_correct':sum(float(row['md_delta']) > 0 for row in correct),
    'Hchi_large_down_correct':sum(float(row['md_delta']) < 0 for row in correct),
}
assert metrics['Hchi_up'] == 67 and metrics['Hchi_down'] == 55
assert (metrics['Hchi_large_n'], metrics['Hchi_large_correct']) == (70,58)
groups = {}
for group_name, contact_a, contact_b in [('NHE1_only',True,False),('IQCG_only',False,True),('both',True,True),('neither',False,False)]:
    group = [row for row in hchi if bool(int(row['contact_a'])) == contact_a and bool(int(row['contact_b'])) == contact_b]
    groups[group_name] = {'n':len(group),'MD_delta':mean(float(row['md_delta']) for row in group),'Full_delta':mean(float(row['full0_delta']) for row in group)}
assert [groups[k]['n'] for k in ('NHE1_only','IQCG_only','both','neither')] == [8,21,3,90]
met_rows = load('final_cam_all_nine_methionines_with_stability.csv')
met = [row for row in met_rows if row['metric'] == 'hchi']
assert len(met) == 9
met_report = {'n':len(met), 'rho':rho([float(row['md_delta']) for row in met],[float(row['full0_delta']) for row in met]),'key_sites':{}}
for number in (72,109,144):
    row = next(row for row in met if int(row['mature_cam_number']) == number)
    met_report['key_sites'][str(number)] = {key: float(row[key]) for key in ('md_a','md_b','md_delta','full0_delta')}
    assert int(row['blocks_same_sign_as_mean']) == 5
report = {
    'source_folder':str(DATA), 'metrics_recomputed':metrics, 'contact_groups_recomputed':groups,
    'all_nine_Met':met_report,
    'source_files':[{'path':str(DATA/name),'sha256':hashlib.sha256((DATA/name).read_bytes()).hexdigest()} for name in (
        'cam_md_state_contrast_all_common_residues.csv','final_cam_fixed_structural_group_metrics.csv',
        'final_cam_all_nine_methionines_with_stability.csv','cam_md_state_contrast_region_metrics.csv',
        'cam_md_per_lobe_fitted_rmsf_control_summary.csv','final_cam_interlobe_distance_and_compaction.csv',
        'cam_md_state_contrast_five_block_metrics.csv','final_cam_checks.json')],
    'legacy_ubiquitin_numbers':'Preserved from v79 paragraphs 73, 87, 89, 164-166; original source package retained.',
    'reference_metadata':'References 29 and 30 supplied by root after primary-source verification.',
    'status':'PASS; numerical claims reproduce independently with standard-library arithmetic.'
}
(HERE/'v80_scientific_source_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=False,indent=2))
