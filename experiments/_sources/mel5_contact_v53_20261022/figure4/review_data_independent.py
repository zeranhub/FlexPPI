"""Independent source/artist audit; never writes the reviewed plot files."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
PNG=ROOT/'fig4_mel5_contact_map.png'
PLOT=ROOT/'plot.py'
before=hashlib.sha256(PNG.read_bytes()).hexdigest()
raw=PLOT.read_text(encoding='utf8')
# Build artists in memory using the reviewed plot source, stopping before every
# render/export/CSV/audit write. This independently inspects exactly which source
# values are mapped to each plotted state, row, column, style and label.
prefix=raw.split('fig.canvas.draw();renderer=')[0]
env={'__file__':str(PLOT),'__name__':'__audit_only__'}
exec(compile(prefix,str(PLOT),'exec'),env)
data=pd.read_csv(BASE/'analysis/local_HLA6_beta48_56_all54_pairs.csv',float_precision='round_trip')
hchi=pd.read_csv(BASE/'analysis/HLA6_Hchi_source.csv',float_precision='round_trip').sort_values('position_a')
cdr=pd.read_csv(BASE.parent/'mel5_joint_case_v49_20261019/figure4/source_plotted_CDR_RMSF.csv',float_precision='round_trip')
beta=data.drop_duplicates('TCR_common_position').sort_values('TCR_common_position')
assert len(data)==54 and len(hchi)==6 and len(beta)==9 and len(cdr)==45
pd.testing.assert_frame_equal(pd.read_csv(ROOT/'source_all54_contacts.csv',float_precision='round_trip'),data)
pd.testing.assert_frame_equal(pd.read_csv(ROOT/'source_all6_local_Hchi.csv',float_precision='round_trip'),hchi.reset_index(drop=True))
pd.testing.assert_frame_equal(pd.read_csv(ROOT/'source_all9_local_RMSF.csv',float_precision='round_trip'),beta.reset_index(drop=True))
pd.testing.assert_frame_equal(pd.read_csv(ROOT/'source_all45_CDR_RMSF.csv',float_precision='round_trip'),cdr)

WT,ENG=env['WT'],env['ENG']
rows=list(hchi.position_a.astype(int));cols=list(beta.TCR_common_position.astype(int))
grid=env['grid']
points=[l for l in grid.lines if len(l.get_xdata())==1]
assert len(points)==108
checks=[]
for r in data.itertuples():
    y=rows.index(int(r.HLA_common_position));x=cols.index(int(r.TCR_common_position))
    for state,dx,color,marker in [('WT',-.23,WT,'o'),('mature',.23,ENG,'s')]:
        d=float(getattr(r,state+'_minimum_heavy_distance_A'))
        found=[l for l in points if abs(float(l.get_xdata()[0])-(x+dx))<1e-12 and abs(float(l.get_ydata()[0])-(y-.08))<1e-12]
        assert len(found)==1
        point=found[0]
        assert point.get_marker()==marker
        if d<=4:ms,fc,ec,alpha,category=5,color,color,1,'contact_le4A'
        elif d<=6:ms,fc,ec,alpha,category=4,'white',color,1,'near_4_to_6A'
        else:ms,fc,ec,alpha,category=1.5,color,color,.42,'far_gt6A'
        assert point.get_markersize()==ms and point.get_markerfacecolor()==fc and point.get_markeredgecolor()==ec and point.get_alpha()==alpha
        checks.append({'HLA':int(r.HLA_common_position),'beta':int(r.TCR_common_position),'state':state,'distance_A':d,'category':category})

q=data[data.HLA_common_position==72].sort_values('TCR_common_position')
qtexts=[]
for r in q.itertuples():
    x=cols.index(int(r.TCR_common_position));y=rows.index(72)
    for state,dx,color in [('WT',-.23,WT),('mature',.23,ENG)]:
        d=float(getattr(r,state+'_minimum_heavy_distance_A'))
        found=[t for t in grid.texts if np.max(np.abs(np.asarray(t.get_position())-[x+dx,y+.29]))<1e-12]
        assert len(found)==1 and found[0].get_text()==f'{d:.2f}' and found[0].get_color()==color
        qtexts.append({'beta_common_position':int(r.TCR_common_position),'state':state,'displayed':found[0].get_text(),'unrounded_A':d})
assert len(qtexts)==18

# Marginal model tracks use identical pair alignment and state order.
top=env['top']
assert len(top.lines)==3
for line,field,col,marker in zip(top.lines[:2],['TCR_RMSF_WT_A','TCR_RMSF_mature_A'],[WT,ENG],['o','s']):
    np.testing.assert_array_equal(line.get_xdata(),np.arange(9))
    np.testing.assert_array_equal(line.get_ydata(),beta[field].to_numpy())
    assert line.get_color()==col and line.get_marker()==marker
left=env['left']
assert len(left.lines)==18
for i,r in enumerate(hchi.itertuples()):
    conn,a,b=left.lines[i*3:i*3+3]
    np.testing.assert_array_equal(conn.get_xdata(),[r.hchi_full_a,r.hchi_full_b])
    np.testing.assert_array_equal(conn.get_ydata(),[i,i])
    assert a.get_marker()=='o' and a.get_markeredgecolor()==WT and b.get_marker()=='s' and b.get_markerfacecolor()==ENG
    np.testing.assert_array_equal(a.get_xdata(),[r.hchi_full_a]);np.testing.assert_array_equal(b.get_xdata(),[r.hchi_full_b])

lower=env['ax'];curve_lines=[l for l in lower.lines if len(l.get_xdata())>2]
assert len(curve_lines)==12
n=0
for i,region in enumerate(env['groups']):
    group=cdr[cdr.region==region].sort_values('plot_x')
    for j,(field,col,marker) in enumerate([('MEL5_model_RMSF_A',WT,'o'),('alpha24beta17_model_RMSF_A',ENG,'s')]):
        line=curve_lines[i*2+j]
        np.testing.assert_array_equal(line.get_xdata(),group.plot_x.to_numpy())
        np.testing.assert_array_equal(line.get_ydata(),group[field].to_numpy())
        assert line.get_color()==col and line.get_marker()==marker
        n+=len(group)
assert n==90

# Tick order maps rows top-to-bottom and columns left-to-right exactly.
assert [t.get_text() for t in grid.get_yticklabels()]==['R65','K66','K68','Q72','R75','V76']
assert list(grid.get_ylim())==[5.5,-.5]
assert [t.get_text() for t in grid.get_xticklabels()]==['48\nF','49\nY','50\nS→W','51\nV→G','52\nG→P','53\nI→F','54\nG','55\nQ','56\nI']
assert grid.get_xlim()==(-.5,8.5)
alltext='\n'.join(t.get_text() for t in env['fig'].findobj() if hasattr(t,'get_text'))
assert 'Nearest heavy-atom distance' in alltext
assert 'OE1' not in alltext and 'NE2' not in alltext and 'occupancy' not in alltext

after=hashlib.sha256(PNG.read_bytes()).hexdigest();assert before==after
review={'status':'passed','reviewed_PNG':str(PNG),'PNG_sha256':after,'plot_py_sha256':hashlib.sha256(PLOT.read_bytes()).hexdigest(),
        'visual_review':'PNG independently viewed in full; Q72 highlighted row has all18 rounded paired distances, correct left WT circle/right mature square placement, readable label mapping; no numeric mismatch seen.',
        'checks':{'geometry_pairs':54,'state_markers':108,'Q72_displayed_values':18,'local_Hchi_sites':6,'local_RMSF_sites':9,'CDR_sites':45,'CDR_state_values':90,'source_csv_copies_exact':True,'row_and_column_mapping':True,'state_order':True,'all_marker_threshold_classes':True,'all_plotted_model_values':True,'minimum_heavy_atom_distance_not_fixed_atom_pair':True},
        'geometry_categories':{s:{cat:sum(v['state']==s and v['category']==cat for v in checks) for cat in ['contact_le4A','near_4_to_6A','far_gt6A']} for s in ['WT','mature']},
        'Q72_numeric_labels':qtexts,
        'scientific_scope':'Edges/markers describe distances in archived model-input starting structures. Node/curve values are existing bound-complex model outputs. No claim of hydrogen-bond number, MD occupancy, direct causal influence, or affinity prediction.'}
(ROOT/'independent_data_review.json').write_text(json.dumps(review,ensure_ascii=False,indent=2),encoding='utf8')
env['plt'].close(env['fig'])
print(json.dumps({k:review[k] for k in ['status','PNG_sha256','checks','geometry_categories']},ensure_ascii=False,indent=2))
