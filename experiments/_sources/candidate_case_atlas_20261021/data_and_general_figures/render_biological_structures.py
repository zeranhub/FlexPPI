"""Matched-camera true-coordinate structure views; no illustrative mechanisms."""
from pathlib import Path
import json,csv,hashlib,sys
import numpy as np
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
ROOT=Path(__file__).resolve().parent
CAT=json.loads((ROOT.parent/'catalog/candidate_inventory.json').read_text(encoding='utf8'))
COLORS=['#CF8C27','#762A83','#2166AC','#1B7837','#B2182B','#666666','#4393C3','#B35806','#9970AB']

def ids(p,role,max_n=None):
    with (ROOT/current/'data'/f'{p}_absolute_state_values.csv').open(encoding='utf-8-sig') as f:
        rows=[r for r in csv.DictReader(f) if r['role']==role]
    if max_n:rows=rows[:max_n]
    return rows
def select(p,role,max_n=None):
    r=ids(p,role,max_n)
    return f"s{p} and chain {r[0]['chain']} and resi "+'+'.join(x['id'] for x in r)

for case in CAT['cases']:
    current=case['case_id']
    if not 5<=int(current[1:])<=17:continue
    out=ROOT/current/'structures';out.mkdir(exist_ok=True)
    marker=out/'molecular_view.json'
    cmd.reinitialize()
    a,b=case['contrasts'][0]
    allp=[a,b]+[p for p in case['pdbs'] if p not in {a,b}]
    if marker.exists() and json.loads(marker.read_text(encoding='utf8'))['all_states']==allp:
        print(current,'reuse existing render',flush=True);continue
    sources=[]
    for k,p in enumerate(allp):
        s=json.loads((ROOT/current/'data'/f'{p}_source_audit.json').read_text(encoding='utf8'))
        cmd.load(s['structure'],'s'+p)
        if p!=a:cmd.align(select(p,'mhc',180)+' and name CA',select(a,'mhc',180)+' and name CA',cycles=0)
        sources.append(s)
        cmd.set_color('state'+str(k),[int(COLORS[k%len(COLORS)].lstrip('#')[q:q+2],16)/255 for q in (0,2,4)])
    cmd.hide('everything','all');cmd.bg_color('white')
    cmd.set('orthoscopic',1);cmd.set('ray_shadows',0);cmd.set('specular',.18);cmd.set('antialias',2)
    cmd.set('ray_trace_mode',0);cmd.set('cartoon_fancy_helices',1);cmd.set('cartoon_smooth_loops',1)
    cmd.set('max_threads',4)
    cmd.set_color('hla',[.72,.79,.79]);cmd.set_color('alpha',[.63,.66,.69]);cmd.set_color('beta',[.39,.43,.47])
    pep=np.asarray(cmd.get_coords(select(a,'peptide')+' and name CA'))
    x=pep[-1]-pep[0];x=x/np.linalg.norm(x)
    tsel='('+select(a,'tcr_alpha',125)+') or ('+select(a,'tcr_beta',125)+')'
    t=np.asarray(cmd.get_coords('('+tsel+') and name CA'))
    z=t.mean(0)-pep.mean(0);z-=x*np.dot(z,x);z=z/np.linalg.norm(z)
    y=np.cross(z,x);zview=.55*z-.84*y;yview=np.cross(zview,x)
    rot=np.stack([x,yview,zview])
    for k,p in enumerate(allp):
        cmd.show('cartoon',select(p,'mhc',180));cmd.color('hla',select(p,'mhc',180))
        for role,col in [('tcr_alpha','alpha'),('tcr_beta','beta')]:
            sel=select(p,role,125);cmd.show('cartoon',sel);cmd.color(col,sel)
        sel=select(p,'peptide')+' and not hydro';cmd.show('sticks',sel);cmd.color('state'+str(k),sel)
        cmd.set('stick_radius',.16,sel)
    cmd.set_view(tuple(rot.T.flatten())+(0,0,-150)+tuple(pep.mean(0)+z*9)+(100,200,0))
    cmd.zoom('s'+a+' and visible',buffer=2,complete=1)
    view=list(cmd.get_view())
    for p in [a,b]:
        cmd.disable('all');cmd.enable('s'+p);cmd.set_view(view)
        cmd.png(str(out/f'{p}_overview.png'),1100,850,dpi=450,ray=1)
    cmd.enable('all');cmd.hide('everything','all')
    for k,p in enumerate(allp):
        sel=select(p,'peptide')+' and not hydro';cmd.show('sticks',sel);cmd.color('state'+str(k),sel);cmd.set('stick_radius',.14,sel)
    cmd.set_view(tuple(rot.T.flatten())+(0,0,-90)+tuple(pep.mean(0))+(50,150,0))
    cmd.zoom('visible',buffer=2,complete=1)
    cmd.png(str(out/'all_state_peptide_superposition.png'),1200,750,dpi=450,ray=1)
    cmd.save(str(out/'structure_views.pse'))
    marker.write_text(json.dumps(dict(status='rendered',case_id=current,reference=a,primary_pair=[a,b],all_states=allp,
           source='Existing per-state source-audit coordinates: ITrace run1 topology or CatchBond experimental-matched PDB; no newly simulated ensemble',sources=sources,
           alignment='MHC first180 C-alpha, PyMOL align cycles=0; common camera for primary pair',
           semantic_colors=COLORS,shared_view=view,no_invented_mechanism=True),ensure_ascii=False,indent=2),encoding='utf8')
    print(current,'molecular views complete',flush=True)
cmd.quit()
