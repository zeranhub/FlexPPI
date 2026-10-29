"""True-coordinate B4 structure layers: reference, prediction, signed error.
No new structure, dynamics, mutation or geometry is inferred by this renderer.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts/tcr_real_states'))
from render_structures import setup, camera_basis, set_camera, project, SPECS, pymol2
OUT=ROOT/'output/tcr_prediction_reference_revision_20260915/structure_assets'
DATA=ROOT/'output/tcr_prediction_reference_revision_20260915/validation_paired_residues.csv'
PDB=ROOT/'output/tcr_real_states_20260915/structure_assets/b423_holo_deposited.pdb'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(DATA,dtype={'author_id':str});df=df[df.evaluation_mask]
    audit=dict(pdb=str(PDB),pdb_sha256=hashlib.sha256(PDB.read_bytes()).hexdigest(),
        data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),assets=[],
        note='All layers use the SAME deposited 5IVX coordinates as a display carrier. MD layers are residue scalars, not an MD mean conformation. Grey: partners/unscored residues.')
    with pymol2.PyMOL() as pm:
        cmd=pm.cmd;setup(cmd);cmd.load(str(PDB),'holo');cmd.hide('everything');cmd.show('cartoon','holo');cmd.dss('holo')
        cmd.set_color('low',[.86,.92,.96]);cmd.set_color('mid',[.42,.68,.84]);cmd.set_color('high',[.03,.19,.42])
        cmd.set_color('neutral',[.81,.83,.84]);cmd.set_color('errblue',[.1294,.4,.6745]);cmd.set_color('errzero',[.9686,.9686,.9686]);cmd.set_color('errred',[.698,.094,.169])
        basis=camera_basis(cmd,SPECS['b423'])
        camera=set_camera(cmd,cmd.get_coords('holo and name CA'),basis,1400,1450,pad=1.04)
        local_coords=cmd.get_coords('holo and byres ((chain F and resi 127-138) around 5)')
        for metric in ['rmsf','hchi']:
            q=df[df.metric==metric]
            maxval=float(max(q.md_reference.max(),q.prediction.max()));upper=float(np.ceil(maxval*2)/2)
            error_upper=float(np.ceil(q.model_error.abs().max()*2)/2)
            for kind,column in [('md','md_reference'),('pred','prediction'),('error','model_error')]:
                cmd.hide('everything');cmd.show('cartoon','holo');cmd.color('neutral','holo');cmd.alter('holo','b=-999')
                for r in q.itertuples():
                    selection=f'holo and chain {r.chain} and resi {r.author_id}'
                    assert cmd.count_atoms(selection+' and name CA')==1
                    cmd.alter(selection,f'b={getattr(r,column):.9f}')
                lo,hi=(-error_upper,error_upper) if kind=='error' else (0.,upper)
                palette='errblue errzero errred' if kind=='error' else 'low mid high'
                cmd.spectrum('b',palette,'holo and b > -998',minimum=lo,maximum=hi)
                if metric=='rmsf':
                    cmd.show('spheres','holo and name CA and ((chain E and resi 196) or (chain F and resi 213))')
                    cmd.set('sphere_scale',.65,'holo')
                cmd.set_view(camera['view']);name=f'b423_{metric}_{kind}'
                cmd.png(str(OUT/(name+'.png')),width=1400,height=1450,dpi=300,ray=1)
                cmd.scene(name,'store')
                landmarks=[]
                for ch,rid in [('E','196'),('F','213'),('F','133'),('F','134')]:
                    xyz=cmd.get_coords(f'holo and chain {ch} and resi {rid} and name CA')[0]
                    landmarks.append(dict(chain=ch,author_id=rid,pixel=project(xyz,camera)))
                audit['assets'].append(dict(name=name,metric=metric,kind=kind,column=column,range=[lo,hi],n=len(q),camera=camera,landmarks=landmarks))
                if metric=='hchi':
                    cmd.hide('cartoon','holo and not (byres ((chain F and resi 127-138) around 5))')
                    cmd.show('cartoon','holo and chain F and resi 127-138')
                    cmd.show('sticks','holo and chain F and resi 133+134')
                    lc=set_camera(cmd,local_coords,basis,1400,950,pad=1.08)
                    lname=f'b423_local_hchi_{kind}'
                    cmd.png(str(OUT/(lname+'.png')),width=1400,height=950,dpi=300,ray=1);cmd.scene(lname,'store')
                    marks=[]
                    for rid in ['133','134']:
                        xyz=cmd.get_coords(f'holo and chain F and resi {rid} and name CA')[0]
                        marks.append(dict(chain='F',author_id=rid,pixel=project(xyz,lc)))
                    audit['assets'].append(dict(name=lname,metric=metric,kind=kind,column=column,range=[lo,hi],camera=lc,landmarks=marks))
                print('rendered',name,flush=True)
        cmd.save(str(OUT/'b423_prediction_reference_scenes.pse'))
    (OUT/'render_provenance.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')

if __name__=='__main__':main()
