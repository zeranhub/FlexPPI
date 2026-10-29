"""Actual simulation/model input structures, mature-CaM mapping, PyMOL rays."""
from pathlib import Path
import os,json,csv,hashlib
import numpy as np
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
if hasattr(os,'add_dll_directory'):
    DLL_HANDLES=[os.add_dll_directory(str(p)) for p in [ENV/'Library/bin',ENV/'DLLs'] if p.exists()]
import pymol2
OUT=Path(__file__).resolve().parent
CORE=Path('E:/Project/FlexPPI/experiments/core_corrected_pipeline')
P=OUT.parent/'partner_screen'
SPECS=[('NHE1','A004C','B'),('IQCG','A007U','A')]
def rgb(h):return [int(h[i:i+2],16)/255 for i in (1,3,5)]
def prepare():
    (OUT/'inputs').mkdir(exist_ok=True);mapping=list(csv.DictReader((P/'cam_md_mature_number_mapping.csv').open()))
    audit=[]
    for label,sid,chain in SPECS:
        mm={r['local_author']:int(r['mature_cam_number']) for r in mapping if r['sample_id']==sid}
        source=CORE/'training_payload/dataset/structures'/(sid+'.pdb');lines=[]
        for line in source.read_text().splitlines():
            if line[:6].strip() not in ['ATOM','HETATM']:continue
            oldchain=line[21:22];author=line[22:26].strip()+line[26:27].strip()
            if oldchain==chain:
                if author not in mm:continue
                line=line[:21]+'A'+f'{mm[author]:4d}'+' '+line[27:]
            else:line=line[:21]+'P'+line[22:]
            lines.append(line)
        path=OUT/'inputs'/(label+'_mature_CaM_and_partner.pdb');path.write_text('\n'.join(lines)+'\nEND\n')
        audit.append(dict(partner=label,sample_id=sid,input_source=str(source),input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            output=str(path),output_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),coordinate_changes=False,
            changes='Only chain/residue identifiers changed; extra CaM expression-tail Ala149 excluded; no atoms rebuilt.'))
    (OUT/'structure_preparation_audit.json').write_text(json.dumps(audit,indent=2))
def project(c,sel,w,h):
    xyz=c.get_coords(sel).mean(0);v=np.array(c.get_view());q=v[:9].reshape(3,3).T@(xyz-v[12:15])+v[9:12]
    hh=-v[11]*np.tan(np.deg2rad(float(c.get('field_of_view'))/2));hw=hh*w/h
    return [float((.5+q[0]/(2*hw))*w),float((.5-q[1]/(2*hh))*h)]
def render(c,name,w,h,landmarks):
    c.png(str(OUT/'structure_renders'/(name+'.png')),width=w,height=h,dpi=400,ray=1,quiet=1)
    return dict(file=name+'.png',width=w,height=h,view=list(c.get_view()),landmarks={k:project(c,s,w,h) for k,s in landmarks.items()})
def style(c,obj):
    c.hide('everything',obj);c.show('cartoon',obj+' and polymer.protein');c.color('N_lobe',obj+' and chain A and resi 4-75')
    c.color('linker',obj+' and chain A and resi 76-81');c.color('C_lobe',obj+' and chain A and resi 82-148');c.color('partner',obj+' and chain P')
    sel=obj+' and chain A and resi 72+109+144';c.show('sticks',sel+' and not name N+C+O');c.color('key',sel+' and elem C');c.color('yellow',sel+' and elem S')
def main():
    prepare();(OUT/'structure_renders').mkdir(exist_ok=True);manifest=dict(assets={})
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.bg_color('white')
        for k,v in dict(ray_opaque_background=0,ray_trace_mode=0,ray_shadows=0,antialias=2,orthoscopic=1,depth_cue=0,specular=.14,ambient=.6,direct=.4,
                        cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=16,max_threads=4,stick_radius=.17,valence=0).items():c.set(k,v)
        for name,col in [('N_lobe','#4C83A1'),('C_lobe','#52A49A'),('linker','#B7BEC2'),('partner','#B87936'),('key','#AA2838')]:c.set_color(name,rgb(col))
        for label,_,_ in SPECS:
            c.load(str(OUT/'inputs'/(label+'_mature_CaM_and_partner.pdb')),label);c.remove(label+' and hydro');c.dss(label);style(c,label)
        ali=c.align('IQCG and chain A and resi 4-75 and name CA','NHE1 and chain A and resi 4-75 and name CA',cycles=0)
        manifest['N_lobe_alignment_rmsd_A']=float(ali[0])
        c.orient('all');c.turn('x',20);c.turn('y',20);c.turn('z',-12);c.zoom('all',buffer=3,complete=1);view=c.get_view()
        for label,_,_ in SPECS:
            c.disable('all');c.enable(label);c.set_view(view)
            marks={'M'+str(k):f'{label} and chain A and resi {k} and name CA' for k in [72,109,144]}
            marks.update(N_lobe=label+' and chain A and resi 4-75 and name CA',C_lobe=label+' and chain A and resi 82-148 and name CA',partner=label+' and chain P and name CA')
            manifest['assets']['whole_'+label]=render(c,'whole_'+label,1400,1500,marks);c.scene('whole_'+label,'store')
            print('RENDERED whole',label,flush=True)
        # Local copies preserve the full-view object coordinates and saved scenes.
        for site in [72,109]:
            for label,_,_ in SPECS:c.create(label+'_local'+str(site),label)
            l='NHE1_local'+str(site);r='IQCG_local'+str(site)
            c.align(f'{r} and chain A and resi {site-5}-{site+5} and name N+CA+C',f'{l} and chain A and resi {site-5}-{site+5} and name N+CA+C',cycles=0)
            focus=f'{l} and chain A and resi {site}';c.orient(focus);c.turn('y',35);c.turn('x',15);c.zoom(focus,buffer=6.5,complete=1);localview=c.get_view()
            for label,obj in [('NHE1',l),('IQCG',r)]:
                c.disable('all');c.enable(obj);c.hide('everything',obj)
                focus=f'{obj} and chain A and resi {site}'
                near=f'{obj} and byres (({obj}) within 9 of ({focus}))'
                c.show('cartoon',near);c.set('cartoon_transparency',.5,obj+' and chain A')
                c.show('sticks',focus+' and not name N+C+O');c.show('sticks',obj+' and chain P and byres (('+obj+' and chain P) within 4.5 of ('+focus+'))')
                c.color('key',focus+' and elem C');c.color('yellow',focus+' and elem S');c.set_view(localview)
                manifest['assets'][f'M{site}_{label}']=render(c,f'M{site}_{label}',900,1250,{'M'+str(site):focus+' and name CE'})
                c.scene(f'local_M{site}_{label}','store');print('RENDERED local',site,label,flush=True)
        c.save(str(OUT/'structure_renders'/'CaM_model_story.pse'))
        manifest['pymol_version']=str(c.get_version()[0]);manifest['scope']='Actual simulation inputs; static snapshots do not encode motion; no ions added.'
    (OUT/'structure_renders'/'render_manifest.json').write_text(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
