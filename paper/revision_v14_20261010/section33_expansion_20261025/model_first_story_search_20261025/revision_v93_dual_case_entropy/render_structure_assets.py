"""Native PyMOL: unchanged real coordinates, cartoon overviews, complete-residue local sticks."""
from pathlib import Path
import os, json, hashlib, shutil, math
import numpy as np
from PIL import Image
ENV=Path('D:/Anaconda/envs/pvb_cpu')
os.environ['PATH']=str(ENV/'Library/bin')+os.pathsep+os.environ.get('PATH','')
DLL_HANDLES=[os.add_dll_directory(str(x)) for x in [ENV/'Library/bin',ENV/'DLLs'] if x.exists()]
import pymol2

ROOT=Path(__file__).resolve().parent
S=ROOT.parent
OUT=ROOT/'pymol';OUT.mkdir(exist_ok=True)
INP=OUT/'inputs';INP.mkdir(exist_ok=True)
CAM_SOURCE=S/'revision_v92_ribbon_overview_local_sticks/pymol/inputs'
for src in [CAM_SOURCE/'1MXE_CaM_native4_147_CaMKI294_318.pdb',S/'revision_v80_redesign/pymol/inputs/A003R_input_structure.pdb']:
    if not (INP/src.name).exists():shutil.copy2(src,INP/src.name)
COLORS=dict(partner='#BAC1C7',target_chain='#8B72A5',key='#3A8D8C',context='#D5DADD',neighbor='#A7B0B6')
FOV=20.;W,H=1600,1400
AA={'GLN':'Q','TRP':'W','MET':'M'}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def rgb(h):return [int(h[i:i+2],16)/255 for i in [1,3,5]]
def setup(c):
    c.reinitialize()
    for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_shadows=0,ray_trace_mode=0,orthoscopic=1,field_of_view=FOV,depth_cue=0,antialias=1,specular=.12,shininess=30,ambient=.72,direct=.30,reflect=.08,light_count=2,cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,cartoon_smooth_loops=1,cartoon_sampling=24,cartoon_loop_radius=.22,cartoon_tube_radius=.22,cartoon_oval_length=1.25,cartoon_oval_width=.25,max_threads=4,valence=0,stick_quality=24,transparency_mode=2,two_sided_lighting=1).items():c.set(k,v)
    for k,v in COLORS.items():c.set_color(k,rgb(v))
    c.bg_color('white')
def view_for(xyz,R,padding=1.5):
    center=xyz.mean(0);q=(xyz-center)@R
    shift=(q.min(0)+q.max(0))/2;center=center+R@shift;q=(xyz-center)@R
    hh=max(abs(q[:,1]).max(),abs(q[:,0]).max()/(W/H))+padding
    d=hh/math.tan(math.radians(FOV/2))
    return list(R.ravel())+[0.,0.,-float(d)]+[float(x) for x in center]+[1.,float(d+200),-20.],hh
def centered_view(center,R,hh):
    d=hh/math.tan(math.radians(FOV/2));return list(R.ravel())+[0.,0.,-float(d)]+list(map(float,center))+[1.,float(d+200),-20.]
def project(xyz,v):
    v=np.asarray(v);q=(np.asarray(xyz)-v[12:15])@v[:9].reshape(3,3)+v[9:12]
    hh=-v[11]*math.tan(math.radians(FOV/2));return np.column_stack([(.5+q[:,0]/(2*hh*W/H))*W,(.5-q[:,1]/(2*hh))*H])
def render(c,name,view,selection,landmarks):
    c.set_view(view);c.viewport(W,H);c.png(str(OUT/f'{name}.png'),width=W,height=H,dpi=800,ray=1,quiet=1);c.save(str(OUT/f'{name}.pse'))
    xyz=c.get_coords(selection);px=project(xyz,view)
    assert (px[:,0]>=0).all() and (px[:,0]<W).all() and (px[:,1]>=0).all() and (px[:,1]<H).all()
    lm={key:project(c.get_coords(sel).mean(0)[None,:],view)[0].tolist() for key,sel in landmarks.items()}
    im=Image.open(OUT/f'{name}.png').convert('RGBA');white=Image.new('RGBA',im.size,'white');white.alpha_composite(im);white.convert('RGB').save(OUT/f'{name}_white.png')
    return dict(png=f'{name}.png',pse=f'{name}.pse',width=W,height=H,view=list(c.get_view()),landmarks=lm,sha256=sha(OUT/f'{name}.png'),selected_atom_count=len(xyz),selected_atoms_in_frame=True)

manifest=dict(colors=COLORS,assets={},inputs={},native_structure_renders=True,no_predicted_conformations=True)
with pymol2.PyMOL() as pm:
    c=pm.cmd
    if (INP/'camera_rotations.json').exists():
        cameras=json.loads((INP/'camera_rotations.json').read_text());UR=np.asarray(cameras['ub']).reshape(3,3);CR=np.asarray(cameras['cam']).reshape(3,3)
    else:
        source=S/'revision_v80_redesign/supplement_ubiquitin/structure_renders/ubiquitin_Q41_Q62_recolored.pse'
        c.load(str(source));c.scene('01_complex','recall');UR=np.asarray(c.get_view()[:9]).reshape(3,3).copy()
        CR=np.asarray(json.loads((CAM_SOURCE/'base_camera_rotation.json').read_text())['rotation']).reshape(3,3)
        (INP/'camera_rotations.json').write_text(json.dumps(dict(ub=UR.ravel().tolist(),cam=CR.ravel().tolist()),indent=2),encoding='utf-8')
    for case,pdb,chain,keys,R in [('cam','1MXE_CaM_native4_147_CaMKI294_318.pdb','B',[303,316],CR),('ub','A003R_input_structure.pdb','A',[41,62],UR)]:
        setup(c);c.load(str(INP/pdb),'complex');c.remove('complex and (hydro or not polymer.protein)');c.dss('complex');c.hide('everything','complex')
        original=c.get_coords('complex').copy();manifest['inputs'][case]=dict(pdb='inputs/'+pdb,sha256=sha(INP/pdb),protein_heavy_atoms=len(original),target_chain=chain)
        c.show('cartoon','complex');c.color('partner','complex');c.color('target_chain',f'complex and chain {chain}')
        c.set('cartoon_transparency',.28,f'complex and not chain {chain}')
        c.color('key',f'complex and chain {chain} and resi '+ '+'.join(map(str,keys)))
        lm={AA[c.get_model(f'complex and chain {chain} and resi {n} and name CA').atom[0].resn]+str(n):f'complex and chain {chain} and resi {n} and name CA' for n in keys}
        view,hh=view_for(original,R)
        result=render(c,case+'_overview',view,'complex',lm);result['stick_atom_count']=c.count_atoms('complex and rep sticks');assert result['stick_atom_count']==0
        manifest['assets'][case+'_overview']=result
        if case=='ub':
            model=c.get_model('complex');atoms=model.atom;xyz=np.array([a.coord for a in atoms]);neigh=[]
            for site in keys:
                sc=c.get_coords(f'complex and chain A and resi {site} and not name N+CA+C+O')
                d=np.linalg.norm(xyz[:,None,:]-sc[None,:,:],axis=2).min(1)
                ids=sorted(set((a.chain,a.resi,a.resn) for a,dist in zip(atoms,d) if dist<=5.0))
                clauses=['(chain '+ch+' and resi '+'+'.join(n for cc,n,rn in ids if cc==ch)+')' for ch in sorted(set(a[0] for a in ids))]
                sel='complex and ('+' or '.join(clauses)+')';near=c.get_coords(sel);v,hh=view_for(near,R,padding=1.3)
                neigh.append(dict(site=site,selection=sel,whole_residues=ids,center=v[12:15],halfheight=hh))
            hh=max(n['halfheight'] for n in neigh)
            for n in neigh:
                site=n['site'];sel=n['selection'];target=f'complex and chain A and resi {site}'
                c.hide('everything','complex');c.show('cartoon','complex');c.color('context','complex');c.set('cartoon_transparency',.88,'complex')
                c.show('sticks',sel);c.color('neighbor',sel);c.color('key',target);c.set('stick_radius',.115,sel);c.set('stick_radius',.20,target)
                view=centered_view(n['center'],R,hh)
                result=render(c,f'ub_Q{site}_local',view,sel,{f'Q{site}':target+' and not name N+CA+C+O'})
                check_coords=c.get_coords('complex');assert np.array_equal(original,check_coords)
                count=c.count_atoms(sel);sticks=c.count_atoms('complex and rep sticks');assert count==sticks
                result.update(neighborhood_cutoff_A=5.0,whole_residues=n['whole_residues'],all_whole_residue_heavy_atoms_shown_as_sticks=True,stick_atom_count=sticks,no_sticks_outside_selection=True,unchanged_source_coordinates=True,shared_optical_halfheight_A=hh,source_atom_count=len(original))
                manifest['assets'][f'ub_Q{site}_local']=result
        assert np.array_equal(original,c.get_coords('complex'))
(OUT/'render_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:dict(png=v['png'],atoms=v['selected_atom_count']) for k,v in manifest['assets'].items()}))
