from pathlib import Path
import json,importlib.util,numpy as np
from PIL import Image
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('render_base',P/'render_mechanism_structures.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
r=json.loads((P/'render_manifest.json').read_text())
with m.pymol2.PyMOL() as p:
 c=p.cmd;c.load(str(P/'1A06_context.pse'));c.disable('all');c.enable('kinase');before=c.get_coords('kinase').copy()
 c.hide('everything','kinase')
 # Local regulatory connection retained. The full observed kinase object remains in the session.
 c.set_color('zoomghost',[.86,.88,.90]);c.show('cartoon','kinase');c.color('zoomghost','kinase');c.set('cartoon_color','zoomghost','kinase');c.set('cartoon_transparency',.85,'kinase')
 c.create('reg_segment','kinase and resi 294-300');c.hide('everything','reg_segment');c.show('cartoon','reg_segment');c.color('peptide','reg_segment');c.set('cartoon_color','peptide','reg_segment');c.set('cartoon_transparency',.05,'reg_segment')
 # Retain transparent local catalytic context, rather than an isolated amino-acid drawing.
 local='kinase and byres (kinase within 9 of (kinase and resi 298)) and not resi 290-316'
 # Full-kinase ghost avoids artificial cut ends in the local cartoon.
 ss='kinase and resi 104+184+210+214+216 and not name N+C+O'
 c.show('sticks',ss);c.set('stick_radius',.14,ss);c.color('contact',ss)
 # G214 has no carbon side chain: explicitly retain its entire backbone, including contacting O.
 c.show('sticks','kinase and resi 214');c.set('stick_radius',.12,'kinase and resi 214');c.color('contact','kinase and resi 214')
 target='reg_segment and resi 298';c.show('sticks',target);c.set('stick_radius',.24,target);c.color('key',target)
 c.viewport(1400,1400);c.set_view(m.VIEW);c.zoom('(kinase and resi 104+184+210+214+216 and not name N+CA+C+O) or (kinase and resi 214) or (reg_segment and resi 298)',buffer=3.2,complete=1)
 v=list(c.get_view());v[:9]=m.VIEW[:9];v[15]=1.;v[16]=-v[11]+160.;c.set_view(v)
 c.png(str(P/'F298_pocket.png'),width=1400,height=1400,dpi=1200,ray=1,quiet=1)
 im=Image.open(P/'F298_pocket.png').convert('RGBA');bb=im.getchannel('A').getbbox();pad=35
 crop=[130,100,1230,1200]  # Square focus retains all six displayed residues; only translucent context extends past the frame.
 im.crop(crop).save(P/'F298_pocket_trim.png')
 sels={'F298':'kinase and resi 298 and name CZ','G214':'kinase and resi 214 and name O','I210':'kinase and resi 210 and name CG2','F104':'kinase and resi 104 and name CE2','Y184':'kinase and resi 184 and name OH','P216':'kinase and resi 216 and name CG'}
 marks={k:m.project(c,s,1400,1400) for k,s in sels.items()}
 a={'file':'F298_pocket.png','trim_file':'F298_pocket_trim.png','width':1400,'height':1400,'dpi':1200,'crop_box':crop,'trim_width':crop[2]-crop[0],'trim_height':crop[3]-crop[1],'view':list(c.get_view()),'landmarks':marks,'trim_landmarks':{k:[x[0]-crop[0],x[1]-crop[1]] for k,x in marks.items()},'geometry_audit':{'atom_count':len(before),'coordinates_exactly_unchanged':bool(np.array_equal(before,c.get_coords('kinase'))),'max_coordinate_change_A':float(abs(before-c.get_coords('kinase')).max())},'scope':'Static heavy-atom proximity only. No hydrogen-bond claim. G214 contact is to its backbone O; glycine has no sidechain.'}
 c.scene('F298_pocket','store');c.save(str(P/'F298_pocket.pse'))
 r['assets']['F298_pocket']=a
(P/'render_manifest.json').write_text(json.dumps(r,indent=2),encoding='utf-8')
(P/'F298_pocket_manifest.json').write_text(json.dumps(a,indent=2),encoding='utf-8')
print(json.dumps(a,indent=2))