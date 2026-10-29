"""Local continuous-lobe detail: complete cartoon selections, closer framing."""
from render_structures_v80 import *

def face_sidechain(c, obj, site):
    ca=c.get_coords(f'{obj} and chain A and resi {site} and name CA')[0]
    ce=c.get_coords(f'{obj} and chain A and resi {site} and name CE')[0]
    x=(ce-ca);x=x/np.linalg.norm(x)
    coords=c.get_coords(f'{obj} and chain A and resi '+('4-75' if site==72 else '82-148')+' and name CA')
    _,_,vh=np.linalg.svd(coords-coords.mean(0))
    y=vh[0]-np.dot(vh[0],x)*x;y=y/np.linalg.norm(y)
    z=np.cross(x,y);z=z/np.linalg.norm(z)
    near=c.get_coords(f'{obj} and chain P and byres (({obj} and chain P) within 5 of ({obj} and chain A and resi {site} and not name N+CA+C+O))')
    if near is not None and np.dot(z,near.mean(0)-ce)>0:
        z=-z;y=-y
    view=list(c.get_view());view[:9]=np.column_stack([x,y,z]).ravel().tolist();c.set_view(view)
    c.turn('z',-15)

def main():
    manifest=json.loads((OUT/'render_manifest.json').read_text())
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(str(OUT/'CaM_continuous_context_v80.pse'))
        c.viewport(2400,2400)
        for site,domain in [(72,'4-81'),(109,'76-148')]:
            left=f'NHE1_M{site}';right=f'IQCG_M{site}'
            face_sidechain(c,left,site)
            c.zoom(f'({left} or {right}) and chain A and resi {domain}',buffer=1.4,complete=1)
            view=c.get_view()
            for label,obj in [('NHE1',left),('IQCG',right)]:
                c.disable('all');c.enable(obj);c.set_view(view)
                c.set('cartoon_transparency',.65,obj+' and chain A')
                c.set('cartoon_transparency',.80,obj+' and chain A and resi '+('82-148' if site==72 else '4-75'))
                c.set('cartoon_transparency',.72,obj+' and chain P')
                c.set('stick_radius',.23,f'{obj} and chain A and resi {site}')
                c.set('stick_transparency',0,f'{obj} and chain A and resi {site}')
                c.show('spheres',f'{obj} and chain A and resi {site} and not name N+C+O')
                c.set('sphere_scale',.17,f'{obj} and chain A and resi {site}')
                c.clip('slab',150)
                focus=f'{obj} and chain A and resi {site}'
                name=f'detail_M{site}_{label}'
                marks={f'M{site}':focus+' and name CE',f'M{site}_CA':focus+' and name CA'}
                # Preserve residue identities for vector labels in composition.
                near=c.get_model(f'{obj} and chain P and byres (({obj} and chain P) within 4.3 of ({focus} and not name N+CA+C+O)) and name CA')
                for a in near.atom:marks[f'{a.resn}{a.resi}']=f'{obj} and chain P and resi {a.resi} and name CA'
                manifest['assets'][name]=render(c,name,2400,2400,marks)
                manifest['assets'][name]['framing']='Continuous complete-chain cartoon selection, framed on full N-lobe4–81 or C-lobe76–148; chain context can continue beyond the image border.'
            paired_trim(manifest['assets'],[f'detail_M{site}_NHE1',f'detail_M{site}_IQCG'])
        c.save(str(OUT/'CaM_continuous_context_v80.pse'))
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
