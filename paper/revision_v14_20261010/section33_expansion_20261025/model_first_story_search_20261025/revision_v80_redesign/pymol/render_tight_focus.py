"""Tight atom-readable local views; full cartoon selection, camera crop only."""
from render_structures_v80 import *
from render_local_details import face_sidechain

def main():
    manifest=json.loads((OUT/'render_manifest.json').read_text())
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(str(OUT/'CaM_continuous_context_v80.pse'))
        c.set_color('focus_deep',rgb('#276F6F'))
        c.set_color('neighbour_bluegrey',rgb('#7D98AD'))
        c.set_color('neighbour_warm',rgb('#C5A16C'))
        c.viewport(2200,2200)
        for site in [72,109]:
            left=f'NHE1_M{site}';right=f'IQCG_M{site}'
            face_sidechain(c,left,site)
            c.zoom(f'({left} or {right}) and chain A and resi {site}',buffer=4.8,complete=1)
            view=c.get_view()
            for label,obj in [('NHE1',left),('IQCG',right)]:
                c.disable('all');c.enable(obj);c.set_view(view)
                c.set('cartoon_transparency',.25,obj+' and chain A')
                c.set('cartoon_transparency',.42,obj+' and chain P')
                c.color('neighbour_bluegrey',obj+' and chain A and elem C')
                c.color('neighbour_warm',obj+' and chain P and elem C')
                c.set('stick_transparency',.18,obj)
                c.set('stick_radius',.11,obj)
                focus=f'{obj} and chain A and resi {site}'
                c.hide('spheres',obj)
                c.color('focus_deep',focus+' and elem C')
                c.color('sulfur',focus+' and elem S')
                c.set('stick_transparency',0,focus)
                c.set('stick_radius',.24,focus)
                c.clip('slab',150)
                name=f'detail_focus_M{site}_{label}'
                marks={f'M{site}':focus+' and name CE',f'M{site}_CA':focus+' and name CA',
                       f'M{site}_SD':focus+' and name SD'}
                manifest['assets'][name]=render(c,name,2200,2200,marks)
                manifest['assets'][name]['framing']='All continuous peptide cartoon selections retained; tight camera crop naturally continues through image border. Approximate 15 A local field of view.'
                manifest['assets'][name]['focus_carbon_color']='#276F6F'
                # Exact white-background QA: raw transparent PNG must alpha-
                # composite over white before assessing the apparent opacity.
                im=Image.open(OUT/(name+'.png')).convert('RGBA')
                Image.alpha_composite(Image.new('RGBA',im.size,'white'),im).convert('RGB').save(OUT/('QA_'+name+'_white.png'))
            paired_trim(manifest['assets'],[f'detail_focus_M{site}_NHE1',f'detail_focus_M{site}_IQCG'])
        c.save(str(OUT/'CaM_continuous_context_v80.pse'))
    (OUT/'render_manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
