"""Exact A003R static input: locate Q41/Q2 without substituting Q62 frames."""
from render_structures_v80 import *

def main():
    src=Path('E:/Project/FlexPPI/experiments/_case_validation_sources/ubiquitin/inputs/A003R_input_structure.pdb')
    shutil.copy2(src,OUT/'inputs'/src.name)
    with pymol2.PyMOL() as pm:
        c=pm.cmd;c.load(str(src),'ubiquitin_UBA');c.remove('hydro')
        c.bg_color('white')
        for k,v in dict(ray_opaque_background=0,opaque_background=0,ray_trace_mode=0,ray_shadows=0,
                        antialias=3,orthoscopic=1,depth_cue=0,specular=.12,ambient=.65,direct=.35,
                        cartoon_fancy_helices=1,cartoon_fancy_sheets=1,cartoon_flat_sheets=1,
                        cartoon_smooth_loops=1,cartoon_sampling=20,max_threads=4,valence=0).items():c.set(k,v)
        for n,col in COLORS.items():c.set_color(n,rgb(col))
        c.hide('everything');c.dss();c.show('cartoon','all')
        c.color('N_lobe','chain A');c.set('cartoon_color','N_lobe','chain A')
        c.color('context','chain B');c.set('cartoon_color','context','chain B')
        c.set('cartoon_transparency',.50,'chain A')
        for site,col in [(41,'partner'),(2,'key')]:
            focus=f'chain A and resi {site}'
            c.show('sticks',focus+' and not name N+C+O')
            c.color(col,focus+' and elem C');c.color('oxygen',focus+' and elem O');c.color('nitrogen',focus+' and elem N')
            c.set('stick_radius',.25,focus)
        c.viewport(3400,3400);c.orient('all');c.turn('x',28);c.turn('y',32);c.turn('z',-20);c.zoom('all',buffer=1.3,complete=1)
        marks={'Q41':'chain A and resi 41 and name CD','Q2':'chain A and resi 2 and name CD',
               'ubiquitin':'chain A and name CA','UBA':'chain B and name CA'}
        e=render(c,'ubiquitin_Q41_Q2_locator',3400,3400,marks)
        paired_trim({'ubiquitin_Q41_Q2_locator':e},['ubiquitin_Q41_Q2_locator'])
        c.save(str(OUT/'ubiquitin_Q41_Q2_locator.pse'))
        out={'asset':e,'input_source':str(src),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
             'scope':'Static exact A003R model input. Complete Ub72/UBA47 cartoon chains. Q41/Q2 are actual residues, no MD frame substitutions or rebuilt atoms.'}
    (OUT/'ubiquitin_Q41_Q2_locator.json').write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
