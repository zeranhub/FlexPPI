from pathlib import Path
import importlib.util,json
P=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('render_base',P/'render_mechanism_structures.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
manifest=json.loads((P/'render_manifest.json').read_text())
with m.pymol2.PyMOL() as p:
 c=p.cmd
 for code,obj,ch in [('1A06','kinase','A'),('1MXE','bound','B')]:
  c.reinitialize();c.load(str(P/(code+'_overview.pse')));c.disable('all');c.enable(obj)
  if code=='1A06':c.set('cartoon_transparency',.30,'kinase and resi 10-293')
  else:
   ss='bound and chain A and resi 116+120+124+127 and not name N+C+O';c.show('sticks',ss);c.set('stick_radius',.12,ss);c.color('contact',ss)
  manifest['assets'][code+'_context']=m.render(c,code+'_context',obj,ch)
(P/'render_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')