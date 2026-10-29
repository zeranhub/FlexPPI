"""Recompose the frozen v94 evidence; only layout and panel references change."""
from pathlib import Path
import shutil
import json

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'revision_v94_restore_per_chi'
for folder in ['figure3/source_data', 'render', 'manuscript', 'audit', 'pymol/inputs']:
    (ROOT / folder).mkdir(parents=True, exist_ok=True)
for p in (OLD / 'figure3/source_data').glob('*'):
    if p.is_file(): shutil.copy2(p, ROOT / 'figure3/source_data' / p.name)
for name in ['render_docx_word_wrapper.py', 'export_word_readonly.ps1']:
    shutil.copy2(OLD / 'render' / name, ROOT / 'render' / name)
manifest = json.loads((OLD / 'pymol/render_manifest.json').read_text(encoding='utf-8'))
manifest['assets'] = {k: v for k, v in manifest['assets'].items() if k in ['cam_overview', 'ub_overview']}
for asset in manifest['assets'].values():
    for key in ['png', 'pse']:
        shutil.copy2(OLD / 'pymol' / asset[key], ROOT / 'pymol' / asset[key])
for source in manifest['inputs'].values():
    shutil.copy2(OLD / 'pymol' / source['pdb'], ROOT / 'pymol' / source['pdb'])
(ROOT / 'pymol/render_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')

s = (OLD / 'figure3/compose_fig3_v94.py').read_text(encoding='utf-8')
def change(old, new):
    global s
    assert old in s, old
    s = s.replace(old, new)

start = s.index('# a/b: only backbone')
stop = s.index("x=cam.CaMKI_native_position", start)
s = s[:start] + '''# a/b: the two intact ribbon structures are grouped across the top row.
panel('a','CaM–CaMKI',1,2)
molecule('cam_overview',10,9,58,50.75,{'W303':(.24,.61),'M316':(.88,.35)})
line(18,61,23,61,color=PURPLE,lw=2.6);txt(25,59.4,'CaMKI',9)
line(49,61,54,61,color=GRAY,lw=2.6);txt(56,59.4,'CaM',9)
panel('b','Ubiquitin–UBA',99,2)
molecule('ub_overview',108,9,58,50.75,{'Q41':(.56,.91),'Q62':(.92,.10)})
line(112,61,117,61,color=PURPLE,lw=2.6);txt(119,59.4,'Ubiquitin',9)
line(150,61,155,61,color=GRAY,lw=2.6);txt(157,59.4,'UBA',9)

# c: all 25 CaMKI positions; use the full width for sequence and contacts.
panel('c','CaMKI · backbone RMSF',1,69)
for xx,label,col,mark,ls in [(116,'MD',INK,'o',':'),(136,'Full',BLUE,'o','-'),(157,'Intra',WARM,'s','--')]:
    line(xx,77.5,xx+4.5,77.5,color=col,lw=1.15,marker=mark,ms=2.8,ls=ls);txt(xx+6,75.9,label,9)
''' + s[stop:]
change("txt(75+(pos-lo)/25*105,15.1", "txt(16+(pos-lo)/25*164,82.1")
change("a=ax(75,20,105,27)", "a=ax(16,87,164,24)")
change("add_values('b','RMSF_A'", "add_values('c','RMSF_A'")
change("a=ax(75,53,105,7)", "a=ax(16,117,164,7)")
change("add_values('b','CaM_contact_residues_4p5A'", "add_values('c','CaM_contact_residues_4p5A'")

start=s.index('# c: all 72 positions')
stop=s.index("a=ax(129,135,50,28)", start)
old_profile=s[start:s.index('# d/e/f:', start)]
old_profile=old_profile.replace("# c: all 72 positions", "# d: all 72 positions")
old_profile=old_profile.replace("panel('c','Ubiquitin · complete side-chain profile',1,74)", "panel('d','Ubiquitin · side-chain Hχ',1,134)")
old_profile=old_profile.replace("[(138,'MD'", "[(50,'MD'").replace("(158,'Full'", "(76,'Full'")
old_profile=old_profile.replace("xx,77,xx+4.5,77", "xx,142,xx+4.5,142").replace("xx+6,75.4", "xx+6,140.4")
old_profile=old_profile.replace("a=ax(16,85,164,25)", "a=ax(16,148,96,23)")
old_profile=old_profile.replace("add_values('c','Hchi_nats'", "add_values('d','Hchi_nats'")
new_gln='''# e: all six same-type residues, with an explicit MD/Full legend.
panel('e','All six Gln positions',118,134)
for xx,label,col,mark,filled in [(130,'MD',INK,'o',False),(155,'Full',BLUE,'s',True)]:
    line(xx,142,xx+2,142,color=col,lw=0,marker=mark,ms=3,mfc=col if filled else 'white',mec=col);txt(xx+4,140.4,label,9)
'''
s=s[:start]+old_profile+new_gln+s[stop:]
change("a=ax(129,135,50,28)", "a=ax(130,148,49,23)")
change("marker=mark,ms=3,mfc=col", "marker=mark,ms=3,markevery=[0],mfc=col")
change("marker=marker,ms=3.6,lw=0", "marker=marker,ms=3.6,markevery=[0],lw=0")
change("marker=mark,ms=4,mfc=BLUE", "marker=mark,ms=4,markevery=[0],mfc=BLUE")
change("{'b':100,'c':132,'e':12,'f':302,'g':332}", "{'c':100,'d':132,'e':12,'f':302,'g':332}")
change("visual_review='pending'", "structure_panels_grouped_at_top=True,local_zoom_panels_removed=True,visual_review='pending'")
confirmation='''# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a,b) native structure overview -> inherited v94 PyMOL overview renders -> native run
# (c,d) absolute residue profiles -> inherited v94 native profiles -> param inherit
# (e,f,g) complete-point scatter -> inherited v94 native plots -> param inherit
# RULE: native run = original PyMOL PNG via Image.open().ax.imshow().
#       param inherit = unchanged data and plotting semantics; layout only.
'''
(ROOT / 'figure3/compose_fig3_v95.py').write_text(confirmation+s, encoding='utf-8')
print('Prepared v95 composition and unchanged source data.')
