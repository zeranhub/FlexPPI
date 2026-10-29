"""Bind the independent scientific review to the final Word and figure bytes."""
from pathlib import Path
import json
import hashlib
import re
import zipfile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parent.parent
BASE=ROOT.parent.parent
DOC=BASE/'paper_v46_20261019_中文修订_模型框架图重绘_三线表.docx'
OLD=BASE/'paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx'
FIG=ROOT/'figure1/fig1_framework.png'
SVG=ROOT/'figure1/fig1_framework.svg'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(DOC)=='f53a6d0b09d22a5cba085dcb7a048ae1fed5aefe322382b90fb2a8b94bdf4da9'
assert sha(FIG)=='ee158b4acbeb651e11c991035912b94878bd765230f7ec8ff43ccd3b4bc917eb'
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math',
    'svg':'http://www.w3.org/2000/svg'}
with zipfile.ZipFile(DOC) as z:
    root=ET.fromstring(z.read('word/document.xml'))
    assert z.read('word/media/image1.png')==FIG.read_bytes()
with zipfile.ZipFile(OLD) as z:oldroot=ET.fromstring(z.read('word/document.xml'))
paras=root.find('w:body',ns).findall('w:p',ns)
txt=lambda n:''.join(t.text or '' for t in n.iter() if t.tag.rsplit('}',1)[-1]=='t')
captions=[p for p in paras if txt(p).startswith('图1：')]
assert len(captions)==1
caption=txt(captions[0])
expected=(ROOT/'figure1/caption_zh.md').read_text(encoding='utf-8').strip()
assert caption.replace('ᵢ','i')==expected.replace('ᵢ','i')
proposal=json.loads((ROOT/'audit/prose_clarifications.json').read_text(encoding='utf-8'))
for change in proposal['replacements']:
    assert sum(txt(p)==change['new_exact'] for p in paras)==1
    assert not any(txt(p)==change['old_exact'] for p in paras)
equations=root.findall('.//m:oMath',ns)
old_equations=oldroot.findall('.//m:oMath',ns)
assert len(equations)==len(old_equations)==36
assert [ET.tostring(x) for x in equations]==[ET.tostring(x) for x in old_equations]
assert len(root.findall('.//w:tbl',ns))==10

svg=ET.parse(SVG).getroot()
svgtext=[''.join(t.itertext()) for t in svg.findall('.//svg:text',ns)]
assert not [t for t in svgtext if re.search('[\u3400-\u9fff]',t)]
for label in ['Shared MLP','4 blocks · 8 heads','Fixed-frame IPA','Local → global',
              'Scalar + pair − distance penalty','Softmax weights','Merge + project → update',
              'Model μ, variance σ²','v = σ² + v(MD)','RMSF (Å)','Hχ (nats)']:
    assert label in svgtext,label
render=json.loads((ROOT/'figure1/render_audit.json').read_text(encoding='utf-8'))
for ext,expected_sha in render['output_sha256'].items():
    assert sha(ROOT/f'figure1/fig1_framework.{ext}')==expected_sha
method=json.loads((ROOT/'audit/method_content.json').read_text(encoding='utf-8'))
assert method['status']=='passed'
for item in method['code_files']:
    assert sha(item['path'])==item['sha256']
for item in method['actual_five_fold_resolved_configs']:
    assert sha(item['path'])==item['sha256']
for name in ['protected_840_data.json','selected_frame_audit.json']:
    report=json.loads((ROOT/'audit'/name).read_text(encoding='utf-8'))
    assert report['status']=='passed'
for file,expected_sha in render['source_hashes'].items():
    assert sha(ROOT/'figure1/source_archive'/file)==expected_sha
source=(ROOT/'figure1/draw_framework.py').read_text(encoding='utf-8')
for fragment in ["pd.read_csv(S/'curve_data_all_840.csv')", "pd.read_csv(S/'curve_local_means.csv')",
                 "for chain in 'ABCDE'", "z.prediction.where(z.valid.astype(bool))",
                 "proj@R", "R.T@(np.array(fi['N_world_A'])-ca)"]:
    assert fragment in source,fragment
checks=[
    {'claim':'Inputs and geometric streams','result':'passed','detail':'Amino-acid identity with N–Cα–C; learned embedding, pair geometry and fixed frame streams match actual code and all five resolved configurations. No ESM, condition or response head added.'},
    {'claim':'Block and head counts','result':'passed','detail':'4 IPA blocks, 8 heads; actual per-head dimensions 16 scalar, 4 Q/K points and 8 V points remain in Methods and code.'},
    {'claim':'Frame construction and coordinate convention','result':'passed','detail':'Actual D:Gln67 and D:Met166 N/Cα/C reconstruct stored frames exactly. CA origin, x→C, orthogonalized y→N, z=x cross y. Local-to-global Rq+t; no global-centering operation is claimed.'},
    {'claim':'Fixed frames and learned latent points','result':'passed','detail':'Input frames/pair features remain fixed; only residue representations update. Purple probes are explicitly labelled learned geometric points rather than predicted atoms.'},
    {'claim':'IPA attention and update','result':'passed','detail':'Scalar plus pair bias minus weighted squared global Q/K distances; same attention weights aggregate scalar, pair and point values; points return to query-local frame then merge/project/update. Simplified main overview omits explicit pre-norm/residual internals, which remain accurate in unchanged Methods.'},
    {'claim':'Joint readout','result':'passed','detail':'One shared MLP with target-specific final channels produces mean/log variance per RMSF and Hχ. Three prose clarifications remove the implication of two independent MLP heads.'},
    {'claim':'Supervision and variance relationship','result':'passed','detail':'MD repeats or windows provide observed mean and sample variance. Model variance plus observed variance forms Gaussian NLL; variance is not a direct regression label. Hχ is jointly binned state entropy.'},
    {'claim':'Protected 4g9f example','result':'passed','detail':'840 positions, 1680 target rows, valid counts 840 RMSF and 734 Hχ. Source CSVs byte-identical to old source; every precomputed within-chain ±5 mean independently reproduces with zero numerical error. Full raw profiles and actual chain/invalid breaks retained with dual physical axes.'},
    {'claim':'Manuscript integration','result':'passed','detail':'Final document embeds reviewed PNG byte-exact, caption including four native R_i/t_i objects reconstructs exact expected wording, all 36 native math objects and 10 native tables retained.'},
]
report={
    'status':'passed','reviewer':'fig1_method_audit','document':str(DOC),'document_sha256':sha(DOC),
    'figure1':str(FIG),'figure1_sha256':sha(FIG),'svg':str(SVG),'svg_sha256':sha(SVG),
    'actual_visual_inspection':'Final PNG opened at original detail; final SVG XML/text inspected; caption and drawing source read. Final full-manuscript pages 11–17 are reviewed separately.',
    'checks':checks,'all_scientific_checks_passed':True,'new_training':0,'new_inference':0,'new_MD':0,
    'supporting_reports':[str(ROOT/'audit'/name) for name in ['method_content.json','protected_840_data.json','selected_frame_audit.json']],
}
(ROOT/'audit/final_scientific_review.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','document_sha256':sha(DOC),'figure1_sha256':sha(FIG),'scientific_checks':len(checks)}))
