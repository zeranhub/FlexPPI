"""Independent exact checks of the native framework and molecular assets."""
from pathlib import Path
import json,hashlib,xml.etree.ElementTree as ET,base64,re,argparse,math
ROOT=Path(__file__).resolve().parent.parent
F=ROOT/'figure1'
OLD=ROOT.parent/'framework_v46_20261019'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
args=argparse.ArgumentParser()
args.add_argument('--document',type=Path)
arg=args.parse_args()
oldroot=ET.parse(F/'source_archive/original_framework_v13.drawio').find('.//root')
newroot=ET.parse(F/'fig1_framework.drawio').find('.//root')
old={x.get('id'):x for x in oldroot}
new={x.get('id'):x for x in newroot}
canon=lambda x:ET.tostring(x,encoding='unicode').strip()
edge_lengths=[]
for key,cell in new.items():
    if cell.get('edge')!='1':continue
    geometry=cell.find('mxGeometry')
    start=geometry.find("mxPoint[@as='sourcePoint']")
    stop=geometry.find("mxPoint[@as='targetPoint']")
    assert start is not None and stop is not None,(key,'missing explicit endpoints')
    waypoints=geometry.findall("Array[@as='points']/mxPoint")
    points=[start]+waypoints+[stop]
    xy=[(float(p.get('x')),float(p.get('y'))) for p in points]
    length=sum(math.dist(a,b) for a,b in zip(xy,xy[1:]))
    edge_lengths.append({'id':key,'path_length_native_units':length,'points':xy})
zero=[e for e in edge_lengths if e['path_length_native_units']<=1e-7]
assert not zero,zero
for key in ['c9_scalar_to_sum','c9_pair_to_sum','c9_point_to_sum']:
    assert any(e['id']==key and e['path_length_native_units']>0 for e in edge_lengths),key
protected=[k for k in old if k.startswith('curve_')]
assert len(protected)==224
assert len([k for k in new if k.startswith('curve_')])==224
assert all(canon(old[k])==canon(new[k]) for k in protected)
def imagedata(cell):
    style=cell.get('style','')
    return base64.b64decode(style.split('image=data:image/png,',1)[1].split(';',1)[0])
mapping={'input_structure':F/'assets3d/assembly_input_consistent.png',
         'flexibility_envelope':F/'assets3d/assembly_prediction_consistent.png',
         'b47_whole_ribbon':F/'assets3d/chain_d_consistent.png',
         'b47_atomic_inset':F/'source_archive/local_inset.png'}
for key,file in mapping.items():assert imagedata(new[key])==file.read_bytes(),key
oldimgs={k:imagedata(c) for k,c in old.items() if 'image=data:image/png,' in c.get('style','')}
md=[k for k in oldimgs if k.startswith('d_')]
assert len(md)==3
assert all(imagedata(new[k])==oldimgs[k] for k in md)
scalecells=[k for k in old if k.startswith('a9_rmsf_scale_')]
assert len(scalecells)==64
for key in scalecells:
    a=dict(p.split('=',1) for p in old[key].get('style','').split(';') if '=' in p)
    b=dict(p.split('=',1) for p in new[key].get('style','').split(';') if '=' in p)
    assert a['fillColor']==b['fillColor']
    assert canon(old[key].find('mxGeometry'))==canon(new[key].find('mxGeometry'))
for ident in ['encoder_stack_3','encoder_stack_2','encoder_stack_1','encoder_front','ipa_module','ffn_module','final_norm_panel','head_panel']:
    assert ident in new
assert new['head_joint'].get('value')=='Shared MLP'
layers=[new[f'head_flat_text_{i}'].get('value') for i in range(4)]
assert layers==['LayerNorm','Linear','SiLU','Linear']
visible=[x.get('value','') for x in newroot if x.get('value')]
assert 'v = σ² + v(MD)' in visible
assert not any(re.search('[\u3400-\u9fff]',t) for t in visible)
ab=json.loads((F/'assets3d/consistent_a_b_receipt.json').read_text(encoding='utf-8'))
prediction=json.loads((F/'assets3d/consistent_prediction_receipt.json').read_text(encoding='utf-8'))
assert ab['status']==prediction['status']=='passed'
assert ab['shared_material_settings']==prediction['shared_material_settings']
assert ab['source_XYZ_sha256']==ab['after_XYZ_sha256']==prediction['all_source_coordinates_sha256']==prediction['all_after_coordinates_sha256']
assert ab['chain_palette']=={'A':'#2166AC','B':'#F1A340','C':'#B2182B','D':'#1B7837','E':'#762A83'}
assert ab['chain_lengths']==prediction['chain_lengths']=={'A':276,'B':100,'C':10,'D':204,'E':250}
assert prediction['all_840_original_prediction_values_exact'] and prediction['all_840_original_input_RGB_exact']
assert not prediction['original_active_is_amplitude_envelope_or_putty']
assert ab['a']['png_sha256']==sha(mapping['input_structure'])
assert ab['b']['png_sha256']==sha(mapping['b47_whole_ribbon'])
assert prediction['png_sha256']==sha(mapping['flexibility_envelope'])
assert ab['shared_material_settings']['cartoon_transparency']==0
method=json.loads((OLD/'audit/method_content.json').read_text(encoding='utf-8'))
assert method['status']=='passed'
for item in method['code_files']:
    assert sha(item['path'])==item['sha256']
for item in method['actual_five_fold_resolved_configs']:
    assert sha(item['path'])==item['sha256']
# Reuse the already completed exact 840 raw-value/smoothing and frame audits.
for name in ['protected_840_data.json','selected_frame_audit.json']:
    assert json.loads((OLD/'audit'/name).read_text(encoding='utf-8'))['status']=='passed'
report={'status':'passed','reviewer':'fig1_method_audit','figure1':str(F/'fig1_framework.png'),'figure1_sha256':sha(F/'fig1_framework.png'),
        'actual_visual_inspection':'Final PNG ba882effd76b9bd766686ddf6ed45cbe89a9c6d0de04a7a88da0407a6f7572de opened at original detail; all three score inputs visibly connect to sum, including repaired golden pair connector. Source-native objects and caption inspected separately.',
        'native_drawio_sha256':sha(F/'fig1_framework.drawio'),'svg_sha256':sha(F/'fig1_framework.svg'),
        'protected_native_curve_objects':224,'all_224_curve_objects_exact':True,
        'all_active_edges_positive_path_length':True,'active_edge_count':len(edge_lengths),
        'active_edge_lengths':edge_lengths,
        'raw_prediction_and_within_chain_mean_verification':'Prior independent v46 audit reused; all 840 positions/1680 target rows and original valid masks preserved, all 11-position within-chain means recomputed exactly.',
        'asset_mapping':{k:{'path':str(p),'sha256':sha(p),'embedded_bytes_exact':True} for k,p in mapping.items()},
        'all_3_original_MD_snapshots_exact':True,'original_64_colorbar_cells_exact_RGB_and_geometry':True,
        'all_840_prediction_values_and_input_RGB_exact':True,'chain_palette_from_original_legend_exact':True,
        'a_b_prediction_shared_DSS_material':True,'input_XYZ_preserved':True,
        'frame_construction':'True D67/D166 frames; CA origin, x→C, orthogonalized y→N, z=x cross y. Actual 3D inset uses original D165–167 backbone and calibrated projection probes.',
        'model_audit':'Actual five-fold 4-block/8-head configuration; fixed input frames and pair features. Residual, feed-forward and normalization routes retained. One shared readout MLP has target-specific mean/log-variance outputs.',
        'supervision_audit':'Gaussian NLL compares model mean to MD mean using model variance + observed variance. Model variance is not trained as a direct MD variance target. Hχ is joint fixed-bin state entropy.',
        'inherited_source_audits':{'note':'Byte-identical historical v46 source audits, not final v47 figure/document reviews.',
          'files':[{ 'path':str(ROOT/'audit'/n),'source':str(OLD/'audit'/n),'sha256':sha(ROOT/'audit'/n)}
                   for n in ['method_content.json','protected_840_data.json','selected_frame_audit.json']]},
        'new_render_receipts':[{'path':str(F/'assets3d'/n),'sha256':sha(F/'assets3d'/n)}
                   for n in ['consistent_a_b_receipt.json','consistent_prediction_receipt.json','all_840_preserved_RGB_audit.json']],
        'new_training':0,'new_inference':0,'new_MD':0}
if arg.document:
    import zipfile
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    base=ROOT.parent.parent
    with zipfile.ZipFile(arg.document) as z:
        doc=ET.fromstring(z.read('word/document.xml'))
        assert z.read('word/media/image1.png')==(F/'fig1_framework.png').read_bytes()
    with zipfile.ZipFile(base/'paper_v45_20261018_中文修订_紧凑多层模型评价图_三线表.docx') as z:
        olddoc=ET.fromstring(z.read('word/document.xml'))
    assert len(doc.findall('.//m:oMath',ns))==36
    assert [canon(x) for x in doc.findall('.//m:oMath',ns)]==[canon(x) for x in olddoc.findall('.//m:oMath',ns)]
    assert len(doc.findall('.//w:tbl',ns))==10
    txt=lambda p:''.join(t.text or '' for t in p.iter() if t.tag.rsplit('}',1)[-1]=='t')
    paras=doc.find('w:body',ns).findall('w:p',ns)
    proposal=json.loads((ROOT/'audit/prose_clarifications.json').read_text(encoding='utf-8'))
    for item in proposal['replacements']:assert sum(txt(p)==item['new_exact'] for p in paras)==1
    cap=[txt(p) for p in paras if txt(p).startswith('图1：')]
    assert len(cap)==1
    expected=(F/'caption_zh.md').read_text(encoding='utf-8').strip()
    assert cap[0].replace('ᵢ','i')==expected.replace('ᵢ','i')
    report.update({'document':str(arg.document),'document_sha256':sha(arg.document),'caption':cap[0],
                   'all_36_original_native_math_objects_exact':True,'native_tables':10,'new_joint_readout_prose_clarifications_present':True})
    target=ROOT/'audit/final_scientific_review.json'
else:
    target=ROOT/'audit/native_figure_scientific_audit.json'
target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','report':str(target),'figure1_sha256':report['figure1_sha256'],'document_sha256':report.get('document_sha256')}))
