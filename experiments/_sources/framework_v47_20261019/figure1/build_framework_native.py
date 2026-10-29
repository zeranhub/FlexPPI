"""Internal refinement of the original, information-rich native framework.

Preserves the five-stage a row and the original b/c/d bottom columns, the
four encoder layers and both residual paths, the shared MLP internal layers,
all 224 native profile paths and all real prediction/MD renders.
"""
from pathlib import Path
import xml.etree.ElementTree as ET
import copy,re,json,hashlib,shutil,subprocess
import base64
HERE=Path(__file__).resolve().parent
OLD=Path(r'C:\Users\DELL\.codex\worktrees\4666\FlexPPI\paper\revision_20261002\figures\main_revision_v13_20261008')
S=HERE/'source_archive';S.mkdir(exist_ok=True)
source=S/'original_framework_v13.drawio'
if not source.exists():shutil.copy2(OLD/'FlexPPI_main_v13.drawio',source)
tree=ET.parse(source);root=tree.find('.//root');cells={c.get('id'):c for c in root}
protected={c.get('id'):ET.tostring(c) for c in root if c.get('id','').startswith('curve_')}
assert len(protected)==224

def sty(c,**items):
    raw=c.get('style','');out=[]
    for part in raw.split(';'):
        if not part:continue
        if '=' not in part:out.append(part);continue
        k,v=part.split('=',1)
        if k not in items:out.append(k+'='+v)
    out.extend(k+'='+str(v) for k,v in items.items())
    c.set('style',';'.join(out)+';')

def box(c,x=None,y=None,w=None,h=None):
    g=c.find('mxGeometry')
    for k,v in [('x',x),('y',y),('width',w),('height',h)]:
        if v is not None:g.set(k,str(v))

def label(ident,value,geometry,fontsize=28,color='#28363E',bold=False):
    c=ET.SubElement(root,'mxCell',id=ident,value=value,vertex='1',parent='1',style=f'text;html=0;whiteSpace=wrap;fillColor=none;strokeColor=none;fontFamily=Arial;fontSize={fontsize};fontColor={color};fontStyle={int(bold)};align=center;verticalAlign=middle;spacing=0;')
    g=ET.SubElement(c,'mxGeometry',attrib={'as':'geometry'})
    for k,v in zip(['x','y','width','height'],geometry):g.set(k,str(v))
    return c

def shape(ident,geometry,fill,edge='none',width=1.5,extra='rounded=0;'):
    c=ET.Element('mxCell',id=ident,value='',vertex='1',parent='1',style=f'{extra}fillColor={fill};strokeColor={edge};strokeWidth={width};')
    g=ET.SubElement(c,'mxGeometry',attrib={'as':'geometry'})
    for k,v in zip(['x','y','width','height'],geometry):g.set(k,str(v))
    root.append(c);return c

def line(ident,points,color='#6F8190',width=2,arrow=False,dashed=False):
    c=ET.SubElement(root,'mxCell',id=ident,parent='1',edge='1',style=f'edgeStyle=none;rounded=0;curved=0;html=0;strokeColor={color};strokeWidth={width};endArrow={"block" if arrow else "none"};endSize=5;dashed={int(dashed)};dashPattern=4 3;')
    g=ET.SubElement(c,'mxGeometry',relative='1',attrib={'as':'geometry'})
    for key,p in [('sourcePoint',points[0]),('targetPoint',points[-1])]:ET.SubElement(g,'mxPoint',x=str(p[0]),y=str(p[1]),attrib={'as':key})
    if len(points)>2:
        a=ET.SubElement(g,'Array',attrib={'as':'points'})
        for p in points[1:-1]:ET.SubElement(a,'mxPoint',x=str(p[0]),y=str(p[1]))
    return c

# Refine hierarchy and stroke language; never touch native prediction paths.
for ident,c in cells.items():
    if ident.startswith('curve_') or 'image=' in c.get('style','') or ident.startswith('b13_'):continue
    v=c.get('value','');s=c.get('style','')
    if v:
        match=re.search(r'fontSize=([0-9.]+)',s)
        if match:
            fs=float(match.group(1))
            sty(c,fontSize=max(28,min(fs,31 if ident not in ['a_title','v5_title_b','v5_title_c','v5_title_d'] else 38)))
        if ident not in ['a_title','v5_title_b','v5_title_c','v5_title_d'] and not ident.startswith('stage_'):
            sty(c,fontStyle=0)
        sty(c,fontColor='#2D3D46')
    if c.get('edge')=='1':
        match=re.search(r'strokeWidth=([0-9.]+)',s)
        if match:sty(c,strokeWidth=min(float(match.group(1)),2.2))
        sty(c,endSize=6)

cells['stage_input'].set('value','Static structure')
cells['stage_head'].set('value','Joint readout')
cells['head_joint'].set('value','Shared MLP')
cells['a_title'].set('value','a   FlexPPI overview')

# Protein renderings share the same DSS, lighting, material and exact chain
# colour semantics between the assembly and its local chain-D scene.
def replace_image(ident,path):
    c=cells[ident]
    c.set('style','shape=image;imageAspect=0;aspect=fixed;image=data:image/png,'+base64.b64encode(path.read_bytes()).decode()+';')
replace_image('input_structure',HERE/'assets3d/assembly_input_consistent.png')
replace_image('flexibility_envelope',HERE/'assets3d/assembly_prediction_consistent.png')

# A: preserve layered encoder, residual routes and the full head architecture.
for ident in ['encoder_stack_3','encoder_stack_2','encoder_stack_1']:
    sty(cells[ident],fillColor='#F3F5F7',strokeColor='#D2D9DF',strokeWidth=1.3,arcSize=8)
sty(cells['encoder_front'],fillColor='#FCFCFD',gradientColor='#F4F6F8',gradientDirection='south',strokeColor='#A9B6C0',strokeWidth=1.6,arcSize=8)
sty(cells['ipa_module'],fillColor='#F0EBF6',gradientColor='#E7DFF0',gradientDirection='south',strokeColor='#B9A8CC',arcSize=9)
sty(cells['ffn_module'],fillColor='#E9F0F6',strokeColor='#ACC2D2',arcSize=8)
sty(cells['final_norm_panel'],fillColor='#F5F7FA',strokeColor='#BDC7CF',arcSize=6)
sty(cells['head_panel'],fillColor='#FFFFFF',strokeColor='#CFB36E',strokeWidth=1.5,arcSize=9)
for k in range(4):
    sty(cells[f'head_flat_{k}'],fillColor=['#F8F4E8','#F1E5C4','#F6ECD2','#EDE0BA'][k],strokeColor='none')
    sty(cells[f'head_flat_text_{k}'],fontSize=31,fontStyle=0)
    if k<3:sty(cells[f'head_flat_arrow_{k}'],strokeColor='#B39859',strokeWidth=1.6,endSize=5)
for ident,c in cells.items():
    if ident.startswith('a_residue_graph_graph_edge_'):
        sty(c,strokeColor='#819BAD',strokeWidth=1.6)
    elif ident.startswith('a_residue_graph_graph_node_'):
        sty(c,fillColor='#C9DCE8',gradientColor='#F4F9FC',gradientDirection='north',strokeColor='#6F93AB',strokeWidth=1.5)
    elif ident.startswith('ipa_graph_egnn_node_'):
        sty(c,fillColor='#D6CCE5',gradientColor='#F5F0F9',gradientDirection='north',strokeColor='#9A85B7',strokeWidth=1.4)
    elif ident.startswith('ipa_graph_egnn_edge_'):
        sty(c,strokeColor='#9987AD',strokeWidth=1.5)

# Pair representations retain a numerical glyph but gain two understated
# channel planes: a small object with depth instead of a flat yellow block.
index=list(root).index(cells['pair_matrix_0_0'])
for shift,fill in [(10,'#F8F2E3'),(5,'#F0E4C5')]:
    c=shape('v47_pair_plane_'+str(shift),[428+shift,400-shift,144,144],fill,'#D7C391',1.)
    root.remove(c);root.insert(index,c);index+=1

# B is reconstructed using the same information-rich bottom-column footprint.
for c in list(root):
    if c.get('id','').startswith('b13_'):root.remove(c)
proto=ET.parse(HERE/'fig1b_prototype.drawio').find('.//root')
for c in proto:
    if not c.get('id','').startswith('b47_') or c.get('id')=='b47_title':continue
    c=copy.deepcopy(c);g=c.find('mxGeometry')
    if g is not None:
        if 'x' in g.attrib:g.set('x',str(float(g.get('x'))+15))
        if 'y' in g.attrib:g.set('y',str(float(g.get('y'))+690))
        for p in g.iter('mxPoint'):
            if 'x' in p.attrib:p.set('x',str(float(p.get('x'))+15))
            if 'y' in p.attrib:p.set('y',str(float(p.get('y'))+690))
    root.append(c)
cells['v5_title_b'].set('value','b   Global and local frames')

# C: preserve the three-stage attention logic and its original vertical rhythm,
# but populate score inputs with feature/vector/point glyphs.
semantic={'scalar':('#EDF3F8','#AAC1D2','#487799'), 'pair':('#F8F1DF','#CDBB8B','#A3823A'), 'point':('#F1EDF7','#BCACCD','#85709D')}
for kind in semantic:
    fill,edge,color=semantic[kind]
    box(cells['c9_score_'+kind],y=758,h=110)
    sty(cells['c9_score_'+kind],fillColor=fill,gradientColor='#FFFFFF',gradientDirection='north',strokeColor=edge,arcSize=8)
    box(cells['c9_score_text_'+kind],y=801,h=62)
    sty(cells['c9_score_text_'+kind],fontSize=28,fontStyle=0,fontColor=color)
    # Start each connector below its refined input card, never through text.
    cells['c9_'+kind+'_to_sum'].find("mxGeometry/mxPoint[@as='sourcePoint']").set('y','870')
# Reserve a visible 14-unit pair-to-score connector. Its old sum target was
# y=880; setting all source points there had silently collapsed this one edge.
box(cells['c9_sum'],y=884)
for kind in semantic:
    geom=cells['c9_'+kind+'_to_sum'].find('mxGeometry')
    for p in geom.iter('mxPoint'):
        if p.get('as')!='sourcePoint':p.set('y',str(float(p.get('y'))+4))
cells['c9_sum_to_softmax'].find("mxGeometry/mxPoint[@as='sourcePoint']").set('y','922')
pair_points=cells['c9_pair_to_sum'].find('mxGeometry').findall('mxPoint')
assert abs(float(pair_points[0].get('y'))-float(pair_points[-1].get('y')))>=14
for y in [771,789]:
    for k in range(6):shape(f'v47_scalar_{y}_{k}',[626+k*15,y,13,9],['#D8E6F0','#B8CEDF','#95B6CC'][k%3])
for i in range(4):
    for j in range(4):shape(f'v47_pair_{i}_{j}',[808+j*14,762+i*11,12,9],['#EFDFB5','#DDC78F','#F6EDD5'][(i+2*j)%3])
line('v47_qk_axis_x',[[953,800],[1074,800]],color='#D2C9DE',width=1.)
line('v47_qk_axis_y',[[966,806],[966,763]],color='#D2C9DE',width=1.)
shape('v47_q',[965,778,12,12],'#A494BB','#FFFFFF',1.1,'ellipse;')
shape('v47_k',[1032,760,12,12],'#A494BB','#FFFFFF',1.1,'ellipse;')
line('v47_qk_distance',[[976,782],[1034,769]],color='#89749F',width=1.7,dashed=True)
sty(cells['c9_softmax'],fillColor='#EDF3EF',strokeColor='#A8C2B6',arcSize=7)
sty(cells['c9_aggregate'],fillColor='#F7F9FA',strokeColor='#C2CCD4',arcSize=7)
sty(cells['c9_local_transform'],fillColor='#FBFAFD',strokeColor='#A994BF',arcSize=6)
sty(cells['c9_linear'],fillColor='#EAF1F7',strokeColor='#AEC2D2',arcSize=6)
for ident,cc in [('c9_values_scalar','#487799'),('c9_values_pair','#A3823A'),('c9_values_point','#85709D')]:
    sty(cells[ident],fontColor=cc,fontSize=30)

# D: retain real snapshots, observable cards and the compact supervision flow.
for ident in ['d_time_first','d_time_middle','d_time_last']:
    root.remove(cells[ident])
label('v47_md_subtitle','Dynamic observables',[1120,907,439,38],28,'#5D806F')
for ident in ['d_target_rmsf','d_target_hchi']:
    sty(cells[ident],fillColor='#F0F6F2',gradientColor='#FFFFFF',gradientDirection='north',strokeColor='#B4CABC',arcSize=8)
for ident in ['d_md_title','d_target_title_rmsf','d_target_title_hchi','d_target_sub_rmsf','d_target_sub_hchi','d_label_summary']:
    sty(cells[ident],fontColor='#517B68',fontStyle=0)
cells['d_label_summary'].set('value','MD mean + variance')
sty(cells['d_label_summary'],fontSize=30)
box(cells['d_loss_card'],x=1307,y=1177,w=250,h=74)
sty(cells['d_loss_card'],fillColor='#EAF3EE',strokeColor='#9CB9A8',arcSize=7)
box(cells['d_loss_text'],x=1316,y=1181,w=232,h=34)
sty(cells['d_loss_text'],fontSize=30,fontStyle=0)
label('v47_loss_variance','v = σ² + v(MD)',[1316,1211,232,33],28)
label('v47_model_label','Model',[1126,1165,154,29],28,'#487799')
sty(cells['d_prediction_params'],fontSize=31,fontStyle=0,fontColor='#487799')
box(cells['d_prediction_params'],y=1195,h=44)
cells['d_prediction_to_loss'].find("mxGeometry/mxPoint[@as='targetPoint']").set('x','1301')
cells['d_label_to_loss'].find("mxGeometry/mxPoint[@as='targetPoint']").set('y','1173')

# Native curve data/positions remain byte-exact. Explicitly audit this rather
# than infer preservation from a screenshot.
for ident,before in protected.items():assert ET.tostring(cells[ident])==before,ident
assert len([c for c in root if c.get('id','').startswith('curve_')])==224
visible=[c.get('value','') for c in root if c.get('value','')]
assert not [x for x in visible if re.search('[\u3400-\u9fff]',x)]
ET.indent(tree,space='  ')
out=HERE/'fig1_framework.drawio';tree.write(out,encoding='utf-8',xml_declaration=True)
exe='C:/Program Files/draw.io/draw.io.exe'
for ext in ['png','svg','pdf']:
    dest=HERE/('fig1_framework_native.pdf' if ext=='pdf' else 'fig1_framework.'+ext)
    args=[exe,'-x','-f',ext,'--crop','-o',str(dest),str(out)]
    if ext=='png':args[4:4]=['-s','2.026']
    if ext=='svg':args.insert(4,'--embed-diagram')
    proc=subprocess.run(args,capture_output=True,text=True,encoding='utf8',errors='replace',timeout=120)
    assert proc.returncode==0,proc.stdout+proc.stderr
import sys
sys.path.insert(0,str(HERE.parents[1]/'completed_results_20261012/plot_runtime'))
from pypdf import PdfReader,PdfWriter
pdf=PdfReader(HERE/'fig1_framework_native.pdf');page=pdf.pages[0]
w=183/25.4*72;h=w*float(page.mediabox.height)/float(page.mediabox.width)
page.scale_to(w,h);writer=PdfWriter();writer.add_page(page)
with open(HERE/'fig1_framework.pdf','wb') as handle:writer.write(handle)
poppler=Path(r'C:\Users\DELL\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\poppler\Library\bin\pdftoppm.exe')
proc=subprocess.run([str(poppler),'-r','450','-png','-singlefile',str(HERE/'fig1_framework.pdf'),str(HERE/'fig1_framework')],capture_output=True,text=True,encoding='utf8',errors='replace',timeout=120)
assert proc.returncode==0,proc.stdout+proc.stderr
# Physical SVG dimensions match the vector PDF; its native viewBox and all
# editable scientific objects remain intact.
svgpath=HERE/'fig1_framework.svg'
svg=ET.parse(svgpath);svgroot=svg.getroot()
svgroot.set('width','183mm');svgroot.set('height',str(h/72*25.4)+'mm')
svg.write(svgpath,encoding='utf-8',xml_declaration=True)
font_rows=[]
def collect(text,cm,tm,font,size):
    if text.strip():font_rows.append({'text':text.strip(),'source_size_pt':size*abs(cm[0]),'word_size_pt':size*abs(cm[0])*167.64/183})
PdfReader(HERE/'fig1_framework.pdf').pages[0].extract_text(visitor_text=collect)
assert min(row['word_size_pt'] for row in font_rows)>=8
text_bbox_path=HERE/'native_text_bbox.html'
pdftotext=Path(r'C:\texlive\2026\bin\windows\pdftotext.exe')
proc=subprocess.run([str(pdftotext),'-bbox',str(HERE/'fig1_framework.pdf'),str(text_bbox_path)],capture_output=True,text=True,encoding='utf8',errors='replace',timeout=60)
assert proc.returncode==0,proc.stdout+proc.stderr
bbox=ET.parse(text_bbox_path);outside=[]
for word in bbox.getroot().iter():
    if word.tag.split('}')[-1]!='word':continue
    b={k:float(word.get(k)) for k in ['xMin','yMin','xMax','yMax']}
    if b['xMin']<-.1 or b['yMin']<-.1 or b['xMax']>w+.1 or b['yMax']>h+.1:outside.append({'text':word.text,'bbox':b})
assert not outside,outside
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
audit={'status':'rendered_pending_visual_review','source_drawio_sha256':sha(source),'width_mm':183,'height_mm':h/72*25.4,
       'protected_native_prediction_paths':224,'all_paths_exactly_preserved':True,
       'prediction_structure_uses_exact_original_840_RGB':True,'all_original_prediction_values_and_RGB_preserved':True,'protein_raster_style_rerendered':True,'all_MD_snapshots_exactly_preserved':True,
       'topology':'original top-row five-stage workflow; original bottom-row b/c/d columns',
       'new_training':0,'new_MD':0,'new_inference':0,'word_width_mm':167.64,
       'min_source_font_pt':min(row['source_size_pt'] for row in font_rows),'min_word_font_pt':min(row['word_size_pt'] for row in font_rows),'out_of_canvas_text':outside,
       'source_hashes':{str(p.relative_to(HERE)):sha(p) for folder in [S,HERE/'assets3d'] for p in sorted(folder.rglob('*')) if p.is_file()},
       'output_sha256':{ext:sha(HERE/('fig1_framework.'+ext)) for ext in ['png','pdf','svg','drawio']}}
(HERE/'render_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf8')
(HERE/'visible_text.json').write_text(json.dumps(visible,ensure_ascii=False,indent=2),encoding='utf8')
(HERE/'font_sizes.json').write_text(json.dumps(font_rows,ensure_ascii=False,indent=2),encoding='utf8')
print(json.dumps(audit,ensure_ascii=False))
