"""Bind the reviewed Word atlas and archive its 32 final figure assets."""
from pathlib import Path
import hashlib,json,zipfile

ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent.parent
QA=ROOT/'document_QA'
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf8')

def main():
    m=read(ROOT/'ATLAS_FIGURE_MAP.json')
    integration=read(QA/'integration_audit.json')
    doc=Path(integration['document'])
    ds=sha(doc)
    b=read(QA/'render_binding.json')
    assert ds==b['document_sha256']==integration['document_sha256']
    assert b['pages']==74
    assert sha(m['main_manuscript'])==m['main_manuscript_sha256']
    for name in ['independent_document_audit.json','final_narrative_review.json']:
        d=read(QA/name)
        assert d['status']=='passed' and d['document_sha256']==ds,name
    scope=read(ROOT/'catalog/final_scope_and_numbers_review.json')
    assert scope['status']=='passed' and scope['final_figure_map_SHA256']==sha(ROOT/'ATLAS_FIGURE_MAP.json')
    visual=[]
    pages={}
    for name in ['page1_13_68_74_visual_QA.json','page14_31_visual_QA.json','page32_49_visual_QA.json','page50_67_visual_QA.json']:
        d=read(QA/name)
        assert d['status']=='passed' and d['document_sha256']==ds,name
        for p in d['pages']:
            n=str(p['page'])
            assert n not in pages
            assert p['status']=='passed' and p['PNG_SHA256']==b['PNG_SHA256'][n]
            assert sha(Path(b['render_directory'])/f'page-{n}.png')==p['PNG_SHA256']
            pages[n]=p
        visual.append({'report':str(QA/name),'sha256':sha(QA/name)})
    assert set(pages)==set(map(str,range(1,b['pages']+1)))
    figures=[]
    for cid,d in sorted(m['figures'].items()):
        for ext in ['png','pdf','svg','caption']:
            assert sha(d[ext])==d['sha256'][ext],(cid,ext)
        figures.append({'case_id':cid,'png_sha256':d['sha256']['png'],'caption_sha256':d['sha256']['caption']})
    archive=BASE/'FlexPPI_32个案例_英文大图_PNG_PDF_SVG_20261021.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for cid,d in sorted(m['figures'].items()):
            for ext in ['png','pdf','svg']:
                z.write(d[ext],f'{cid}/{cid}_case_figure.{ext}')
            z.write(d['caption'],f'{cid}/caption_zh.txt')
        for source,name in [(ROOT/'catalog/all_discovery_pair_index.json','catalog/all_discovery_pair_index.json'),(ROOT/'evidence/case_assessment_table.json','catalog/case_assessment_table.json')]:
            z.write(source,name)
        z.writestr('README_zh.txt','本包含32例最终英文大图，每例提供450 dpi PNG、PDF与SVG，以及中文图注。C01–C30为TCR候选系列，C31为KIX，C32为泛素。定量位点、对照及论证范围详见独立Word图册；并非32个成功验证。SVG保留矢量文字和曲线，结构视图作为图内栅格图像。原始候选配对索引包括343个结构对，不等于343个独立机制案例。\n')
        z.writestr('figure_content_hashes.json',json.dumps(figures,ensure_ascii=False,indent=2))
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None and len([n for n in z.namelist() if n.endswith('.png')])==32
    b['status']='passed';b['all_pages_visually_reviewed']=True
    write(QA/'render_binding.json',b)
    integration.update(status='passed',pages=b['pages'],all_pages_visually_reviewed=True,visual_QA_reports=visual)
    write(QA/'integration_audit.json',integration)
    report={'status':'complete','document':str(doc),'document_sha256':ds,'pages':b['pages'],'case_count':32,'figure_count':32,'TCR_series':30,'non_TCR_cases':2,'native_three_line_tables':34,'references':37,'raw_discovery_pairs':343,'figure_archive':str(archive),'figure_archive_sha256':sha(archive),'figure_map':str(ROOT/'ATLAS_FIGURE_MAP.json'),'figure_map_sha256':sha(ROOT/'ATLAS_FIGURE_MAP.json'),'all_pages_visually_reviewed':True,'visual_QA_reports':visual,'main_manuscript':m['main_manuscript'],'main_manuscript_sha256':m['main_manuscript_sha256'],'main_manuscript_unchanged':True,'new_training':0,'new_inference':0,'new_MD':0,'main_manuscript_delivery_pointer_updated':False}
    write(ROOT/'delivery_status.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
