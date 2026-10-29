from pathlib import Path
from zipfile import ZipFile
from lxml import etree
from PIL import Image
import io,json
base=Path(r'E:\Project\FlexPPI\paper\revision_v14_20261010');v=base/'section33_expansion_20261025/model_first_story_search_20261025/revision_v83_model_story';out=v/'palette_audit';out.mkdir(exist_ok=True)
ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
with ZipFile(base/'paper_v82_20261025_中文修订_图3统一视角与局部并排.docx') as z:
    ps=etree.fromstring(z.read('word/document.xml')).xpath('//w:body/w:p',namespaces=ns);rels=etree.fromstring(z.read('word/_rels/document.xml.rels'));relmap={n.get('Id'):n.get('Target') for n in rels}
    found=[]
    for i,p in enumerate(ps):
        ids=p.xpath('.//a:blip/@r:embed',namespaces=ns)
        if not ids:continue
        nxt=''.join(ps[i+1].xpath('.//w:t/text()',namespaces=ns)) if i+1<len(ps) else ''
        if not nxt.startswith('图'):continue
        for rid in ids:
            name=relmap[rid];b=z.read('word/'+name);im=Image.open(io.BytesIO(b));path=out/f'paragraph_{i}_{Path(name).name}';path.write_bytes(b);found.append({'index':i,'file':str(path),'size':im.size,'caption':nxt[:250]})
(out/'figure_inventory.json').write_text(json.dumps(found,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(found,ensure_ascii=False,indent=2))