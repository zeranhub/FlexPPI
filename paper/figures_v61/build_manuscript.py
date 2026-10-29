"""Insert revised Figure 1 and make local caption edits preserving native math."""
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
from lxml import etree as ET
from PIL import Image

ROOT=Path(__file__).resolve().parent
BASE=Path('E:/Project/FlexPPI/paper/revision_v14_20261010/paper_v60_20261023_中文修订_图1配色轻调.docx')
OUTPUT=ROOT.parent/'paper_v61_20261023_中文修订_图1输入输出标注修正.docx'
FIGURE=ROOT/'01_模型框架/fig1_framework.png'
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
    'm':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
sha=lambda data:hashlib.sha256(data).hexdigest()

def main():
    assert Image.open(FIGURE).size==(3802,2288)
    with ZipFile(BASE) as source:
        original={name:source.read(name) for name in source.namelist()}
        doc=ET.fromstring(original['word/document.xml'])
        captions=[p for p in doc.findall('.//w:p',NS) if ''.join(p.xpath('.//w:t/text()',namespaces=NS)).startswith('图1：')]
        assert len(captions)==1
        caption=captions[0]
        equations=[ET.tostring(m) for m in caption.findall('.//m:oMath',NS)]
        assert len(equations)==4
        replacements={'氨基酸类型及':'残基序列及',
            '（c）共享读出层':'（c）Hχ为全部可定义χ角的联合状态熵；共享读出层',
            '推理输入为静态主链几何与类型特征':'推理输入为静态主链几何与残基序列'}
        counts={s:0 for s in replacements}
        for text in caption.findall('.//w:t',NS):
            for before,after in replacements.items():
                if before in (text.text or ''):
                    text.text=text.text.replace(before,after)
                    counts[before]+=1
        assert all(n==1 for n in counts.values()),counts
        assert [ET.tostring(m) for m in caption.findall('.//m:oMath',NS)]==equations
        changed={'word/media/image1.png':FIGURE.read_bytes(),
            'word/document.xml':ET.tostring(doc,encoding='UTF-8',xml_declaration=True,standalone=True)}
        with ZipFile(OUTPUT,'w') as destination:
            for info in source.infolist():
                destination.writestr(info,changed.get(info.filename,original[info.filename]))
    with ZipFile(OUTPUT) as result:
        differences=[name for name,data in original.items() if result.read(name)!=data]
        assert set(differences)==set(changed)
    (ROOT/'document_QA').mkdir(exist_ok=True)
    audit={'status':'local_edit_verified_render_pending','baseline':str(BASE),'output':str(OUTPUT),
        'modified_members':differences,'caption_replacements':replacements,
        'native_caption_math_objects_preserved':4,'all_other_members_byte_identical':True,
        'document_sha256':sha(OUTPUT.read_bytes()),'figure1_sha256':sha(FIGURE.read_bytes())}
    (ROOT/'document_QA/document_edit_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(OUTPUT)

if __name__=='__main__':main()
