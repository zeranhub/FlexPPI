from pathlib import Path
from lxml import etree
from collections import Counter
from docx.oxml.ns import qn
import zipfile
r = Path(__file__).resolve().parents[2]
with zipfile.ZipFile(r / 'paper_v18_20261010_中文修订_复用KIX单轨迹_三线表.docx') as a, zipfile.ZipFile(r / 'paper_v19_20261012_中文修订_补充实验完成_三线表.docx') as b:
    x, y = a.read('word/comments.xml'), b.read('word/comments.xml')
    texts = lambda raw: [e.text for e in etree.fromstring(raw).iter(qn('w:t'))]
    anchors = lambda z: Counter((e.tag, e.get(qn('w:id'))) for e in etree.fromstring(z.read('word/document.xml')).iter() if e.tag in {qn('w:commentRangeStart'), qn('w:commentRangeEnd'), qn('w:commentReference')})
    print('lengths', len(x), len(y))
    print('comment_text_identical', texts(x) == texts(y))
    print('anchors_identical', anchors(a) == anchors(b))
    print('media_identical', all(a.read(n) == b.read(n) for n in a.namelist() if n.startswith('word/media/') and not n.endswith('/')))
    print('canonical_identical', etree.tostring(etree.fromstring(x), method='c14n') == etree.tostring(etree.fromstring(y), method='c14n'))
    print('original_prefix', x[:160])
    print('saved_prefix', y[:160])
