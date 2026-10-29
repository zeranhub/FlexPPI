from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from PIL import Image
import json,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parents[2]
SRC=BASE/'paper_v97_20261027_中文修订_期刊参考组图重排.docx'
DEST=BASE/'paper_v98_20261027_中文修订_结构区域放大.docx'
sha=lambda b:hashlib.sha256(b).hexdigest()
assert sha(SRC.read_bytes())=='3e22f3ee12b4825e13d7400f8674750de0667b9788432863b3d5069d20d334cc'
with ZipFile(SRC) as z:parts={n:z.read(n) for n in z.namelist()}
old=parts.copy()
png=ROOT/'figure3/fig3_dual_case_entropy.png'
assert Image.open(png).size==Image.open(ROOT.parent/'revision_v97_journal_composition/figure3/fig3_dual_case_entropy.png').size
parts['word/media/image5.png']=png.read_bytes()
assert [n for n in parts if old[n]!=parts[n]]==['word/media/image5.png']
with ZipFile(DEST,'w',ZIP_DEFLATED) as z:
    for n,data in parts.items():z.writestr(n,data)
with ZipFile(DEST) as z:assert z.testzip() is None
for ext in ['png','pdf','svg']:
    (BASE/f'figure3_v98_20261027_结构区域放大.{ext}').write_bytes((ROOT/f'figure3/fig3_dual_case_entropy.{ext}').read_bytes())
audit=dict(source=str(SRC),source_sha256=sha(SRC.read_bytes()),output=str(DEST),output_sha256=sha(DEST.read_bytes()),changed_zip_parts=['word/media/image5.png'],all_document_text_caption_references_styles_tables_and_other_media_unchanged=True,figure_physical_size_unchanged=True)
(ROOT/'manuscript/build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False))
