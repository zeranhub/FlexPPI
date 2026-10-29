from pathlib import Path
import json,zipfile
root=Path(__file__).resolve().parent;base=root.parents[2]
failure=json.loads((root/'generation_01/failure.json').read_text(encoding='utf-8'))
assert failure['status']==402 and failure['stage']=='Retriever'
assert not list((root/'generation_01').glob('*final.png'))
status={'actual_PaperBanana_called':True,'dry_run_passed':True,'completed':False,'failure':failure,
        'new_layout_generated':False,'Word_or_Figure3_changed':False,'input_ready_for_resume':True,
        'credentials_packaged':False}
(root/'run_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2),encoding='utf-8')
package=base/'figure3_paperbanana_v111_20261027_新版面输入与调用诊断.zip'
with zipfile.ZipFile(package,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(root.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(root).as_posix())
with zipfile.ZipFile(package) as z:assert z.testzip() is None
print(json.dumps({'package':str(package),'bytes':package.stat().st_size,'actual_api_status':402,'new_PaperBanana_layout_generated':False},ensure_ascii=False))
