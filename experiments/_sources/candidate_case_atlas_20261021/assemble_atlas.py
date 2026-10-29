"""Bind the complete case atlas to its actual final sources and figure assets."""
from pathlib import Path
import json,hashlib,re
ROOT=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2),encoding='utf8')

def collect_figures():
    f={}
    d=ROOT/'figures/C01'
    cap='MEL5亲和成熟。a，α24β17复合物4JFF的真实结构及Q72与βPhe53/βGly54主链接触，亲和力来自原研究SPR。b，六个完整CDR的45个位置的模型绝对RMSF，共同编号中工程构建体输入整体加1。c，全部23个有效且氨基酸类型相同的pMHC界面位置的模型Hχ。d，已有ITrace三重复轨迹中Q72–Gly54近接触及Q72两端同时符合氢键几何的帧占比，点为独立重复，短竖线为等权均值。近接触阈值4 Å；氢键几何为供体–受体距离≤3.5 Å及供体–氢–受体角≥150°。两个模型输出与轨迹接触分别承担动态定位和机制佐证。'
    (d/'caption_zh.txt').write_text(cap,encoding='utf8')
    f['C01']={'case_id':'C01',**{ext:str(d/f'fig4_mel5_joint.{ext}') for ext in ['png','svg','pdf']},
              'caption':str(d/'caption_zh.txt'),'audit':str(d/'data_audit.json')}
    p=ROOT/'figure_design/three_case_figures_manifest.json'
    for x in read(p):
        f[x['case_id']]={'case_id':x['case_id'],**{k:x[k.upper()] for k in ['png','pdf','svg']},
          'caption':str(Path(x['directory'])/'caption_zh.md'),'audit':str(Path(x['directory'])/'data_audit.json')}
    for rel in ['figure_design/biological_series_manifest.json','data_and_general_figures/general_figure_manifest.json','evidence/non_tcr_figures/non_tcr_figure_manifest.json']:
        p=ROOT/rel
        if not p.exists():continue
        data=read(p);items=data if isinstance(data,list) else data['figures']
        for x in items:
            cid=x.get('case_id',x.get('caseid'))
            record={'case_id':cid}
            for key in ['png','pdf','svg','caption','audit']:
                value=x.get(key,x.get(key.upper()))
                if value is not None:record[key]=str(value)
            if 'caption' not in record:
                directory=Path(record['png']).parent
                choices=[directory/'caption_zh.txt',directory/'caption_zh.md']
                record['caption']=str(next(p for p in choices if p.exists()))
            f[cid]=record
    expected={f'C{i:02}' for i in range(1,33)}
    assert set(f)==expected,{'missing':sorted(expected-set(f)),'extra':sorted(set(f)-expected)}
    for cid,record in f.items():
        for key in ['png','pdf','svg','caption']:
            assert Path(record[key]).is_file(),(cid,key,record.get(key))
        record['sha256']={key:sha(record[key]) for key in ['png','pdf','svg','caption']}
    return {k:f[k] for k in sorted(f)}

def main():
    figures=collect_figures()
    narratives=read(ROOT/'evidence/narratives.json')
    byid={x['catalog_case_id']:x for x in narratives['cases']}
    assert set(byid)==set(figures)
    references=[];refkeys={};case_sources={}
    for cid,c in sorted(byid.items()):
        ids=[]
        for s in c['sources']:
            key=(s.get('doi') or s.get('url') or s['title']).lower()
            if key not in refkeys:
                refkeys[key]=len(references)+1
                source=dict(s);source['number']=refkeys[key]
                if source.get('local') and Path(source['local']).exists():source['current_local_sha256']=sha(source['local'])
                references.append(source)
            ids.append(refkeys[key])
        case_sources[cid]=ids
    manifest={'status':'assets_bound_pending_document_QA','cases':32,'figures':figures,
              'references':references,'case_references':case_sources,
              'narratives':str(ROOT/'evidence/narratives.json'),'narratives_sha256':sha(ROOT/'evidence/narratives.json'),
              'main_manuscript':str(ROOT.parent.parent/'paper_v49_20261019_中文修订_MEL5双动态机制实例_三线表.docx'),
              'main_manuscript_sha256':'5f79ddf7d657dc9b782c9fb432e373a1e7c170c1c88d019e528fbf73c366c222'}
    assert sha(manifest['main_manuscript'])==manifest['main_manuscript_sha256']
    write(ROOT/'ATLAS_FIGURE_MAP.json',manifest)
    print(json.dumps({'status':manifest['status'],'cases':32,'figures':len(figures),'references':len(references)}))

if __name__=='__main__':main()
