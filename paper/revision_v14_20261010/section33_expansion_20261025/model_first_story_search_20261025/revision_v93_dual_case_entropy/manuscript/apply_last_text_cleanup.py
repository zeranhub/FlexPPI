from pathlib import Path
import json,subprocess,sys
HERE=Path(__file__).resolve().parent

script=HERE/'prepare_global_sync.py'
text=script.read_text(encoding='utf-8')
old='仅按预测确定的高值组实际角区平均宽7.76°/角'
new='仅按预测划分的高值组，其实际角区较低值组平均宽7.76°/角'
assert text.count(old)==1
script.write_text(text.replace(old,new),encoding='utf-8')
subprocess.run([sys.executable,str(script)],check=True)

path=HERE/'final_text_contract.json'
cfg=json.loads(path.read_text(encoding='utf-8'))
old='泛素逐角分布仅作为这种对应的局部诊断，移至图S4c、d，不将MD自由能曲线误写为模型输出。'
new='泛素逐角分布的局部核查见图S4c、d；其中自由能曲线由MD占据统计得到，不是模型输出。'
assert cfg['section34']['paragraphs'][-1].count(old)==1
cfg['section34']['paragraphs'][-1]=cfg['section34']['paragraphs'][-1].replace(old,new)
cap=Path(cfg['supplement4_caption_file'])
text=cap.read_text(encoding='utf-8')
old='面板b来自原pipeline留出预测，非本轮新家族划分的类别分析。'
new='面板b采用原体系划分的折外预测。'
assert text.count(old)==1
text=text.replace(old,new)
cap.write_text(text,encoding='utf-8')
cfg['supplement4_caption']=text.strip()
cfg['paragraph_replacements']=json.loads((HERE/'global_paragraph_replacements_proposed.json').read_text(encoding='utf-8'))
path.write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf-8')
print('Three final text corrections synchronized in builder source, contracts, global proposal and S4 caption source.')
