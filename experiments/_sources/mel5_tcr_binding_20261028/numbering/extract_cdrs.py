import json
from pathlib import Path
p=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\numbering')
r=json.loads((p/'anarcii_raw_numbering.json').read_text())
s={}
for name,x in r.items():
 assert x['error'] is None
 def part(lo,hi):return ''.join(aa for (n,ins),aa in x['numbering'] if lo<=n<=hi and aa!='-')
 s[name]={'chain_type':x['chain_type'],'CDR1_IMGT27_38':part(27,38),'CDR2_IMGT56_65':part(56,65),'CDR3_IMGT105_117':part(105,117),'CDR3_with_C104_F118':part(104,118),'score':x['score']}
(p/'verified_CDR_sequences.json').write_text(json.dumps(s,indent=2),encoding='utf-8')
print(json.dumps(s,indent=2))
