import sys,json,csv,hashlib
from pathlib import Path
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028')
sys.path.insert(0,str(BASE/'numbering/deps'))
from Bio import SeqIO
from anarcii import Anarcii
ROOT=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_absolute_affinity_20261028\affinity\structures')
seqs={}
for pdb in ['3hg1','4jff','4jfd','4jfe']:
 for rec in SeqIO.parse(str(ROOT/f'{pdb}_crystal_raw.pdb'),'pdb-seqres'):
  chain=rec.annotations['chain']
  if chain in ['D','E']:seqs[pdb+'_'+chain]=str(rec.seq)
(BASE/'numbering/crystal_SEQRES_sequences.json').write_text(json.dumps(seqs,indent=2),encoding='utf-8')
numberer=Anarcii(seq_type='tcr',mode='accuracy',cpu=True,ncpu=2,verbose=True)
results=numberer.number(seqs)
(BASE/'numbering/anarcii_raw_numbering.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
