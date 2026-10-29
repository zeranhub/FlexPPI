from pathlib import Path
from pypdf import PdfReader
from pypdf.generic import ContentStream

p=PdfReader(r'E:\Project\FlexPPI\paper\revision_v14_20261010\figure2_relayout_20261024\figS3_internal_PR.pdf').pages[0]
c=ContentStream(p.get_contents(),p.pdf)
for i,(a,o) in enumerate(c.operations):
    if o in [b'Tj',b'TJ']:
        print(i,repr(a),o)
        if str(a) in ["[['a']]", "[['b']]"]:
            print('CONTEXT', c.operations[max(i-6,0):i+3])
