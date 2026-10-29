from pathlib import Path
import pandas as pd,json
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
for tag in ['downsample','full']:
 d=pd.read_csv(r/'training'/('train_'+tag+'.csv'),sep=';')
 print(tag, d.shape, d.columns.tolist())
 for cdr in ['CAWSETGLGTGELFF','CAWSETGLGMGGWQF']:
  q=d[d.cdr3==cdr];print('CDR3',cdr,'rows',q.to_dict('records'))
 for pep in ['ELAGIGILTV','ELAAIGILTV','ELAGIGALTV','ILAKFLHWL','SLLMWITQC','YLEPGPVTA']:
  q=d[d['antigen.epitope']==pep];print('epitope',pep,len(q),q.y.value_counts().to_dict())
