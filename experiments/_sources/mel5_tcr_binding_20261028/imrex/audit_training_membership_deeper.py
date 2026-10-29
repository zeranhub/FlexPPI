from pathlib import Path
import pandas as pd,json,re,hashlib
r=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028\imrex')
q='CAWSETGLGMGGWQF'
variants={'full':q,'without_initial_C':q[1:],'without_terminal_F':q[:-1],'without_both_anchors':q[1:-1]}
def basic(s):return re.sub('[^A-Z]','',str(s).upper())
def core(s):
 s=basic(s)
 if s.startswith('C'):s=s[1:]
 if s.endswith('F') or s.endswith('W'):s=s[:-1]
 return s
def lev(a,b):
 p=list(range(len(b)+1))
 for i,x in enumerate(a,1):
  c=[i]
  for j,y in enumerate(b,1):c.append(min(c[-1]+1,p[j]+1,p[j-1]+(x!=y)))
  p=c
 return p[-1]
audit={'query_mature_CDR3_beta':q,'query_length':len(q),'variants':variants,'normalization':'Basic strips non-uppercase amino-acid symbols/whitespace; core removes at most one initial C and one terminal F or W, applied equally to query and training. Does not remove internal residues.','scope':'Only released actual checkpoint train CSVs, not unknown external or pretraining databases. No V/J or full receptor sequence available in these files; exact absence does not establish family exclusion.','checkpoints':{}}
rows=[]
for tag in ['downsample','full']:
 p=r/'training'/('train_'+tag+'.csv');d=pd.read_csv(p,sep=';');seq=sorted(set(d.cdr3.astype(str)));d['norm_cdr3']=d.cdr3.map(basic);d['core_cdr3']=d.cdr3.map(core)
 nearest=sorted([(lev(q,basic(s)),s) for s in seq],key=lambda x:(x[0],x[1]))[:5]
 near=[]
 for rank,(distance,s) in enumerate(nearest,1):
  subset=d[d.cdr3==s]
  pairs=subset[['y','antigen.epitope']].drop_duplicates().to_dict('records')
  one={'rank':rank,'cdr3':s,'edit_distance':distance,'length':len(s),'training_rows':len(subset),'pairs':pairs}
  near.append(one)
  for pair in pairs:rows.append({'checkpoint':tag,'rank':rank,'query':q,'training_cdr3':s,'edit_distance':distance,'training_label':int(pair['y']),'peptide':pair['antigen.epitope'],'training_rows_for_cdr3':len(subset)})
 ep={}
 for pep in ['ELAAIGILTV','ELAGIGALTV','ELAGIGILTV']:
  sub=d[d['antigen.epitope'].map(basic)==pep]
  ep[pep]={'rows':len(sub),'unique_cdr3':sub.cdr3.nunique(),'label_counts':{str(k):int(v) for k,v in sub.y.value_counts().items()}}
 exact={k:{'raw_rows':int((d.cdr3==v).sum()),'basic_normalized_rows':int((d.norm_cdr3==v).sum())} for k,v in variants.items()}
 audit['checkpoints'][tag]={'train_csv':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'training_rows':len(d),'unique_raw_CDR3':len(seq),'query_variant_matches':exact,'both_sides_core_normalization_rows':int((d.core_cdr3==core(q)).sum()),'G4A_I7A_and_ELA_peptide_presence':ep,'nearest_5_distinct_CDR3_full_sequence':near,'min_core_distance':min(lev(core(q),core(s)) for s in seq)}
 # Also retain all exact normalized mature matches, should any appear.
 d[d.core_cdr3==core(q)].to_csv(r/('normalized_mature_training_matches_'+tag+'.csv'),index=False)
pd.DataFrame(rows).to_csv(r/'mature_nearest_training_CDR3.csv',index=False)
(r/'training_membership_deeper_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
