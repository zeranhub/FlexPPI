from pathlib import Path
import json
import numpy as np,pandas as pd
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028')
p=pd.read_csv(R/'paired_all_residue_predictions.csv',dtype={'author':str})
n=pd.read_csv(R/'paired_numbering.csv',dtype={'author':str})
p=p.merge(n,on=['state_id','chain','author','aa'],how='left',validate='many_to_one')
p['author_number']=p.author.str.extract(r'(-?\d+)').astype(int)
cons=json.loads((R/'pair_screen_contrasts.json').read_text());all_d=[];summ=[]
for con in cons:
 for v in ['full','intrachain']:
  b=p[p.state_id.eq(con['before'])&p.variant.eq(v)];a=p[p.state_id.eq(con['after'])&p.variant.eq(v)]
  d=b.merge(a,on=['chain','author'],suffixes=('_before','_after'),validate='one_to_one')
  d['same_aa']=d.aa_before.eq(d.aa_after);d['contrast']=con['id'];d['variant']=v;d['counterfactual']=con['counterfactual'];d['comparison_type']=con['type']
  d['delta_rmsf']=d.rmsf_after-d.rmsf_before;d['delta_hchi']=d.hchi_after-d.hchi_before
  d['rmsf_eligible']=d.model_valid_before&d.model_valid_after&d.geometry_valid_before&d.geometry_valid_after
  d['hchi_eligible']=d.rmsf_eligible&d.hchi_valid_before&d.hchi_valid_after&d.same_aa
  masks={'all_shared':np.ones(len(d),bool),'all_unchanged':d.same_aa,'complete_variable_domains':d.imgt_number_before.notna()&d.imgt_number_after.notna()}
  for ch in ['H','L','A']:
   masks['chain_'+ch]=d.chain.eq(ch)
   if ch!='A':
    for reg in ['CDR1','CDR2','CDR3']:masks[ch+'_'+reg]=d.chain.eq(ch)&d.region_before.eq(reg)&d.region_after.eq(reg)
  net=np.zeros(len(d),bool);around=np.zeros(len(d),bool)
  for ch,num in [('H',31),('H',97),('A',67)]:
   net|=(d.chain.eq(ch)&d.author_number_before.eq(num)).to_numpy()
   around|=(d.chain.eq(ch)&d.author_number_before.between(num-2,num+2)).to_numpy()
  masks['predefined_network_three_sites']=net;masks['predefined_network_complete_neighbors_plusminus2']=around
  all_d.append(d)
  for scope,mask in masks.items():
   for metric in ['rmsf','hchi']:
    z=d[mask&d[metric+'_eligible']]
    if z.empty:continue
    summ.append({'contrast':con['id'],'variant':v,'comparison_type':con['type'],'counterfactual':con['counterfactual'],'scope':scope,'metric':metric,'n':len(z),
      'mean_before':z[metric+'_before'].mean(),'mean_after':z[metric+'_after'].mean(),'mean_delta':z['delta_'+metric].mean(),'mean_absolute_delta':z['delta_'+metric].abs().mean(),
      'n_increase':int((z['delta_'+metric]>0).sum()),'n_decrease':int((z['delta_'+metric]<0).sum())})
d=pd.concat(all_d,ignore_index=True);s=pd.DataFrame(summ)
d.to_csv(R/'paired_complete_position_deltas.csv',index=False);s.to_csv(R/'paired_predefined_region_deltas.csv',index=False)
network=d[(d.chain.eq('H')&d.author_number_before.isin([31,97]))|(d.chain.eq('A')&d.author_number_before.eq(67))]
assert (network.aa_before==network.aa_after).all()
network.to_csv(R/'paired_predefined_network_all_outputs.csv',index=False)
print(s[s.variant.eq('full')&s.scope.isin(['H_CDR3','predefined_network_three_sites','complete_variable_domains'])&s.contrast.isin(['P103A_minus_WT_copy1','P103A_minus_WT_copy2','WT_copy2_minus_copy1','P103A_copy2_minus_copy1'])].to_string(index=False))
print('\nNetwork sites, all observed comparisons (not validated entropy/affinity truth)')
print(network[network.variant.eq('full')&network.contrast.isin(['P103A_minus_WT_copy1','P103A_minus_WT_copy2'])][['contrast','chain','author','aa_before','rmsf_before','rmsf_after','delta_rmsf','hchi_before','hchi_after','delta_hchi']].to_string(index=False))
