from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
B=Path(r'E:\Project\FlexPPI');P=B/'paper/revision_v14_20261010/section33_expansion_20261025';S=P/'model_first_story_search_20261025';O=S/'class_level_feasibility_20261026/inventory_agent';D=B/'experiments/non_tcr_case_screen_20261025';F=B/'experiments/_sources/family_folds'
fam=pd.read_csv(F/'record_family_folds.csv'); co=pd.read_csv(P/'matched_complementarity/cohort83_record_family_inventory.csv');records=pd.read_csv(D/'all_616_record_screen.csv').merge(fam[['sample_id','component_id','short_chains']],on='sample_id',validate='one_to_one'); chain=pd.read_csv(D/'all_chain_metrics.csv'); seq=pd.read_csv(B/'experiments/_sources/baselines/cohort/chain_manifest.csv'); feats=pd.read_parquet(P/'feature_baseline/features_all.parquet',columns=['sample_id','chain','residue_index','amino_acid','partner_ca_contacts_8A','interface_ca8','valid_rmsf','valid_hchi']);ca=feats.groupby(['sample_id','chain']).agg(n_interface=('interface_ca8','sum'),n_valid_rmsf=('valid_rmsf','sum'),n_valid_hchi=('valid_hchi','sum'),contact_count=('partner_ca_contacts_8A','sum')).reset_index();ch=chain.merge(seq[['sample_id','chain','sequence']],on=['sample_id','chain'],validate='one_to_one').merge(ca,on=['sample_id','chain'],validate='one_to_one').merge(records[['sample_id','component_id','n_chains','n_protein_identities']],on='sample_id',validate='many_to_one')
records.to_csv(O/'all616_records_with_family.csv',index=False);ch.to_csv(O/'all1535_chains_with_sequence_contacts_identity_family.csv',index=False)
short=[]
for sid,g in ch.groupby('sample_id',sort=False):
 for _,a in g[g.n_residues.between(5,40)&g.n_interface.gt(0)&g.n_valid_rmsf.gt(0)&g.n_valid_hchi.gt(0)].iterrows():
  b=g[g.n_residues.ge(70)&g.n_interface.gt(0)&g.n_valid_rmsf.gt(0)&g.n_valid_hchi.gt(0)&g.identity.ne(a.identity)]
  if b.empty:continue
  short.append(dict(sample_id=sid,pdb_id=a.pdb_id,family_id=a.component_id,short_chain=a.chain,short_length=int(a.n_residues),short_description=a.description,short_identity=a.identity,short_sequence=a.sequence,short_interface_n=int(a.n_interface),short_hchi_valid_n=int(a.n_valid_hchi),receptor_chains='|'.join(b.chain),receptor_descriptions=' | '.join(b.description),receptor_identities=' | '.join(b.identity),receptor_lengths='|'.join(b.n_residues.astype(str)),two_chain_unambiguous=int(a.n_chains)==2))
short=pd.DataFrame(short);short.to_csv(O/'short_peptide_5to40_candidates.csv',index=False)
# receptor exactly matching sequence >=70 aa; protein identity of other chains used only as membership labels, never model score.
contexts=[]
for sid,g in ch.groupby('sample_id',sort=False):
 for _,a in g[(g.n_residues>=70)&(g.n_interface>0)&(g.n_valid_rmsf>0)&(g.n_valid_hchi>0)].iterrows():
  pp=g[g.chain.ne(a.chain)&g.identity.ne(a.identity)&g.n_interface.gt(0)]
  if pp.empty:continue
  contexts.append(dict(sample_id=sid,pdb_id=a.pdb_id,family_id=a.component_id,receptor_chain=a.chain,receptor_n=int(a.n_residues),receptor_description=a.description,receptor_identity=a.identity,receptor_sequence=a.sequence,sequence_sha256=hashlib.sha256(a.sequence.encode()).hexdigest(),partner_signature=' || '.join(sorted(set(pp.identity))),partner_descriptions=' || '.join(sorted(set(pp.description))),two_chain_unambiguous=int(a.n_chains)==2))
ctx=pd.DataFrame(contexts);gr=[]
for hs,g in ctx.groupby('sequence_sha256',sort=False):
 if g.pdb_id.nunique()<2 or g.partner_signature.nunique()<2:continue
 gr.append(dict(sequence_sha256=hs,family_id='|'.join(sorted(set(g.family_id))),receptor_identity=' || '.join(sorted(set(g.receptor_identity))),receptor_description=' || '.join(sorted(set(g.receptor_description))),receptor_n=int(g.receptor_n.iloc[0]),n_records=g.sample_id.nunique(),n_pdb=g.pdb_id.nunique(),n_partner_signatures=g.partner_signature.nunique(),sample_ids='|'.join(sorted(set(g.sample_id))),partner_descriptions=' ||| '.join(sorted(set(g.partner_descriptions))),all_two_chain=bool(g.two_chain_unambiguous.all())))
groups=pd.DataFrame(gr);groups.to_csv(O/'exact_receptor_different_partner_candidates.csv',index=False);ctx[ctx.sequence_sha256.isin(groups.sequence_sha256)].to_csv(O/'exact_receptor_context_members.csv',index=False)
# Conservative two-chain counterpart excludes multichain signature ambiguity.
gr2=[]
for hs,g in ctx[ctx.two_chain_unambiguous].groupby('sequence_sha256',sort=False):
 if g.pdb_id.nunique()<2 or g.partner_signature.nunique()<2:continue
 gr2.append(dict(sequence_sha256=hs,family_id='|'.join(sorted(set(g.family_id))),receptor_identity=' || '.join(sorted(set(g.receptor_identity))),receptor_description=' || '.join(sorted(set(g.receptor_description))),receptor_n=int(g.receptor_n.iloc[0]),n_records=g.sample_id.nunique(),n_pdb=g.pdb_id.nunique(),n_partner_signatures=g.partner_signature.nunique(),sample_ids='|'.join(sorted(set(g.sample_id))),partner_descriptions=' ||| '.join(sorted(set(g.partner_descriptions)))))
groups2=pd.DataFrame(gr2);groups2.to_csv(O/'two_chain_exact_receptor_different_partner_candidates.csv',index=False)
# descriptions-only enzyme/protein inhibitor candidate; separate from broad protease keyword.
ep=[]
for sid,g in ch.groupby('sample_id',sort=False):
 enz=g[g.description.str.contains(r'protease|proteinase|trypsin|thrombin|chymotrypsin|elastase|subtilisin|thermolysin|collagenase|stromelysin|ribonuclease|amylase',case=False,na=False)&~g.description.str.contains('inhibitor|ATP-binding|adaptor',case=False,na=False)]
 inh=g[g.description.str.contains(r'inhibitor|barstar|tendamistat|elafin|eglin',case=False,na=False)]
 for _,a in enz.iterrows():
  b=inh[inh.chain.ne(a.chain)]
  if b.empty:continue
  ep.append(dict(sample_id=sid,pdb_id=a.pdb_id,family_id=a.component_id,enzyme_chain=a.chain,enzyme_description=a.description,inhibitor_chains='|'.join(b.chain),inhibitor_descriptions=' | '.join(b.description),min_inhibitor_length=int(b.n_residues.min()),n_chains=int(a.n_chains),classification_status='description_keyword_candidate_manual_confirmation_needed'))
ep=pd.DataFrame(ep);ep.to_csv(O/'enzyme_protein_inhibitor_description_candidates.csv',index=False)
cc=co.merge(fam,on='sample_id',suffixes=('','_fam'));cc.to_csv(O/'cohort83_with_chain_counts.csv',index=False)
summary=dict(family_eval_rows=12252,family_eval_records=3063,all_families=int(fam.component_id.nunique()),multichain_records=int((fam.chains>1).sum()),cohort83_records=len(co),cohort83_families=int(co.family_id.nunique()),cohort83_multichain_records=int((cc.chains>1).sum()),cohort83_multichain_families=int(cc[cc.chains>1].family_id.nunique()),cohort83_by_dataset=cc.groupby('dataset').agg(records=('sample_id','size'),families=('family_id','nunique'),multichain_records=('chains',lambda z:int((z>1).sum()))).to_dict('index'),nonTCR616_records=len(records),nonTCR616_PDBs=int(records.pdb_id.nunique()),nonTCR616_families=int(records.component_id.nunique()),annotated_chains=len(ch),short_peptide_candidate_chains=len(short),short_peptide_records=int(short.sample_id.nunique()),short_peptide_PDBs=int(short.pdb_id.nunique()),short_peptide_families=int(short.family_id.nunique()),short_peptide_twochain_records=int(short[short.two_chain_unambiguous].sample_id.nunique()),short_peptide_twochain_families=int(short[short.two_chain_unambiguous].family_id.nunique()),exact_receptor_context_groups=len(groups),exact_receptor_context_families=int(ctx[ctx.sequence_sha256.isin(groups.sequence_sha256)].family_id.nunique()),exact_receptor_twochain_context_groups=len(groups2),exact_receptor_twochain_context_families=int(ctx[ctx.sequence_sha256.isin(groups2.sequence_sha256)].family_id.nunique()),enzyme_inhibitor_candidate_records=int(ep.sample_id.nunique()),enzyme_inhibitor_candidate_families=int(ep.family_id.nunique()),caveats=['All categories are membership-first feasibility screens, without model-quality selection.','Short chains may be native peptides or crystallographic constructs; exact chain-to-chain contacts in multichain complexes still need confirmation. CA 8A contact label only.','Different partner identity signature may reflect additions in a multicomponent complex; conservative two-chain subset is separately saved.','family_id is the all-chain connected-component split unit, not a biological family name. Antibody/immune connected supercomponent can merge biologically different enzymes.','83 uniformly labelled window cohort does not equal all 3063 original-label family OOF records.','The copy named all_616_record_screen.csv inside case_selection_audit has only 3 records; actual616 file is experiments/non_tcr_case_screen_20261025/all_616_record_screen.csv.'])
(O/'inventory_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2));print('EXACT TWOCHAIN\n',groups2.to_string(index=False));print('SHORT EXAMPLES\n',short[['sample_id','pdb_id','short_description','receptor_descriptions']].head(20).to_string(index=False))
