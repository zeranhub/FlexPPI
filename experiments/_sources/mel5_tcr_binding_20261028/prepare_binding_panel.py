from pathlib import Path
import json,csv,hashlib
root=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_tcr_binding_20261028')
cd=json.loads((root/'numbering/verified_CDR_sequences.json').read_text())
original=Path(r'E:\Project\FlexPPI\experiments\_sources\mel5_absolute_affinity_20261028')
with (original/'dynamics/all_four_states_complete_chain_sequences.csv').open(encoding='utf-8-sig',newline='') as f:
 atom={(r['pdb'],r['chain']):r['sequence'] for r in csv.DictReader(f)}
for k,info in cd.items():
 pdb,ch=k.split('_')
 for label in ['CDR1_IMGT27_38','CDR2_IMGT56_65','CDR3_with_C104_F118']:
  assert info[label] in atom[(pdb,ch)], (k,label)
for pdb in ['4jfd','4jfe']:
 for ch in ['D','E']:
  for label in ['CDR1_IMGT27_38','CDR2_IMGT56_65','CDR3_with_C104_F118']:
   assert cd[pdb+'_'+ch][label]==cd['4jff_'+ch][label]
rows=[
 ('3hg1','3hg1','MEL5','ELA','ELAGIGILTV',18000,1,'core_four_states','Table2'),
 ('4jff','4jff','alpha24beta17','ELA','ELAGIGILTV',0.6,1,'core_four_states','Table2'),
 ('4jfd','4jfd','alpha24beta17','G4A','ELAAIGILTV',36000,1,'core_four_states','Table2'),
 ('4jfe','4jfe','alpha24beta17','I7A','ELAGIGALTV',31000,1,'core_four_states','Table2'),
 ('mature_E1A','','alpha24beta17','E1A','ALAGIGILTV',140000,1,'published_alanine_scan','Table2'),
 ('mature_I5A','','alpha24beta17','I5A','ELAGAGILTV','',0,'published_alanine_scan','Table2_NB'),
 ('mature_G6A','','alpha24beta17','G6A','ELAGIAILTV',41000,1,'published_alanine_scan','Table2'),
 ('mature_L8A','','alpha24beta17','L8A','ELAGIGIATV',21000,1,'published_alanine_scan','Table2'),
 ('mature_T9A','','alpha24beta17','T9A','ELAGIGILAV',37000,1,'published_alanine_scan','Table2'),
 ('mature_control_ILAKFLHWL','','alpha24beta17','negative_control','ILAKFLHWL','',0,'published_negative_control','Results_SPR_controls_data_not_shown'),
 ('mature_control_SLLMWITQC','','alpha24beta17','negative_control','SLLMWITQC','',0,'published_negative_control','Results_SPR_controls_data_not_shown'),
 ('mature_control_YLEPGPVTA','','alpha24beta17','negative_control','YLEPGPVTA','',0,'published_negative_control','Results_SPR_controls_data_not_shown')]
result=[]
for case,pdb,tcr,mut,pep,kd,binding,group,evidence in rows:
 ref='3hg1' if tcr=='MEL5' else '4jff'
 a,b=cd[ref+'_D'],cd[ref+'_E']
 r=dict(case_id=case,pdb=pdb,tcr=tcr,peptide_label=mut,peptide=pep,HLA='HLA-A*02:01',experimental_KD_nM=kd,experimental_binding=binding,panel_group=group,evidence=evidence,CDR3_alpha=a['CDR3_with_C104_F118'],CDR3_beta=b['CDR3_with_C104_F118'],TRAV='TRAV12-2',TRAJ='',TRBV='TRBV30',TRBJ='',T_cell_type='CD8',CDR_source_pdb=ref)
 for prefix,vals in [('A',a),('B',b)]:
  r[prefix+'1']=vals['CDR1_IMGT27_38'];r[prefix+'2']=vals['CDR2_IMGT56_65'];r[prefix+'3']=vals['CDR3_IMGT105_117'];r[prefix+'3_with_anchors']=vals['CDR3_with_C104_F118']
 result.append(r)
with (root/'binding_panel_inputs.csv').open('w',encoding='utf-8-sig',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
(root/'binding_panel_inputs.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
audit=dict(input_sha256=hashlib.sha256((root/'binding_panel_inputs.csv').read_bytes()).hexdigest(),numbering='ANARCII2.0.8 TCR accuracy CPU, original RCSB SEQRES',CDRs_also_found_in_processed_atom_sequences=True,mature_CDRs_identical_in_all_three_crystals=True,V_gene_source='TCR3d classI metadata matched all four PDBs: TRAV12-2, TRBV30',J_genes='Not reliably annotated; left blank for official UNK handling instead of guessed names',experimental_source='Madura et al JBC2013 doi10.1074/jbc.M113.464560 Table2 and Results paragraph on SPR negative controls',label_definition='Detected SPR binding vs reported NB/no binding at original assay conditions. NB is not mathematical zero affinity.',panel_prespecified_before_model_results=True,rows=12,positive=8,negative=4,core_four_states_all_positive=True,evaluation='Descriptive within-case scores, uncalibrated fixed 0.5 reference classification only if no model-specific cutoff, mutation rank vs ELA separate from binding label. No custom KD cutoff or score calibration.',training_overlap='Audited separately per model; native MEL5+ELA known exact overlap for NetTCR2.2/ImRex. NetTCR2.2 known-overlap-excluding partition0 CV model summary fixed before prediction, all20 output retained.')
(root/'panel_design_and_numbering_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
