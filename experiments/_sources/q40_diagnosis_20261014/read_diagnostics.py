from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parent
s=ROOT/'structure';t=ROOT/'trajectory'
p=pd.read_csv(s/'polar_atom_pair_persistence.csv')
cols=['grid','site','source_atom_name','target_chain','target_resid','target_resname','target_atom_name','direct_sequence_neighbor','opposite_N_O_elements','contact_fraction','mean_distance_A']
print('Q40 polar atom pairs all5000')
print(p[(p.site=='Q40')&(p.grid=='original5000')].sort_values('contact_fraction',ascending=False)[cols].head(12).to_string(index=False))
print('Q40 chemically opposite N/O contacts bywindow, contact fraction >= .02')
q=p[(p.site=='Q40')&p.grid.str.startswith('window')&p.opposite_N_O_elements]
print(q.pivot_table(index=['source_atom_name','target_chain','target_resid','target_resname','target_atom_name'],columns='grid',values='contact_fraction').query('window1>=.02 or window2>=.02 or window3>=.02 or window4>=.02 or window5>=.02').to_string())
r=pd.read_csv(s/'static_backbone_neighborhood.csv');print('static backbone neighborhoods');print(r.to_string(index=False))
e=pd.read_csv(t/'original120_entropy_by_window_and_scope.csv')
print('Q40 entropy scopes');print(e[e.residue=='Q40'][['scope','joint_Hchi_nats','Hchi1_marginal_nats','Hchi2_marginal_nats','Hchi3_marginal_nats','total_correlation_nats','joint_states_cover90','largest_joint_state_probability']].to_string(index=False))
