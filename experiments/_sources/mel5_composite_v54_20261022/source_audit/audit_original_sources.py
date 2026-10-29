"""Independent read-only verification for user-requested composite Figure 4."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser

OUT=Path(__file__).resolve().parent
EXP=OUT.parent.parent
INV=EXP/'opposite_affinity_rmsf_cases_20261022'/'inventory'
OLD=EXP/'mel5_joint_case_v49_20261019'/'figure4'
NPZ=EXP/'tcr_case_v40_20261017'/'family_predictions'/'full'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda n:pd.read_csv(INV/n, float_precision='round_trip')
states=['3hg1','4jff','4jfd','4jfe']
pairs=[('3hg1','4jff'),('4jff','4jfd'),('4jff','4jfe')]
frames=read('MEL5_four_states_all_residues.csv')
pep=read('MEL5_four_states_complete_peptide.csv')
cdr=read('MEL5_three_contrasts_all_six_CDR.csv')
iface=read('MEL5_three_contrasts_four_state_union_pMHC_interface.csv')
elig=read('MEL5_four_state_union_all30_Hchi_eligibility.csv')
old_cdr=pd.read_csv(OLD/'source_plotted_CDR_RMSF.csv',float_precision='round_trip')
old_h=pd.read_csv(OLD/'source_plotted_interface_Hchi.csv',float_precision='round_trip')
checks=[]
geometry_checks=[]
parser=PDBParser(QUIET=True)
for p in states:
    q=frames[frames.pdb.eq(p)].sort_values('index').reset_index(drop=True)
    z=np.load(NPZ/f'itrace_{p}.npz',allow_pickle=False)
    assert np.allclose(q.rmsf_full,z['mean'][:,0],rtol=0,atol=1e-14)
    assert np.allclose(q.hchi_full,z['mean'][:,1],rtol=0,atol=1e-14)
    assert np.array_equal(q.chain.astype(str),z['chain'].astype(str))
    assert np.array_equal(q.id.astype(str),z['author_residue_id'].astype(str))
    ps=pep[pep.pdb.eq(p)].sort_values('id')
    assert len(ps)==10 and ps.rmsf_valid.all() and set(ps.id)==set(range(1,11))
    assert np.array_equal(ps.rmsf_full, q[q.role.eq('peptide')].sort_values('id').rmsf_full)
    checks.append({'pdb':p,'npz_sha256':sha(NPZ/f'itrace_{p}.npz'),'n_residues':len(q),'all_rmsf_and_hchi_match_npz':True,'peptide_sequence':''.join(ps.aa),'complete_peptide_n':len(ps)})
    structure=Path(f'E:/data/ITrace/trajectories/{p}_run1_topology.pdb')
    model=parser.get_structure(p,structure)[0]
    coordinates={}
    for r in q.itertuples():
        hits=[v for v in model[r.chain] if str(v.id[1])+str(v.id[2]).strip()==str(r.id)]
        assert len(hits)==1 and 'CA' in hits[0]
        coordinates[r.Index]=hits[0]['CA'].coord.astype(float)
    tcrcoords=np.array([coordinates[r.Index] for r in q[q.role.isin(['tcr_alpha','tcr_beta'])].itertuples()])
    calculated=[]
    for r in q[q.role.isin(['mhc','peptide'])].itertuples():
        distance=float(np.linalg.norm(tcrcoords-coordinates[r.Index],axis=1).min())
        assert (distance<=8)==bool(r.interface),(p,r.chain,r.id,distance,r.interface)
        if distance<=8:calculated.append((r.chain,int(r.id)))
    geometry_checks.append({'pdb':p,'topology_sha256':sha(structure),'all_pMHC_CA8_interface_flags_recomputed':True,'interface_positions':calculated})

groups=['CDR1α','CDR2α','CDR3α','CDR1β','CDR2β','CDR3β']
pair_checks=[]
ref_map={}
ref=frames[frames.pdb.eq('3hg1')]
for p in states:
    g=frames[frames.pdb.eq(p)]
    for role in ['tcr_alpha','tcr_beta']:
        a=ref[ref.role.eq(role)].reset_index(drop=True)
        b=g[g.role.eq(role)].reset_index(drop=True)
        al=PairwiseAligner();al.mode='global';al.match_score=2;al.mismatch_score=-1;al.open_gap_score=-6;al.extend_gap_score=-.5
        al.target_end_gap_score=0;al.query_end_gap_score=0
        hit=al.align(''.join(a.aa),''.join(b.aa))[0]
        for i,j in zip(*hit.indices):
            if i>=0 and j>=0:
                ref_map[(p,b.iloc[j].chain,int(b.iloc[j].id))]=(a.iloc[i].chain,int(a.iloc[i].id),a.iloc[i].region)

cdrrange=[];hrange=[]
for i,(a,b) in enumerate(pairs):
    x=cdr[cdr.pdb_before.eq(a)&cdr.pdb_after.eq(b)]
    h=iface[iface.pdb_before.eq(a)&iface.pdb_after.eq(b)]
    assert len(x)==45 and len(h)==30 and x.rmsf_valid.all() and x.paired.all()
    for table,metric in [(x,'rmsf'),(h,'hchi')]:
        for r in table.itertuples():
            qa=frames[frames.pdb.eq(a)&frames.chain.eq(r.chain_before)&frames.id.eq(r.id_before)].iloc[0]
            qb=frames[frames.pdb.eq(b)&frames.chain.eq(r.chain_after)&frames.id.eq(r.id_after)].iloc[0]
            assert abs(getattr(r,metric+'_full_before')-qa[metric+'_full'])<1e-14
            assert abs(getattr(r,metric+'_full_after')-qb[metric+'_full'])<1e-14
            assert abs(getattr(r,metric+'_change_after_minus_before')-(qb[metric+'_full']-qa[metric+'_full']))<1e-14
    common=[ref_map[(a,r.chain_before,int(r.id_before))] for r in x.itertuples()]
    assert len(set((r[0],r[1]) for r in common))==45
    counts={g:sum(k[2]==g for k in common) for g in groups}
    assert list(counts.values())==[6,6,9,6,5,13]
    e=h.hchi_valid & h.same_aa
    invalid=h[~e]
    eligible=h[e]
    cdrrange.extend(x.rmsf_change_after_minus_before.tolist());hrange.extend(eligible.hchi_change_after_minus_before.tolist())
    pair_checks.append({'before':a,'after':b,'delta_definition':'after minus before','CDR_n':len(x),'CDR_region_counts':counts,'interface_n':len(h),'interface_Hchi_eligible_n':int(e.sum()),'interface_Hchi_ineligible':[{'chain':r.chain_before,'position':int(r.id_before),'before_aa':r.aa_before,'after_aa':r.aa_after,'both_valid':bool(r.hchi_valid),'same_aa':bool(r.same_aa)} for r in invalid.itertuples()],'Q72_Hchi_delta':float(h[h.chain_before.eq('A')&h.id_before.eq(72)].hchi_change_after_minus_before.iloc[0])})

# The ranked Hchi values and absolute CDR profiles use precisely the same maturation contrast.
m=cdr[cdr.pdb_before.eq('3hg1')&cdr.pdb_after.eq('4jff')]
for r in old_cdr.itertuples():
    chain='D' if 'α' in r.region else 'E'
    q=m[m.chain_before.eq(chain)&m.id_before.eq(r.common_TCR_position)]
    assert len(q)==1
    q=q.iloc[0]
    assert q.id_after==r.derivative_input_position
    assert abs(q.rmsf_full_before-r.MEL5_model_RMSF_A)<1e-14
    assert abs(q.rmsf_full_after-r.alpha24beta17_model_RMSF_A)<1e-14
ih=iface[iface.pdb_before.eq('3hg1')&iface.pdb_after.eq('4jff')]
assert len(old_h)==23
ranked=[]
for r in old_h.itertuples():
    q=ih[ih.chain_before.eq(r.chain_a)&ih.id_before.eq(r.position_a)]
    assert len(q)==1
    q=q.iloc[0]
    assert q.hchi_valid and q.same_aa
    assert abs(q.hchi_full_before-r.hchi_full_a)<1e-14
    assert abs(q.hchi_full_after-r.hchi_full_b)<1e-14
    ranked.append({'chain':r.chain_a,'position':int(r.position_a),'AA':r.aa_a,'delta_Hchi_nats':r.hchi_full_b-r.hchi_full_a})
ranked=sorted(ranked,key=lambda x:x['delta_Hchi_nats'])
assert ranked[0]['chain']=='A' and ranked[0]['position']==72
same_all=elig[elig.same_AA_valid_Hchi_all_four]
assert len(same_all)==22
assert [p['interface_Hchi_eligible_n'] for p in pair_checks]==[23,23,22]
assert max(abs(v) for v in cdrrange)<.45 and max(abs(v) for v in hrange)<.45
sources=[INV/n for n in ['MEL5_four_states_all_residues.csv','MEL5_four_states_complete_peptide.csv','MEL5_three_contrasts_all_six_CDR.csv','MEL5_three_contrasts_four_state_union_pMHC_interface.csv','MEL5_four_state_union_all30_Hchi_eligibility.csv','MEL5_four_states_sequence_and_KD.csv']]
sources += [OLD/'source_plotted_CDR_RMSF.csv',OLD/'source_plotted_interface_Hchi.csv',OLD/'4jff_overview.png',OLD/'4jff_Q72_F53_G54.png',OLD/'plot.py',INV.parent/'plot_case_series.py']
audit={'status':'passed','scope':'Read-only source and numeric audit for composite Figure 4; no manuscript or plot mutation.','model':{'training_records':2430,'ATLAS':1866,'DynaRepo':564,'ITrace':0,'checkpoint_sha256':'ff2974406d5d1ee07e6dc7118e43e57c0acc9849cc8feba7f5290476fca15ec5'},'new_training_inference_MD':False,'NPZ_checks':checks,'contrasts':pair_checks,'old_vs_new_sources_exact':True,'ranked_Hchi_sites':ranked,'Hchi_decreased_maturation':sum(r['delta_Hchi_nats']<0 for r in ranked),'fixed_all_four_Hchi_same_AA_n':22,'data_ranges':{'CDR_delta_RMSF_A':[min(cdrrange),max(cdrrange)],'Hchi_delta_nats':[min(hrange),max(hrange)],'heatmaps_within_plus_minus_0_45':True},'interface_definition':'Union across four states of pMHC residues with nearest TCR C-alpha distance <=8 Angstrom; no atomic contact or MD occupancy implied.','missingness':'Seven maturation/G4A positions are A or G without side-chain chi angles. I7A additionally changes peptide I7 to Ala; numeric network outputs at invalid sites are never interpretable Hchi.','reference_numbering':'All TCR columns mapped by sequence alignment to MEL5 3HG1. Engineered TCR author/input IDs are +1 at highlighted beta53/54, but x-axis remains MEL5 numbering. All45 CDR positions included; gaps between CDR segments are different loops, not chain breaks.','sources':[{'path':str(p),'sha256':sha(p)} for p in sources]}
audit['direct_geometry_checks']=geometry_checks
(OUT/'independent_source_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
pd.DataFrame(ranked).to_csv(OUT/'ranked23_Hchi_verified.csv',index=False,encoding='utf-8-sig')
print(json.dumps({'status':'passed','output':str(OUT/'independent_source_audit.json'),'contrasts':[(x['before'],x['after'],x['interface_Hchi_eligible_n']) for x in pair_checks],'ranges':audit['data_ranges']},ensure_ascii=False))
