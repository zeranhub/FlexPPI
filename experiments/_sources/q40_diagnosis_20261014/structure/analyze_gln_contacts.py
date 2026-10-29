"""Read-only six-GLN structural and existing-MD contact diagnosis. CPU only."""
from __future__ import annotations
import sys
sys.dont_write_bytecode=True
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
os.environ['CUDA_VISIBLE_DEVICES']=''
from pathlib import Path
from datetime import datetime, timezone
import json, hashlib, tempfile, shutil, time, warnings, copy
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
RUNTIME=BASE/'label_sensitivity/runtime'
sys.path.insert(0,str(RUNTIME))
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

SOURCE_PDB=Path('D:/data/dynarepo_non_nucleic_full/A003R/structure.pdb')
SOURCE_XTC=Path('D:/data/dynarepo_non_nucleic_full/A003R/trajectory.xtc')
PROVENANCE=Path('E:/Project/FlexPPI/experiments/_case_validation_sources/ubiquitin/provenance.json')
OLD_CONTACT=PROVENANCE.parent/'reference_environment_contacts.csv'
OLD_SASA=PROVENANCE.parent/'reference_solvent_accessibility.csv'
SITES=(2,31,40,41,49,62)
SIDE_NAMES=('CB','CG','CD','OE1','NE2')
BACK_NAMES=('N','CA','C','O')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean(x) for x in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value

def js(path,value):
    Path(path).write_text(json.dumps(clean(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def check_sources(sources):
    for path,digest in sources.items():assert sha(path)==digest,path

def freeze():
    HERE.mkdir(parents=True,exist_ok=True)
    provenance=json.loads(PROVENANCE.read_text(encoding='utf-8'))
    sources={str(SOURCE_PDB):provenance['source_structure_sha256'],str(SOURCE_XTC):provenance['source_trajectory_sha256'],
             str(PROVENANCE):sha(PROVENANCE),str(OLD_CONTACT):sha(OLD_CONTACT),str(OLD_SASA):sha(OLD_SASA)}
    check_sources(sources)
    plan={
      'protocol':'six_Gln_A003R_structural_contact_diagnosis_v1','created_utc':datetime.now(timezone.utc).isoformat(),
      'analysis_code_sha256':sha(__file__),'sources_sha256':sources,'sites':['Q'+str(s) for s in SITES],
      'source':'A003R full protein complex, chainA ubiquitin72 and chainB UBA47; all6 valid GLN, none omitted',
      'static':'Original inputPDB coordinates; compare sidechainheavy CB/CG/CD/OE1/NE2 and backbone N/CA/C/O to all proteinheavy atoms',
      'contact':'Sidechain-to-partner-residue minimum heavyatomdistance <=4.0angstrom. Exclude self and direct preceding/following residue in same chain by sequence index. Count any/frame and each partner residue/frame.',
      'polar':'OE1/NE2 to other proteinN/O heavyatoms <=3.5angstrom. Exclude selfresidue, retain adjacent residues with explicit flag. All N/O closepairs are distance-only; opposite-elements flag distinguishes O-N from O-O/N-N. Not assigned hydrogen bonds.',
      'primary_grid':{'indices':list(range(0,1001,2)),'frames':501,'range_ns':[0,100],'spacing_ps':200},
      'full_grid':{'indices':list(range(5000)),'frames':5000,'range_ns':[0,499.9],'spacing_ps':100},
      'five_original_windows':'frames0:1000,1000:2000,2000:3000,3000:4000,4000:5000; each1000 observations,0-99.9,...,400-499.9ns. Sourceframe5000 at500ns is not in original labels.',
      'late_change':'All6window4→5 changes retained, no contact chosen because it supports Q40story; topcontacts ranked by full5000frameoccupancy, ties by atom/residue identity.',
      'PBC':'Use trajectory unitcell minimumimage distances when valid positive lengths/angles; report dimensions. StaticPDB uses direct Cartesian distances, no box supplied.',
      'SASA':'Original heavyatomcomplex only, BiopythonShrakeRupley probe1.4angstrom n_points960. SidechainSASA is sum5sideatom areas. Also isolatedsameconformerreference: sidechainarea in isolatedsame residue, not standard maximalASA. No MDaveragedSASA claimed.',
      'scientific_limits':'One existing trajectory, exploratory contactdiagnosis; noindependentrepeats, no causalerrorattribution, no hydrogenbondangles/energies; modelusesbackboneinput, sidechaincontacts are observations rather than model input.',
    }
    js(HERE/'analysis_plan.json',plan)
    return plan

def atom_row(atom,prefix=''):
    return {prefix+'atom_index':int(atom.index),prefix+'chain':str(atom.chainID),prefix+'resid':int(atom.resid),
            prefix+'resname':str(atom.resname),prefix+'atom_name':str(atom.name),prefix+'element':str(atom.element)}

def main():
    sys.stdout.reconfigure(encoding='utf-8')
    started=time.perf_counter()
    plan=freeze()
    source_stats={path:(Path(path).stat().st_size,Path(path).stat().st_mtime_ns) for path in plan['sources_sha256']}
    stage=HERE/'_staging';stage.mkdir(exist_ok=True)
    static_rows,static_polar,backbone_density,frame_rows,atom_mapping=[],[],[],[],[]
    contact_counts=np.zeros((7,6,119),dtype=np.int32)
    contact_min=np.full((7,6,119),np.inf)
    contact_sum=np.zeros((7,6,119),dtype=float)
    polar_counts=polar_min=polar_sum=None
    old_contact_checks=[]
    with tempfile.TemporaryDirectory(prefix='verified_existing_md_',dir=stage) as folder:
        folder=Path(folder)
        assert folder.resolve().is_relative_to(stage.resolve())
        top,traj=folder/'structure.pdb',folder/'trajectory.xtc'
        shutil.copy2(SOURCE_PDB,top);shutil.copy2(SOURCE_XTC,traj)
        assert sha(top)==plan['sources_sha256'][str(SOURCE_PDB)]
        assert sha(traj)==plan['sources_sha256'][str(SOURCE_XTC)]
        u=mda.Universe(str(top),str(traj))
        try:
            assert len(u.trajectory)==5001
            assert len(u.residues)==119
            assert [(str(x),int(y)) for x,y in zip(*np.unique(u.atoms.chainIDs,return_counts=True))]
            residues=u.residues
            chain_res={chain:[r for r in residues if str(r.atoms[0].chainID)==chain] for chain in ['A','B']}
            assert len(chain_res['A'])==72 and len(chain_res['B'])==47
            heavy=u.atoms[~np.isin(np.char.upper(u.atoms.elements.astype(str)),['H','D'])]
            assert len(heavy)==935
            heavy_ri=heavy.resindices.astype(int)
            assert np.array_equal(np.unique(heavy_ri),np.arange(119))
            starts=np.r_[0,np.flatnonzero(np.diff(heavy_ri))+1]
            assert len(starts)==119
            no=heavy[np.isin(heavy.elements,['N','O'])]
            heavy_no_indices=no.indices
            side_indices=[];polar_indices=[];eligible=[];direct_neighbor=[]
            static_positions=u.atoms.positions.copy()
            for si,site in enumerate(SITES):
                rg=u.select_atoms(f'chainID A and resid {site}').residues
                assert len(rg)==1 and rg[0].resname=='GLN'
                r=rg[0];lookup={a.name:int(a.index) for a in r.atoms}
                assert len(lookup)==len(r.atoms)
                assert r.ix==site-1
                sides=[lookup[n] for n in SIDE_NAMES];side_indices.extend(sides)
                ends=[lookup[n] for n in ['OE1','NE2']];polar_indices.extend(ends)
                same_chain_neighbors={x.ix for x in chain_res['A'] if abs(chain_res['A'].index(x)-chain_res['A'].index(r))==1}
                allowed=~np.isin(np.arange(119),[r.ix,*same_chain_neighbors])
                eligible.append(allowed)
                direct_neighbor.append(np.isin(no.resindices,list(same_chain_neighbors)))
                atom_mapping.append({'site':'Q'+str(site),'residue_index':int(r.ix),'side_atoms':[atom_row(u.atoms[i]) for i in sides],
                                     'terminal_polar_atoms':[atom_row(u.atoms[i]) for i in ends],
                                     'excluded_residue_indices_for_4A':[int(r.ix),*sorted(map(int,same_chain_neighbors))]})
                for group,names in [('sidechain',SIDE_NAMES),('backbone',BACK_NAMES)]:
                    ids=[lookup[n] for n in names]
                    distances=distance_array(static_positions[ids],static_positions[heavy.indices])
                    by_res=np.minimum.reduceat(distances,starts,axis=1).min(axis=0)
                    for partner in residues:
                        if partner.ix==r.ix:continue
                        aa=partner.atoms[0]
                        same=str(aa.chainID)=='A'
                        static_rows.append({'site':'Q'+str(site),'source_group':group,'partner_chain':str(aa.chainID),
                          'partner_resid':int(partner.resid),'partner_resname':str(partner.resname),
                          'samechain':same,'direct_sequence_neighbor':bool(partner.ix in same_chain_neighbors),
                          'eligible_nonadjacent_4A':bool(allowed[partner.ix]),'minimum_distance_angstrom':float(by_res[partner.ix]),
                          'within4A_nonadjacent':bool(allowed[partner.ix] and by_res[partner.ix]<=4.)})
                    if group=='sidechain':
                        for oldrow in pd.read_csv(OLD_CONTACT).query("site == @('Q'+str(site))").to_dict('records') if False else []:pass
                pdist=distance_array(static_positions[ends],static_positions[no.indices])
                for ei,source_index in enumerate(ends):
                    for ti,target in enumerate(no):
                        if target.resindex==r.ix:continue
                        if pdist[ei,ti] > 5.:continue
                        static_polar.append({'site':'Q'+str(site),**atom_row(u.atoms[source_index],'source_'),**atom_row(target,'target_'),
                          'samechain':str(target.chainID)=='A','direct_sequence_neighbor':bool(direct_neighbor[-1][ti]),
                          'opposite_N_O_elements':u.atoms[source_index].element!=target.element,
                          'target_backbone':target.name in ['N','CA','C','O','OXT'],
                          'distance_angstrom':float(pdist[ei,ti]),'within3p5A':bool(pdist[ei,ti]<=3.5)})
                ca=lookup['CA'];caatoms=u.select_atoms('name CA')
                cad=distance_array(static_positions[[ca]],static_positions[caatoms.indices])[0]
                valid=allowed[caatoms.resindices]
                backbone_density.append({'site':'Q'+str(site),'source':'original_input_PDB','CA_neighbors8A_samechain_nonadjacent':int((valid&(caatoms.chainIDs=='A')&(cad<=8)).sum()),
                  'CA_neighbors8A_interchain':int((valid&(caatoms.chainIDs=='B')&(cad<=8)).sum()),
                  'nearest_nonself_UBA_CA_A':float(cad[caatoms.chainIDs=='B'].min())})
            eligible=np.array(eligible,bool)
            side_indices=np.array(side_indices,int);polar_indices=np.array(polar_indices,int)
            source_site_res=np.repeat(np.array(SITES)-1,2)
            polar_valid=no.resindices[None,:]!=source_site_res[:,None]
            polar_counts=np.zeros((7,len(polar_indices),len(no)),np.int32)
            polar_min=np.full(polar_counts.shape,np.inf)
            polar_sum=np.zeros(polar_counts.shape,float)
            primary_count=0;boxes=[];max_static_first_difference=None
            for frame,ts in enumerate(u.trajectory[:5000]):
                assert abs(ts.time-frame*100.)<=.05
                assert np.isfinite(ts.positions).all()
                if frame==0:max_static_first_difference=float(np.max(np.abs(ts.positions-static_positions)))
                box=ts.dimensions
                valid_box=box is not None and np.isfinite(box).all() and (box[:3]>0).all() and (box[3:]>0).all()
                if frame in [0,1000,2000,3000,4000,4999]:boxes.append({'frame':frame,'time_ns':float(ts.time/1000),'box':None if box is None else list(map(float,box)),'minimum_image_used':bool(valid_box)})
                dist=distance_array(ts.positions[side_indices],ts.positions[heavy.indices],box=box if valid_box else None)
                per_res=np.minimum.reduceat(dist,starts,axis=1).reshape(6,5,119).min(axis=1)
                per_res[~eligible]=np.inf
                touch=per_res<=4.
                pdist=distance_array(ts.positions[polar_indices],ts.positions[no.indices],box=box if valid_box else None)
                pdist[~polar_valid]=np.inf
                ptouch=pdist<=3.5
                grids=[0,2+frame//1000]
                primary=(frame<=1000 and frame%2==0)
                if primary:grids.append(1);primary_count+=1
                for gi in grids:
                    contact_counts[gi]+=touch
                    contact_min[gi]=np.minimum(contact_min[gi],per_res)
                    contact_sum[gi]+=np.where(eligible,per_res,0.)
                    polar_counts[gi]+=ptouch
                    polar_min[gi]=np.minimum(polar_min[gi],pdist)
                    polar_sum[gi]+=np.where(polar_valid,pdist,0.)
                for si,site in enumerate(SITES):
                    chain_a=np.arange(119)<72;chain_b=~chain_a
                    psel=slice(2*si,2*si+2);pnear=ptouch[psel]
                    no_a=no.chainIDs=='A';no_b=~no_a
                    nonadj=~direct_neighbor[si]
                    frame_rows.append({'site':'Q'+str(site),'frame_zero_based':frame,'time_ns':float(ts.time/1000),'original_window':frame//1000+1,
                      'in_primary501':primary,'nonadjacent_4A_samechain_any':bool(touch[si,chain_a].any()),
                      'nonadjacent_4A_interchain_any':bool(touch[si,chain_b].any()),
                      'nonadjacent_4A_samechain_residue_count':int(touch[si,chain_a].sum()),
                      'nonadjacent_4A_interchain_residue_count':int(touch[si,chain_b].sum()),
                      'nearest_samechain_nonadjacent_heavy_A':float(per_res[si,chain_a].min()),
                      'nearest_interchain_heavy_A':float(per_res[si,chain_b].min()),
                      'polar3p5_samechain_any':bool(pnear[:,no_a].any()),'polar3p5_interchain_any':bool(pnear[:,no_b].any()),
                      'polar3p5_samechain_atom_pair_count':int(pnear[:,no_a].sum()),'polar3p5_interchain_atom_pair_count':int(pnear[:,no_b].sum()),
                      'polar3p5_nonadjacent_samechain_any':bool(pnear[:,no_a&nonadj].any()),
                      'polar3p5_nonadjacent_interchain_any':bool(pnear[:,no_b&nonadj].any())})
                if frame%1000==0:print(json.dumps({'phase':'existing_MD_contacts','processed_frame':frame,'seconds':round(time.perf_counter()-started,1)}),flush=True)
            assert primary_count==501
            grid_names=['original5000','primary501']+['window'+str(i) for i in range(1,6)]
            denominators=[5000,501]+[1000]*5
            persistence=[];polar_persistence=[]
            for gi,grid in enumerate(grid_names):
                for si,site in enumerate(SITES):
                    for target in residues:
                        if not eligible[si,target.ix]:continue
                        persistence.append({'grid':grid,'site':'Q'+str(site),'partner_chain':str(target.atoms[0].chainID),
                          'partner_resid':int(target.resid),'partner_resname':str(target.resname),'samechain':str(target.atoms[0].chainID)=='A',
                          'frames':denominators[gi],'contact_frames':int(contact_counts[gi,si,target.ix]),
                          'contact_fraction':float(contact_counts[gi,si,target.ix]/denominators[gi]),
                          'minimum_distance_A':float(contact_min[gi,si,target.ix]),'mean_minimum_distance_A':float(contact_sum[gi,si,target.ix]/denominators[gi])})
            for ei,source_index in enumerate(polar_indices):
                si=ei//2;source=u.atoms[source_index]
                for ti,target in enumerate(no):
                    if not polar_valid[ei,ti] or polar_counts[0,ei,ti]==0:continue
                    for gi,grid in enumerate(grid_names):
                        polar_persistence.append({'grid':grid,'site':'Q'+str(SITES[si]),**atom_row(source,'source_'),**atom_row(target,'target_'),
                          'samechain':str(target.chainID)=='A','direct_sequence_neighbor':bool(direct_neighbor[si][ti]),
                          'opposite_N_O_elements':source.element!=target.element,'target_backbone':target.name in ['N','CA','C','O','OXT'],
                          'frames':denominators[gi],'contact_frames':int(polar_counts[gi,ei,ti]),'contact_fraction':float(polar_counts[gi,ei,ti]/denominators[gi]),
                          'minimum_distance_A':float(polar_min[gi,ei,ti]),'mean_distance_A':float(polar_sum[gi,ei,ti]/denominators[gi])})
            frame_df=pd.DataFrame(frame_rows)
            summaries=[]
            avgcols=[c for c in frame_df if c.endswith('_any') or c.endswith('_count')]
            for grid in grid_names:
                if grid=='original5000':section=frame_df
                elif grid=='primary501':section=frame_df[frame_df.in_primary501]
                else:section=frame_df[frame_df.original_window.eq(int(grid[-1]))]
                for site,g in section.groupby('site',sort=False):
                    row={'grid':grid,'site':site,'frames':len(g),'first_time_ns':g.time_ns.min(),'last_time_ns':g.time_ns.max()}
                    row.update({c+'_mean':float(g[c].mean()) for c in avgcols})
                    row.update({'nearest_nonadjacent_samechain_heavy_mean_A':float(g.nearest_samechain_nonadjacent_heavy_A.mean()),
                                'nearest_interchain_heavy_mean_A':float(g.nearest_interchain_heavy_A.mean())})
                    summaries.append(row)
            static_df=pd.DataFrame(static_rows)
            old_df=pd.read_csv(OLD_CONTACT)
            for r in old_df.itertuples(index=False):
                v=static_df[static_df.site.eq(r.site)&static_df.source_group.eq('sidechain')&static_df.partner_chain.eq(r.partner_chain)&static_df.partner_resid.eq(r.partner_input_residue_id)]
                assert len(v)==1
                difference=abs(v.iloc[0].minimum_distance_angstrom-r.min_sidechain_to_partner_heavy_atom_A)
                assert difference<.0006
                old_contact_checks.append({'site':r.site,'partner_chain':r.partner_chain,'partner_resid':int(r.partner_input_residue_id),'absolute_rounding_difference_A':float(difference)})
            output_rows={
              'static_contacts_by_residue.csv':static_rows,'static_polar_atom_contacts_within5A.csv':static_polar,
              'static_backbone_neighborhood.csv':backbone_density,'trajectory_contact_timeseries_all5000.csv':frame_rows,
              'trajectory_contact_timeseries_primary501.csv':frame_df[frame_df.in_primary501].to_dict('records'),
              'contact_summary_all6_by_grid_window.csv':summaries,'residue_contact_persistence.csv':persistence,
              'polar_atom_pair_persistence.csv':polar_persistence,
            }
            for name,rows in output_rows.items():pd.DataFrame(rows).to_csv(HERE/name,index=False)
            js(HERE/'atom_mapping.json',atom_mapping)
            js(HERE/'trajectory_read_audit.json',{'source_frames':5001,'source_atoms':len(u.atoms),'protein_heavy_atoms':len(heavy),
               'protein_N_O_heavy_atoms':len(no),'chain_residue_counts':{c:len(v) for c,v in chain_res.items()},
               'frames_read':5000,'primary_frames':primary_count,'time_spacing_ps':100.,'sourceframe5000_at500ns_excluded_from_labels':True,
               'minimum_image_boxes_checked':boxes,'max_sourcePDB_vs_XTCframe0_coordinate_difference_A':max_static_first_difference,
               'all_17_archived_Q41_Q62_static_distance_rows_reproduced':old_contact_checks,'MDAnalysis':mda.__version__})
        finally:u.trajectory.close()
    # Original input SASA only. No PyMOL installation or new package is needed.
    structure=PDBParser(QUIET=True).get_structure('A003R',str(SOURCE_PDB))
    for residue in list(structure.get_residues()):
        for atom in list(residue):
            if str(atom.element).upper() in ['H','D']:residue.detach_child(atom.id)
    assert len(list(structure.get_atoms()))==935
    sr=ShrakeRupley(probe_radius=1.4,n_points=960)
    sr.compute(structure,level='A')
    sasa=[]
    for site in SITES:
        residue=structure[0]['A'][(' ',site,' ')]
        atom_sasa={a.name:float(a.sasa) for a in residue}
        isolated=copy.deepcopy(residue)
        sr.compute(isolated,level='A')
        own=float(sum(a.sasa for a in isolated if a.name in SIDE_NAMES))
        area=float(sum(atom_sasa[n] for n in SIDE_NAMES))
        sasa.append({'site':'Q'+str(site),'source':'original_inputPDB_not_MDaverage','sidechain_SASA_A2':area,
          'whole_residue_SASA_A2':sum(atom_sasa.values()),'isolated_same_conformer_sidechain_SASA_A2':own,
          'fraction_of_isolated_sameconformer_SASA':area/own,'probe_radius_A':1.4,'sphere_points':960,
          'heavy_complex_atoms':935,'software':'BiopythonShrakeRupley','atom_SASA_A2_json':json.dumps(atom_sasa,sort_keys=True)})
    pd.DataFrame(sasa).to_csv(HERE/'static_SASA_all6.csv',index=False)
    check_sources(plan['sources_sha256'])
    assert all((Path(k).stat().st_size,Path(k).stat().st_mtime_ns)==v for k,v in source_stats.items())
    js(HERE/'run_audit.json',{'status':'complete','created_utc':datetime.now(timezone.utc).isoformat(),'elapsed_seconds':time.perf_counter()-started,
       'plan_sha256':sha(HERE/'analysis_plan.json'),'code_sha256':sha(__file__),'sources_sha256':plan['sources_sha256'],
       'all_source_bytes_size_mtime_unchanged':True,'source_trajectory_read_only_staging_copy':True,'all6GLN_preserved':True,
       'no_new_MD_or_model_training':True,'primary501_and_all5000_separate':True,
       'outputs':{f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in HERE.iterdir() if f.is_file() and f.name not in ['run_audit.json','analysis_plan.json']},
       'software':{'MDAnalysis':mda.__version__,'numpy':np.__version__,'pandas':pd.__version__}})
    print(json.dumps({'status':'complete','elapsed_seconds':round(time.perf_counter()-started,1),'directory':str(HERE)},ensure_ascii=False),flush=True)

if __name__=='__main__':main()
