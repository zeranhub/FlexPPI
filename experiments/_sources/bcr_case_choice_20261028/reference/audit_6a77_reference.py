from pathlib import Path
import sys,types,json,hashlib,collections
import numpy as np
R=Path(r'E:\Project\FlexPPI\experiments\_sources\bcr_case_choice_20261028');P=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline\training_payload')
for name,path in [('flexppi',P/'src/flexppi'),('flexppi.data',P/'src/flexppi/data')]:
 m=types.ModuleType(name);m.__path__=[str(path)];sys.modules[name]=m
from flexppi.data.structures import parse_pdb
from flexppi.data.records import StateRecord
from flexppi.data.labels import compute_dynarepo_labels
import MDAnalysis as mda
D=Path(r'D:\data\dynarepo_non_nucleic_full\A00LZ');raw=D/'structure.pdb';canonical=P/'dataset/structures/A00LZ.pdb';cached=P/'dataset/labels/A00LZ.npz';legacy=Path(r'E:\Project\FlexPPI\artifacts\data\dynarepo_labels\A00LZ.npz')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
st=parse_pdb(str(canonical));oldst=parse_pdb(str(raw));meta=json.loads((D/'dynarepo_metadata.json').read_text())
def atoms(path):
 rows=[]
 for line in Path(path).read_text().splitlines():
  if line.startswith('ATOM  '):rows.append((line[21],line[22:27].strip(),line[17:20],line[12:16].strip(),[float(line[30:38]),float(line[38:46]),float(line[46:54])]))
 return rows
ra,ca=atoms(raw),atoms(canonical);assert len(ra)==len(ca)
assert all(a[:2]==b[:2] for a,b in zip(ra,ca));maxcoords=float(np.max(np.abs(np.array([a[4] for a in ra])-np.array([a[4] for a in ca]))))
name_changes=collections.Counter((a[2].strip()+' '+a[3],b[2].strip()+' '+b[3]) for a,b in zip(ra,ca) if a[2:4]!=b[2:4])
u=mda.Universe(str(canonical),str(D/'trajectory.xtc'));n=len(u.trajectory);start=float(u.trajectory[0].time);end=float(u.trajectory[-1].time);dt=float(u.trajectory.dt)
print('loaded',st.num_residues,'atoms',len(u.atoms),'frames',n,'times',start,end,'dt',dt,'positionsdiff',maxcoords,flush=True)
out=R/'reference/A00LZ_canonical_recomputed.npz'
record=StateRecord(sample_id='A00LZ',dataset='dynarepo',structure_path=str(canonical),trajectory_paths=[str(D/'trajectory.xtc')],temperature=300.,window_ns=100.,target_chains=st.unique_chains(),source_group='dynabench_complex')
result=compute_dynarepo_labels(record,str(out),window_frames=1000,frame_stride=1,include_sasa=False,max_windows=None)
a=np.load(cached);b=np.load(out);o=np.load(legacy)
assert a['values'].shape==(519,3) and b['values'].shape==(519,3)
valuesdiff=[float(np.nanmax(np.abs(a['values'][:,j]-b['values'][:,j]))) if np.any(a['mask'][:,j]) else None for j in range(3)]
vardiff=[float(np.nanmax(np.abs(a['variances'][:,j]-b['variances'][:,j]))) if np.any(a['mask'][:,j]) else None for j in range(3)]
report={'status':'complete','prediction_blind':True,'sample':'A00LZ','pdb':'6A77','chain_roles':{'A':'Robo1 antigen Ig5 domain','B':'antibody light Fab chain','C':'antibody heavy Fab chain'},'canonical_chain_lengths':{c:len(st.indices_for_chain(c)) for c in st.unique_chains()},'raw_parsed_chain_lengths':{c:len(oldst.indices_for_chain(c)) for c in oldst.unique_chains()},'canonical_chains':[{'chain':c,'sequence':st.chain_sequence(c),'authors':[st.author_residue_ids[i] for i in st.indices_for_chain(c)]} for c in st.unique_chains()],
'raw_and_canonical_atom_count':len(ra),'max_atom_coordinate_difference':maxcoords,'atom_name_changes':[{ 'from':a,'to':b,'atoms':n} for (a,b),n in name_changes.items()],
'metadata_declared_replicates':meta['mds'],'local_xtc_files':[p.name for p in D.glob('*.xtc')],'frames':n,'xtc_recorded_time_ps_start':start,'xtc_recorded_time_ps_end':end,'xtc_recorded_dt_ps':dt,'xtc_recorded_span_ns':(end-start)/1000.,'timing_caution':'Header times span 0.999 ns if read literally, while metadata declares 100 ns per replica. Without original frame mapping, actual physical sampling duration is not independently established by local files; do not claim three 100 ns independent replicas.',
'recompute':result,'canonical_label_shape':list(a['values'].shape),'canonical_valid_counts':a['mask'].sum(axis=0).tolist(),'legacy_shape':list(o['values'].shape),'legacy_valid_counts':o['mask'].sum(axis=0).tolist(),'masks_identical_after_recompute':bool(np.array_equal(a['mask'],b['mask'])),'max_abs_label_difference_by_metric':dict(zip(a['metric_names'].tolist(),valuesdiff)),'max_abs_label_variance_difference_by_metric':dict(zip(a['metric_names'].tolist(),vardiff)),
'reference_scope':'Absolute-state MD labels for one Fab-antigen complex, not complete membrane BCR; not paired apo/bound or a mutation-response reference. No new MD performed.',
'paths_and_sha256':{str(p):sha(p) for p in [raw,canonical,D/'trajectory.xtc',cached,legacy,out,P/'src/flexppi/data/labels.py',P/'src/flexppi/constants.py']},'runtime':{'python':sys.version,'numpy':np.__version__,'MDAnalysis':mda.__version__}}
(R/'audits/6a77_reference_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:report[k] for k in ['status','recompute','masks_identical_after_recompute','max_abs_label_difference_by_metric','xtc_recorded_span_ns']},indent=2),flush=True)
