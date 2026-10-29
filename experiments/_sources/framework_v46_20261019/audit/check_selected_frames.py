"""Independent reconstruction of the two actual PDB-backed frames drawn in b."""
from pathlib import Path
import json
import hashlib
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
S=ROOT/'figure1/source_archive'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
data=json.loads((S/'local_frame_candidates.json').read_text(encoding='utf-8'))
assert data['source_sha256']==sha(S/'itrace_4g9f.pdb')
frames={x['residue_number']:x for x in data['all_204_local_frames']}
atoms={}
for line in (S/'itrace_4g9f.pdb').read_text().splitlines():
    if line.startswith('ENDMDL'):break
    if not line.startswith('ATOM  ') or line[21]!='D':continue
    atom=line[12:16].strip()
    if atom not in {'N','CA','C'}:continue
    number=int(line[22:26])
    atoms[(number,atom)]=np.array([float(line[k:k+8]) for k in (30,38,46)],dtype=np.float32).astype(np.float64)
out=[]
unit=lambda a:a/np.linalg.norm(a)
for number in [67,166]:
    f=frames[number]
    n,ca,c=[atoms[number,k] for k in ['N','CA','C']]
    for atom,arr in [('N',n),('CA',ca),('C',c)]:
        assert np.array_equal(arr,np.array(f[atom+'_world_A']))
    x=unit(c-ca)
    y=unit((n-ca)-np.dot(n-ca,x)*x)
    z=unit(np.cross(x,y))
    R=np.column_stack((x,y,z))
    recorded=np.array(f['R_columns_xyz'])
    err=float(np.max(np.abs(R-recorded)))
    assert err<2e-12
    assert abs(np.linalg.det(R)-1)<1e-12
    assert np.max(np.abs(R.T@R-np.eye(3)))<1e-12
    local=R.T@(n-ca)
    assert abs(local[2])<1e-12 and local[1]>0
    probes=np.array([[.2,1.1,-.3],[2.4,.6,.8],[-.7,-.2,1.9]])
    global_points=(R@probes.T).T+ca
    recovered=(R.T@(global_points-ca).T).T
    rt=float(np.max(np.abs(probes-recovered)))
    assert rt<1e-12
    out.append({'chain':'D','residue':number,'residue_name':f['residue_name'],
                'origin_CA_matches_real_PDB':True,'axes_reconstructed_from_real_N_CA_C':True,
                'maximum_R_error':err,'determinant':float(np.linalg.det(R)),
                'transform_round_trip_max_error_A':rt,'N_local_coordinates_A':local.tolist()})
report={'status':'passed','source_PDB':str(S/'itrace_4g9f.pdb'),'source_sha256':sha(S/'itrace_4g9f.pdb'),
        'selected_frames':out,'frame_definition':'CA origin, x toward C; y orthogonalized toward N; z=x cross y; R has basis columns.',
        'coordinate_units':'Angstrom','physical_backbone_atom_input':True}
(ROOT/'audit/selected_frame_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','selected_frames':[67,166],'maximum_R_error':max(x['maximum_R_error'] for x in out)}))
