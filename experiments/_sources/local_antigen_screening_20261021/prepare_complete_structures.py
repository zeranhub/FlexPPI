"""Rebuild missing heavy atoms while preserving every observed backbone coordinate.

All candidates and parents use one fixed algorithm and seed. This is a static
side-chain completion experiment, not a prediction of a mutant backbone or MD.
"""
from pathlib import Path
import sys,os,json,hashlib,time,argparse,traceback
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'tool_runtime'))
sys.path.insert(0,'E:/Project/FlexPPI/src')
os.environ['OPENMM_CPU_THREADS']='2'
_handles=[]
if os.name=='nt':_handles.append(os.add_dll_directory(str(Path(sys.prefix)/'Library/bin')))
import numpy as np
import openmm
from openmm import app,unit
from pdbfixer import PDBFixer
from flexppi.data.structures import parse_pdb
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf8')

def main(limit=None):
    manifest=ROOT/'structures/fixed_backbone_manifest.json'
    d=json.loads(manifest.read_text(encoding='utf8'))
    rows=d['candidates']
    if limit:rows=rows[:limit]
    outroot=ROOT/'structures/static_complete';outroot.mkdir(parents=True,exist_ok=True)
    outmanifest=ROOT/'structures/static_complete_manifest.json'
    result={'status':'running','mode':'fixed_backbone_with_rebuilt_heavy_sidechains',
        'method':'PDBFixer1.12.0 addMissingAtoms seed17; OpenMM8.6.1 CPU; original atoms mass0; no missing-residue loops, no hydrogens, no global relaxation, no MD.',
        'input_manifest_sha256':sha(manifest),'preparation_script_sha256':sha(__file__),
        'new_model_outputs_seen_at_method_selection':False,'candidates':[]}
    for row in rows:
        t=time.perf_counter()
        cid=row['candidate_id'];folder=outroot/cid;folder.mkdir(exist_ok=True)
        # The input may already be a prepared PDB or may refer to the template
        # plus requested sequence. The prepare-only audit holds the actual path.
        prepared=ROOT/'structures/fixed_backbone'/cid/'model_input_fixed_backbone.pdb'
        assert prepared.exists(),prepared
        output=folder/'complete_heavy_atoms.pdb'
        auditpath=folder/'preparation_audit.json'
        if output.exists() and auditpath.exists():
            old=json.loads(auditpath.read_text(encoding='utf8'))
            if old.get('status')=='passed' and old.get('input_sha256')==sha(prepared) and old.get('output_sha256')==sha(output):
                result['candidates'].append({**row,'input_pdb':str(output),'static_pdb':str(output),'input_mode':'fixed-backbone heavy-atom completion','chain_sequences':{},'preparation_audit':str(auditpath)})
                continue
        before=parse_pdb(str(prepared))
        fixer=PDBFixer(filename=str(prepared),platform=openmm.Platform.getPlatformByName('CPU'))
        fixer.removeHeterogens(False)
        fixer.missingResidues={}
        fixer.findMissingAtoms()
        missing=[{'chain':r.chain.id,'residue_id':r.id,'name':r.name,'atoms':[a.name for a in atoms]} for r,atoms in fixer.missingAtoms.items()]
        terminal=[{'chain':r.chain.id,'residue_id':r.id,'atoms':atoms} for r,atoms in fixer.missingTerminals.items()]
        fixer.addMissingAtoms(seed=17)
        with output.open('w',encoding='utf8') as f:
            f.write('REMARK 900 SIDECHAINS REBUILT AT FIXED OBSERVED BACKBONE; NOT RELAXED BACKBONE\n')
            app.PDBFile.writeFile(fixer.topology,fixer.positions,f,keepIds=True)
        after=parse_pdb(str(output))
        assert before.sequence==after.sequence
        assert before.chain_ids==after.chain_ids
        assert before.author_residue_ids==after.author_residue_ids
        assert np.array_equal(before.backbone_mask,after.backbone_mask),'Backbone mask changed'
        assert np.array_equal(before.backbone,after.backbone),'Backbone coordinate changed'
        audit={'status':'passed','input':str(prepared),'input_sha256':sha(prepared),'output':str(output),'output_sha256':sha(output),
            'residues':after.num_residues,'new_sidechain_atoms':missing,'new_terminal_atoms':terminal,
            'all_observed_backbone_values_bit_identical':True,'sidechain_rebuilt':True,'backbone_rebuilt':False,
            'OpenMM':openmm.__version__,'platform':'CPU','seed':17,'MD':False,'elapsed_seconds':time.perf_counter()-t}
        write(auditpath,audit)
        result['candidates'].append({**row,'input_pdb':str(output),'static_pdb':str(output),'input_mode':'fixed-backbone heavy-atom completion','chain_sequences':{},'preparation_audit':str(auditpath)})
        write(outmanifest,result)
        print(json.dumps({'id':cid,'status':'passed','seconds':round(audit['elapsed_seconds'],2),'n_rebuilt_residues':len(missing)}),flush=True)
    result['status']='passed'
    write(outmanifest,result)
    print('completed',len(result['candidates']),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int);a=p.parse_args()
    main(a.limit)
