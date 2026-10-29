"""Local identity/sequence audit for new multimode panels, without inference."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import argparse,hashlib,json,sys

HERE=Path(__file__).resolve().parent
BASE=Path(r'E:\Project\FlexPPI\experiments\_sources')
REF=Path(r'E:\Project\FlexPPI\experiments\core_corrected_pipeline')
sys.path.insert(0,r'E:\Project\FlexPPI\src')
from flexppi.data.structures import parse_pdb
from Bio.Align import PairwiseAligner
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text(encoding='utf8'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',required=True)
    parser.add_argument('--output',default=str(HERE/'new_inputs_training_overlap.json'))
    args=parser.parse_args()
    trainpath=REF/'strict_external_payload/splits/strict_external_train.jsonl'
    assert sha(trainpath)=='e0c84c566dea09b436dd7d07e9e6fd8ce1a8cf0c72be20931179e6689c8b44bd'
    records=[json.loads(l) for l in trainpath.read_text().splitlines() if l.strip()]
    trainids={r['sample_id'] for r in records};assert len(trainids)==2431
    inventory=read(REF/'sequence_audit/chain_inventory.json')
    trainseq=defaultdict(list);allseq=defaultdict(list)
    for c in inventory:
        allseq[c['sequence']].append(c)
        if c['sample_id'] in trainids:trainseq[c['sequence']].append(c)
    groups=read(BASE/'family_folds/family_components.json')
    recordgroup={s:g for g in groups for s in g['members']}
    data=read(args.manifest)
    templates={t['template_id']:t for t in data['templates']}
    candidates=data.get('candidates',[])
    aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
    cache={}
    def sequence_audit(seq):
        if seq in cache:return cache[seq]
        hits=[]
        if len(seq)>=20:
            for target,cc in trainseq.items():
                if min(len(seq),len(target))/max(len(seq),len(target))<.7:continue
                if seq==target:
                    hits.append(dict(identity=1.,query_coverage=1.,target_coverage=1.,training_chains=[dict(sample_id=c['sample_id'],chain=c['chain'],system=c['biological_system']) for c in cc]));continue
                a=aligner.align(seq,target)[0];qb,tb=a.aligned
                if not len(qb):continue
                qc=int(qb[-1,1]-qb[0,0])/len(seq);tc=int(tb[-1,1]-tb[0,0])/len(target)
                if min(qc,tc)<.7:continue
                identical=sum(sum(aa==bb for aa,bb in zip(seq[int(qs):int(qe)],target[int(ts):int(te)])) for (qs,qe),(ts,te) in zip(qb,tb))
                x,y=a.coordinates;cols=sum(max(int(x[i+1]-x[i]),int(y[i+1]-y[i])) for i in range(len(x)-1))
                identity=identical/cols
                if identity>=.3:hits.append(dict(identity=identity,query_coverage=qc,target_coverage=tc,training_chains=[dict(sample_id=c['sample_id'],chain=c['chain'],system=c['biological_system']) for c in cc]))
        inherited={}
        for c in allseq.get(seq,[]):
            if c['sample_id'] in recordgroup:
                g=recordgroup[c['sample_id']]
                inherited[g['component_id']]=dict(component_id=g['component_id'],via_exact_inventory_sample=c['sample_id'],training_members=sorted(trainids.intersection(g['members'])))
        out=dict(sequence_id='seq_'+hashlib.sha256(seq.encode()).hexdigest(),length=len(seq),sequence=seq,
            exact_training_chains=[dict(sample_id=c['sample_id'],chain=c['chain'],system=c['biological_system']) for c in trainseq.get(seq,[])],
            independent_local_threshold_hits=hits,inherited_exact_sequence_component_memberships=list(inherited.values()),
            short_chain_rule='exact_only' if len(seq)<20 else 'full_chain_local30pct_identity_and70pct_span_both',
            fully_family_disjoint_established=False)
        cache[seq]=out;return out
    template_rows=[]
    for tid,t in templates.items():
        p=Path(t['input_pdb']);st=parse_pdb(str(p));roles=t['roles'];rev={v:k for k,v in roles.items()}
        if t.get('input_sha256'):assert sha(p)==t['input_sha256']
        for chain in st.unique_chains():
            seq=st.chain_sequence(chain)
            template_rows.append(dict(template_id=tid,chain=chain,role=rev.get(chain,'unassigned'),source_path=str(p),source_sha256=sha(p),
                scope='full_chain',**sequence_audit(seq)))
            if rev.get(chain,'').startswith('tcr_') and len(seq)>120:
                template_rows.append(dict(template_id=tid,chain=chain,role=rev[chain],source_path=str(p),source_sha256=sha(p),scope='first120_variable_domain_screen',**sequence_audit(seq[:120])))
            if rev.get(chain) in {'hla_heavy','mhc'} and len(seq)>180:
                template_rows.append(dict(template_id=tid,chain=chain,role=rev[chain],source_path=str(p),source_sha256=sha(p),scope='first180_peptide_binding_domain_screen',**sequence_audit(seq[:180])))
        print('template audited',tid,flush=True)
    # Exact candidate sequence equality is cheap for all new substitutions.
    # Direct homology screens inherit the template only for same-length small
    # substitutions; this is explicitly not proof of complete family isolation.
    candidate_rows=[]
    for c in candidates:
        assert c['template_id'] in templates
        for chain,seq in c.get('chain_sequences',{}).items():
            exact=trainseq.get(seq,[])
            candidate_rows.append(dict(candidate_id=c['candidate_id'],biological_candidate_id=c.get('biological_candidate_id',c['candidate_id']),panel_id=c.get('panel_id'),template_id=c['template_id'],chain=chain,length=len(seq),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                exact_training_chains=[dict(sample_id=x['sample_id'],chain=x['chain'],system=x['biological_system']) for x in exact]))
    direct=[r for r in template_rows if r['exact_training_chains'] or r['independent_local_threshold_hits']]
    inherited=[r for r in template_rows if any(g['training_members'] for g in r['inherited_exact_sequence_component_memberships'])]
    out=dict(status='complete_local_overlap_audit',reviewed_at_utc=datetime.now(timezone.utc).isoformat(),
        manifest_path=args.manifest,manifest_sha256=sha(args.manifest),training_manifest_path=str(trainpath),training_manifest_sha256=sha(trainpath),training_records=len(records),sources=dict(Counter(r['dataset'] for r in records)),
        n_template_chains_and_domains=len(template_rows),n_direct_overlap_queries=len(direct),n_known_transitive_overlap_queries=len(inherited),n_candidate_chain_exact_hits=sum(bool(r['exact_training_chains']) for r in candidate_rows),
        templates=template_rows,candidate_changed_chains=candidate_rows,
        interpretation='Direct or known transitive overlap does not erase independence of new functional labels, but makes an unseen-protein claim invalid. Report assay-transfer within known protein families separately. Absence of the stated local screen is not proof of zero remote homology.',
        source_hashes={str(p):sha(p) for p in [REF/'sequence_audit/chain_inventory.json',BASE/'family_folds/family_components.json',Path(__file__)]},
        inference_run=False,remote_access=False,model_training_run=False)
    Path(args.output).write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps({k:out[k] for k in ['status','n_template_chains_and_domains','n_direct_overlap_queries','n_known_transitive_overlap_queries','n_candidate_chain_exact_hits']},ensure_ascii=False))

if __name__=='__main__':main()
