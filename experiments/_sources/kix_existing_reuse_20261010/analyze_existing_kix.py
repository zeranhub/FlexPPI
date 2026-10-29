"""Read-only CPU reuse of historical KIX rep01; never runs MD or model inference."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import struct
import numpy as np
import pandas as pd

PROJECT = Path('E:/Project/FlexPPI')
ROOT = PROJECT/'artifacts/kix_md_20260914'
OLD = ROOT/'pilot_results_20260916'
QUALITY = OLD/'quality'
FIG = PROJECT/'experiments/_case_validation_sources'
OUT = PROJECT/'experiments/_sources/kix_existing_reuse_20261010'
STATES = ['apo','myb','mll','both']


def require(condition,message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4*1024*1024),b''):
            digest.update(block)
    return digest.hexdigest()


def corr(left,right):
    return float(np.corrcoef(pd.Series(left).rank(method='average'),pd.Series(right).rank(method='average'))[0,1])


def dcd_header(path):
    with path.open('rb') as handle:
        opening = handle.read(4)
        endian = '<' if struct.unpack('<i',opening)[0] == 84 else '>'
        require(struct.unpack(endian+'i',opening)[0] == 84,'Unexpected DCD opening record')
        header = handle.read(84)
        require(header[:4] == b'CORD','Not a DCD coordinate trajectory')
        nframes = struct.unpack(endian+'i',header[4:8])[0]
        require(struct.unpack(endian+'i',handle.read(4))[0] == 84,'DCD header record corrupted')
        title_bytes = struct.unpack(endian+'i',handle.read(4))[0]
        handle.read(title_bytes)
        require(struct.unpack(endian+'i',handle.read(4))[0] == title_bytes,'DCD title record corrupted')
        require(struct.unpack(endian+'i',handle.read(4))[0] == 4,'DCD atom count record corrupted')
        atoms = struct.unpack(endian+'i',handle.read(4))[0]
    return nframes,atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify-trajectory-hashes',action='store_true')
    args = parser.parse_args()
    summary_dir = OUT/'summary'
    summary_dir.mkdir(parents=True,exist_ok=True)
    sources,raw_rows = {},[]
    qa = json.loads((FIG/'kix/kix_v13_refresh_QA.json').read_text(encoding='utf-8'))
    require(qa['status'] == 'pass' and qa['A00AQ_in_fold1_test'] and not qa['A00AQ_in_fold1_train'],'Wrong retained case checkpoint protocol')
    require(qa['checkpoint_sha256'] == '63b48e271413a91173c5763a1607feb365daf1bc48818c237588718903355ebb','Wrong fixed fold1 checkpoint')
    profile_path = FIG/'data/kix_fold1_all_residues.csv'
    metric_path = FIG/'data/kix_fold1_case_metrics.csv'
    for path in [profile_path,metric_path]:
        require(sha(path) == qa['source_tables'][path.name],'Corrected figure source was altered: '+str(path))
        sources[str(path)] = sha(path)
    profile = pd.read_csv(profile_path)
    metrics = pd.read_csv(metric_path).set_index('state')
    require(len(profile) == 348 and set(profile.metric) == {'rmsf'} and set(profile.model_fold) == {1},'Expected four corrected RMSF-only profiles')
    qsummary = json.loads((QUALITY/'summary.json').read_text(encoding='utf-8'))
    for name,entry in qsummary['outputs'].items():
        path = QUALITY/name
        require(sha(path) == entry['sha256'],'Historical quality source hash mismatch: '+name)
        sources[str(path)] = sha(path)
    sources[str(QUALITY/'summary.json')] = sha(QUALITY/'summary.json')
    sources[str(FIG/'kix/kix_v13_refresh_QA.json')] = sha(FIG/'kix/kix_v13_refresh_QA.json')
    window = pd.read_csv(QUALITY/'residue_window_metrics.csv')
    trace = pd.read_csv(QUALITY/'trajectory_trace_binned.csv')
    thermo = pd.read_csv(QUALITY/'thermodynamics_all_samples.csv')
    state_rows,quality_rows,stability_rows = [],[],[]
    for state in STATES:
        previous_path = OLD/state/'summary.json'
        previous = json.loads(previous_path.read_text(encoding='utf-8'))
        detail = next(row for row in qsummary['states'] if row['state'] == state)
        md_path = OLD/state/'residue_metrics.csv'
        require(sha(previous_path) == detail['prior_summary_sha256'],'Historical analyzer summary changed')
        require(sha(md_path) == detail['original_residue_metrics_sha256'],'Historical MD label source changed')
        expected = next(row for row in qa['sources'] if row['state'] == state)
        require(sha(md_path) == expected['md_sha256'],'Current figure MD reference differs from quality reference')
        require(detail['dcd_hashes_verified'] and previous['n_frames'] == 100000,'Incomplete historical quality analysis')
        require(detail['full_repeat_metric_max_abs_difference'] == {'rmsf_angstrom':0.0,'chi_entropy_nats':0.0},'Windows do not reproduce original full labels')
        protein = PROJECT/'artifacts/external/kix_md_20260914'/state/'protein.pdb'
        require(sha(protein) == previous['protein_pdb_sha256'],'Original protein topology changed')
        production = protein.parent/'production'
        paths = sorted(production.glob('segment_*.dcd'))
        require(len(paths) == 100 and len(previous['segments']) == 100,'Incomplete production segments')
        for path,entry in zip(paths,previous['segments']):
            require(path.name == entry['filename'] and path.stat().st_size == entry['bytes'],'Original production file identity mismatch')
            nframes,atoms = dcd_header(path)
            require(nframes == entry['num_frames'] == 1000 and atoms == previous['protein_atom_count'],'DCD frame/atom count mismatch')
            actual = sha(path) if args.verify_trajectory_hashes else None
            require(actual is None or actual == entry['sha256'],'Original DCD SHA256 mismatch: '+str(path))
            raw_rows.append(dict(state=state,path=str(path),bytes=path.stat().st_size,n_frames=nframes,protein_atoms=atoms,
                expected_sha256=entry['sha256'],actual_sha256=actual,hash_rechecked_now=bool(actual)))
        sources[str(previous_path)],sources[str(md_path)],sources[str(protein)] = sha(previous_path),sha(md_path),sha(protein)
        md = pd.read_csv(md_path)
        md = md[md.chain.eq('B')].set_index('author_residue_number').sort_index()
        current = profile[profile.state.eq(state)].set_index('author_residue_number').sort_index()
        require(list(current.index) == list(md.index) == list(range(586,673)),'Incorrect KIX residue identity/order')
        require(list(current.amino_acid) == list(md.amino_acid),'Incorrect KIX residue chemistry')
        require(np.allclose(current.md_value,md.rmsf_angstrom,atol=1e-12,rtol=0),'Corrected figure reference differs from original MD values')
        rho = corr(current.md_value,current.model_value)
        mae = float(np.mean(np.abs(current.md_value-current.model_value)))
        require(abs(rho-metrics.loc[state,'spearman']) < 1e-12,'Figure rho does not reproduce')
        require(abs(mae-metrics.loc[state,'mae']) < 1e-6,'Figure MAE does not reproduce within CSV precision')
        loop = current.loc[614:621]
        state_rows.append(dict(state=state,independent_trajectories=1,equilibration_ns=10,production_ns=100,frames=100000,
            KIX_residues=87,RMSF_MAE_angstrom=float(metrics.loc[state,'mae']),RMSF_spearman=float(metrics.loc[state,'spearman']),
            top20_overlap=int(metrics.loc[state,'top_overlap']),top20_nominated=18,
            top20_hit_fraction=float(metrics.loc[state,'top_overlap']/18),
            KIX_MD_mean_RMSF_angstrom=float(current.md_value.mean()),KIX_fold1_mean_RMSF_angstrom=float(current.model_value.mean()),
            L12_residues=8,L12_MD_mean_RMSF_angstrom=float(loop.md_value.mean()),L12_fold1_mean_RMSF_angstrom=float(loop.model_value.mean()),
            L12_MD_mean_percentile=float(loop.md_percentile.mean()),L12_fold1_mean_percentile=float(loop.model_percentile.mean()),
            L12_Hchi_available_including_Pro=int(md.loc[614:621].hchi_available.sum()),
            L12_MD_mean_Hchi_including_Pro_nats=float(md.loc[614:621].chi_entropy_nats.mean())))
        for metric in ['KIX_CA_RMSD_angstrom','L12_CA_RMSD_after_KIX_fit_angstrom']:
            values = trace[trace.state.eq(state)&trace.metric.eq(metric)]
            require(len(values) == 1000 and set(values.n_frames) == {100},'Incorrect all-frame RMSD bins')
            first = values[values.end_ns_inclusive.le(50)]
            second = values[values.start_ns_exclusive.ge(50)]
            quality_rows.append(dict(state=state,metric=metric,mean=float(values['mean'].mean()),
                minimum=float(values.minimum.min()),maximum=float(values.maximum.max()),
                first_half_mean=float(first['mean'].mean()),second_half_mean=float(second['mean'].mean()),
                first_10ns_mean=float(values[values.end_ns_inclusive.le(10)]['mean'].mean()),
                final_10ns_mean=float(values[values.start_ns_exclusive.ge(90)]['mean'].mean()),
                interpretation='All-frame descriptive displacement from production frame1 after KIX CA fit; not convergence proof'))
        for partner,entry in detail['contacts'].items():
            values = trace[trace.state.eq(state)&trace.metric.eq('partner_'+partner+'_native_residue_contact_fraction')]
            require(len(values) == 1000,'Missing native-contact trace')
            quality_rows.append(dict(state=state,metric='partner_'+partner+'_native_residue_contact_fraction',
                mean=entry['mean_residue_pair_retention'],minimum=float(values.minimum.min()),maximum=float(values.maximum.max()),
                first_half_mean=float(values[values.end_ns_inclusive.le(50)]['mean'].mean()),
                second_half_mean=float(values[values.start_ns_exclusive.ge(50)]['mean'].mean()),
                first_10ns_mean=float(values[values.end_ns_inclusive.le(10)]['mean'].mean()),
                final_10ns_mean=float(values[values.start_ns_exclusive.ge(90)]['mean'].mean()),
                native_residue_pairs=entry['native_residue_pairs'],interpretation='Retention of fixed initial <=4.5A contacts; excludes new contacts'))
        temperature = thermo[thermo.state.eq(state)]
        require(len(temperature) == 10000,'Missing reporter samples')
        quality_rows.append(dict(state=state,metric='temperature_kelvin',mean=float(temperature.temperature_kelvin.mean()),
            minimum=float(temperature.temperature_kelvin.min()),maximum=float(temperature.temperature_kelvin.max()),
            first_half_mean=float(temperature[temperature.production_ns.le(50)].temperature_kelvin.mean()),
            second_half_mean=float(temperature[temperature.production_ns.gt(50)].temperature_kelvin.mean())))
        for start in [0,50]:
            half = window[window.state.eq(state)&window.chain.eq('B')&window.window_kind.eq('disjoint_half')&window.start_ns_exclusive.eq(start)]
            require(len(half) == 87,'Missing disjoint-half RMSF profile')
            half = half.set_index('author_residue_number').sort_index()
            stability_rows.append(dict(state=state,start_ns_exclusive=start,end_ns_inclusive=start+50,
                KIX_mean_RMSF_angstrom=float(half.rmsf_angstrom.mean()),L12_mean_RMSF_angstrom=float(half.loc[614:621].rmsf_angstrom.mean()),
                rho_RMSF_vs_full=float(corr(half.rmsf_angstrom,md.rmsf_angstrom)),
                L12_mean_Hchi_including_Pro_nats=float(half.loc[614:621].chi_entropy_nats.mean())))
    state_table = pd.DataFrame(state_rows).set_index('state')
    deltas = []
    for left,right in itertools.combinations(STATES,2):
        md_delta = state_table.loc[right,'L12_MD_mean_RMSF_angstrom']-state_table.loc[left,'L12_MD_mean_RMSF_angstrom']
        model_delta = state_table.loc[right,'L12_fold1_mean_RMSF_angstrom']-state_table.loc[left,'L12_fold1_mean_RMSF_angstrom']
        deltas.append(dict(contrast=right+'-'+left,positive_definition=right+' minus '+left,L12_residues=8,
            MD_delta_RMSF_angstrom=float(md_delta),fold1_delta_RMSF_angstrom=float(model_delta),
            absolute_delta_error_angstrom=float(abs(md_delta-model_delta)),
            direction_agrees=bool(np.sign(md_delta)==np.sign(model_delta))))
    outputs = {'state_summary.csv':pd.DataFrame(state_rows),'l12_delta_summary.csv':pd.DataFrame(deltas),
        'quality_summary.csv':pd.DataFrame(quality_rows),'half_window_stability.csv':pd.DataFrame(stability_rows),
        'original_trajectory_inventory.csv':pd.DataFrame(raw_rows)}
    for name,data in outputs.items():
        data.to_csv(summary_dir/name,index=False)
    audit = dict(status='complete',analysis='Existing historical single trajectories only; no MD or model inference performed',
        original_trajectories=4,original_production_ns=400,original_equilibration_ns=40,
        raw_DCD_files=400,raw_frames=400000,raw_file_hashes_rechecked_now=args.verify_trajectory_hashes,
        corrected_fold1_checkpoint_sha256=qa['checkpoint_sha256'],model_source=str(profile_path),
        figure_CSV_max_expected_numeric_precision_difference=1e-6,source_hashes=sources,
        output_hashes={name:sha(summary_dir/name) for name in outputs},script_sha256=sha(__file__),
        limitations=['One independent trajectory per state; no independent-repeat CI, p value or statistical-significance claim.',
            'RMSF windows are descriptive and do not establish convergence; cumulative windows overlap.',
            'Native contacts refer to fixed initial heavy-atom pairs; low retention can also reflect interface rearrangement.',
            'Current corrected figure CSV is RMSF only; old74-site Hchi or old five-model deltas are not used.',
            'Fold1 record A00AQ was held out; this case is not an all-family exclusion evaluation.'])
    (OUT/'reuse_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(pd.DataFrame(state_rows)[['state','RMSF_MAE_angstrom','RMSF_spearman','top20_overlap','L12_MD_mean_RMSF_angstrom','L12_fold1_mean_RMSF_angstrom']].to_string(index=False))
    print(pd.DataFrame(deltas).to_string(index=False))
    print(pd.DataFrame(quality_rows)[['state','metric','mean','first_half_mean','second_half_mean']].to_string(index=False))
    print(pd.DataFrame(stability_rows).to_string(index=False))
    print('Audit',OUT/'reuse_audit.json')


if __name__ == '__main__':
    main()
