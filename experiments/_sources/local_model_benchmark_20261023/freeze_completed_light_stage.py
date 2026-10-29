"""Freeze a reviewed 14-entry stage. Default is a read-only prerequisite check.

Use --freeze only after the actual 14-entry PNGs have been visually inspected
and qa_current.py has bound that inspection to the current image hashes.
Standard library only: arrays are never deserialized and no models are loaded.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import shutil
import sys
import uuid

ROOT = Path(__file__).resolve().parent
DESTINATION = ROOT/'snapshots/completed_14_entries'
BASE_SCOPE = ROOT/'snapshots/completed_11_entries_exposure_audited/SNAPSHOT_SCOPE.json'
EXPECTED = {
    'Full','Intra','Static','AA','DynaProt-M','ANM','GNM','FastProtFlex',
    'backflip-2.1 (assembly)','backflip-2.1 (chains)',
    'backflip-1.0 (assembly)','backflip-1.0 (chains)',
    'backflip-1.0-seq (assembly)','backflip-1.0-seq (chains)',
}
GROUPS = {'single_chain':2193,'other_multichain':626,'tcr_pmhc_itrace':244}
QA_STATUS = 'PASS_CURRENT_SNAPSHOT_WITH_DOCUMENTED_VALIDATOR_LIMITATIONS'
MEMBERSHIP_SHA = '72521f5abc8609ce7bd4dcf10cf3e1e6c0d16626017f06e9a7be98fdf6942e5b'
EXPOSURE_SHA = '83b0fbbddd486633cf2bb2c3dd34de8da4d896d6fd05378ff2f9498836454575'
RESOURCE_POLICY = 'lifetime_v2_rss115_reserve256_commit1024_v2'
SELECTED_NUMERICAL_AUDIT = 'fastprotflex/memory_lifetime_v2/runs/smoke_01/validation_audit.json'
SELECTED_NUMERICAL_SHA = 'a283aceaaf2eb4f5e1027c77418e61903de5c85136a157bb3c48849ebee1df81'
SELECTED_EXTREME_AUDIT = 'fastprotflex/memory_lifetime_v2/resource_probes/extreme_03/resource_extension_audit.json'
SELECTED_EXTREME_SHA = '2e4b2ba36c6b5bbae876aaadce22496381d1b4bf2fbd13225d037b9d821ea146'
SHA_OBSERVATIONS = {}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def remember_sha(path, digest):
    key = str(Path(path).resolve()).casefold()
    previous = SHA_OBSERVATIONS.setdefault(key,digest)
    require(previous==digest, f'File changed after its first verified read: {path}')
    return digest


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024**2), b''):
            digest.update(block)
    return remember_sha(path,digest.hexdigest())


def load_json(path):
    require(Path(path).is_file(), f'Required audit missing: {path}')
    raw = Path(path).read_bytes()
    remember_sha(path,hashlib.sha256(raw).hexdigest())
    return json.loads(raw.decode('utf-8-sig'))


def utc(value):
    result = datetime.fromisoformat(value.replace('Z','+00:00'))
    require(result.tzinfo is not None, f'UTC timestamp lacks timezone: {value}')
    return result.astimezone(timezone.utc)


def exact_methods(values, label):
    require(len(values)==14 and set(values)==EXPECTED,
            f'{label} must contain exactly the 14 completed entries: {values}')


def relative_path(value):
    value = str(value).replace('\\','/')
    path = PurePosixPath(value)
    require(not path.is_absolute() and '..' not in path.parts and ':' not in value,
            f'Unsafe snapshot-relative path: {value}')
    require(path.parts and str(path) not in ('.',''), f'Empty path: {value}')
    return path.as_posix()


def checked_source(relative):
    source = ROOT/relative_path(relative)
    require(source.is_file(), f'Required source file missing: {source}')
    require(source.resolve().is_relative_to(ROOT.resolve()), f'Source escapes task root: {source}')
    return source


def resource_requirements(calibration, n_atoms):
    require(calibration['resource_policy_version']==RESOURCE_POLICY and
            calibration['safety_multiplier']==1.15 and
            calibration['additional_available_RAM_reserve_MiB']==256.0 and
            calibration['additional_available_commit_reserve_MiB']==1024.0 and
            calibration['prior_additional_available_RAM_reserve_MiB']==512.0 and
            calibration.get('policy_change_reason'), 'Unknown/unapproved v2 resource policy')
    estimate = max(calibration['minimum_peak_floor_MiB'], calibration['quadratic_intercept_MiB']+
                   calibration['quadratic_slope_MiB_per_atom_squared']*n_atoms**2)
    require(math.isfinite(estimate) and estimate>0, 'Invalid predicted resource peak')
    return dict(estimated_R_process_tree_RSS_MiB=estimate,
        required_available_RAM_MiB=1.15*estimate+256.0,
        previous_required_available_RAM_MiB=1.15*estimate+512.0,
        required_available_commit_MiB=estimate+1024.0,safety_multiplier=1.15,
        additional_available_RAM_reserve_MiB=256.0,additional_available_commit_reserve_MiB=1024.0,
        resource_policy_version=RESOURCE_POLICY,policy_change_reason=calibration['policy_change_reason'])


def verify_launch_gate(gate, calibration, n_atoms):
    for key,value in resource_requirements(calibration,n_atoms).items():
        actual = gate[key]
        if isinstance(value,(int,float)):
            require(math.isfinite(float(actual)) and abs(actual-value)<1e-8, f'Actual resource threshold differs: {key}')
        else:
            require(actual==value, f'Actual resource policy differs: {key}')
    available = gate.get('available_RAM_before_launch_MiB',gate.get('available_RAM_before_launch_mib'))
    commit = gate['available_system_commit_MiB']
    require(available is not None and math.isfinite(float(available)) and
            available>=gate['required_available_RAM_MiB'] and
            commit is not None and math.isfinite(float(commit)) and commit>=gate['required_available_commit_MiB'] and
            gate['system_commit_counter']=='Windows_GetPerformanceInfo' and
            gate['system_commit_counter_error'] is None and
            math.isfinite(float(gate['system_commit_limit_MiB'])) and gate['system_commit_limit_MiB']>=commit and
            not gate['blocking_jobs'], 'Actual physical/Windows commit launch gate did not pass')


def verify_optional_resource_measurements(audit):
    """Unknown sampling stays unknown; successful launch still needs known commit."""
    for status_key,known_key,unknown_key,value_key in (
        ('private_bytes_counter_status','known_private_samples','unknown_private_samples','peak_R_process_tree_private_bytes'),
        ('system_commit_sampling_status','known_commit_samples','unknown_commit_samples','minimum_system_available_commit_MiB')):
        known,unknown = audit[known_key],audit[unknown_key]
        require(isinstance(known,int) and isinstance(unknown,int) and known>=0 and unknown>=0,
                'Resource sampling counts are invalid')
        expected = 'unknown' if not known else 'available' if not unknown else 'partially_available'
        require(audit[status_key]==expected, f'Resource sampling status/counts differ: {status_key}')
        measured = audit[value_key]
        require((not known and measured is None) or
                (known>0 and measured is not None and math.isfinite(float(measured)) and measured>=0),
                f'Unknown/known resource samples are misrepresented: {value_key}')
    minimum = audit['minimum_system_available_physical_MiB']
    require(minimum is None or (math.isfinite(float(minimum)) and minimum>=0), 'Invalid sampled physical RAM minimum')


def verify_fast_protocol(protocol, full_audit, add, requested_audits):
    """Bind only the actually adopted entry point and its numeric validation."""
    def bound(value, digest, label):
        path = Path(value)
        source = path if path.is_absolute() else ROOT/path
        require(source.is_file() and source.resolve().is_relative_to(ROOT.resolve()),
                f'{label}: required bound file is missing/outside task root: {value}')
        require(sha(source)==digest, f'{label}: file SHA differs from actual protocol: {value}')
        add(source.resolve().relative_to(ROOT.resolve()).as_posix())
        return source.resolve()

    require(protocol['regression_coefficients_unchanged'] is True and
            protocol['no_model_fitting'] is True and protocol['complete_structure_input'] is True,
            'FastProtFlex scientific protocol changed')
    require(protocol['full_schedule_sha256']==sha(ROOT/'fastprotflex/full_schedule.csv') and
            protocol['input_cohort_manifest_sha256']==sha(ROOT/'cohort/records.csv'),
            'FastProtFlex protocol does not bind the frozen cohort/schedule')
    version = protocol.get('algorithm_lifecycle_version','official_v1')
    require(version in ('official_v1','lifetime_v2'), f'Unknown FastProtFlex lifecycle version: {version}')
    for name, digest in protocol['upstream_sha256'].items():
        require(name in ('predict.r','FUNCTION_GDV.r','LICENSE','README.md'), 'Unexpected upstream source name')
        bound('fastprotflex/upstream/'+name,digest,'Official source mirror')
    if version=='official_v1':
        require(not requested_audits, 'Cannot attach inactive v2 candidates to an official_v1 snapshot')
        bound('fastprotflex/run_local.py',protocol['wrapper_sha256'],'Official v1 entry point')
        bound('fastprotflex/capture_full_precision.R',protocol['capture_wrapper_sha256'],'Capture wrapper')
        require(all(r.get('algorithm_lifecycle_version','official_v1')=='official_v1'
                    for r in full_audit['results']), 'v2 record exists without an adopted v2 protocol')
        return dict(algorithm_lifecycle_version=version,actual_validation_audits=[],external_sources=[])

    entry = bound(protocol['execution_entrypoint_path'],protocol['execution_entrypoint_sha256'],'Actual v2 entry point')
    require(protocol['resource_policy_version']==RESOURCE_POLICY and protocol.get('policy_change_reason'),
            'Adopted v2 protocol does not identify the authorized physical/commit resource policy')
    require(bound(protocol['wrapper_path'],protocol['wrapper_sha256'],'Wrapper alias')==entry and
            protocol['resume_entrypoint_sha256']==protocol['execution_entrypoint_sha256'],
            'v2 entry point aliases disagree')
    bound(protocol['capture_wrapper_path'],protocol['capture_wrapper_sha256'],'Capture wrapper')
    candidate = bound(protocol['candidate_function_path'],protocol['candidate_function_sha256'],'Adopted candidate function')
    official = Path(protocol['official_algorithm_source_path'])
    require(official.is_file() and sha(official)==protocol['official_algorithm_source_sha256'],
            'Actual official algorithm source SHA mismatch')
    bound('fastprotflex/upstream/FUNCTION_GDV.r',protocol['official_algorithm_source_sha256'],'Official algorithm mirror')
    archive = bound(protocol['official_v1_protocol_archive_path'],protocol['official_v1_protocol_archive_sha256'],'Archived v1 protocol')
    original = load_json(archive)
    require(original.get('algorithm_lifecycle_version','official_v1')=='official_v1', 'Archive is not the original v1 protocol')
    bound('fastprotflex/run_local.py',original['wrapper_sha256'],'Preserved v1 entry point')
    require(original['capture_wrapper_sha256']==protocol['capture_wrapper_sha256'] and
            original['upstream_sha256']==protocol['upstream_sha256'], 'v2 changed original scientific source/capture hashes')
    audit_path = bound(protocol['lifecycle_validation_audit_path'],protocol['lifecycle_validation_audit_sha256'],
                       'Actual three-source audit')
    selected = audit_path.relative_to(ROOT.resolve()).as_posix()
    require(selected==SELECTED_NUMERICAL_AUDIT and sha(audit_path)==SELECTED_NUMERICAL_SHA,
            'v2 numeric evidence is not the selected actual smoke_01 audit')
    require(not requested_audits or {relative_path(p) for p in requested_audits}=={selected},
            'Supplied v2 audit must be exactly the protocol-bound actual audit')
    audit = load_json(audit_path)
    require(audit['status']=='three_sources_exact_all_atom_match' and
            audit['cases_requested']==audit['cases_attempted']==3,
            'Actual v2 numeric validation is incomplete; static preparation is insufficient')
    provenance_path = ROOT/'fastprotflex/memory_lifetime_v2/lifecycle_provenance.json'
    provenance = load_json(provenance_path)
    require(audit['provenance']==provenance and provenance['non_cleanup_program_text_identical_in_order'] is True and
            provenance['all_diff_executable_lines_are_only_rm_or_gc'] is True,
            'v2 source derivation differs from the numerically validated provenance')
    require(provenance['candidate_function_sha256']==protocol['candidate_function_sha256'] and
            provenance['official_function_sha256']==protocol['official_algorithm_source_sha256'] and
            provenance['capture_wrapper_sha256']==protocol['capture_wrapper_sha256'] and
            provenance['official_predict_sha256']==protocol['upstream_sha256']['predict.r'],
            'v2 protocol/provenance scientific hashes disagree')
    for relative in ('lifecycle_provenance.json','prepare_static.py','README.md'):
        add('fastprotflex/memory_lifetime_v2/'+relative)
    bound('fastprotflex/memory_lifetime_v2/validate_candidate.py',audit['validation_script_sha256'],'Actual validator source')
    bound('fastprotflex/memory_lifetime_v2/capture_and_validate.R',audit['R_validation_wrapper_sha256'],'Actual R validator')
    bound('fastprotflex/memory_lifetime_v2/FUNCTION_GDV_lifecycle_only.diff',provenance['diff_sha256'],'Lifecycle-only diff')
    bound('fastprotflex/memory_lifetime_v2/original/FUNCTION_GDV.r',provenance['official_function_sha256'],'Derivation original')
    bound('fastprotflex/memory_lifetime_v2/wrapper/predict.r',provenance['official_predict_sha256'],'Unchanged official predict wrapper')
    bound('fastprotflex/memory_lifetime_v2/wrapper/capture_full_precision.R',provenance['capture_wrapper_sha256'],'Unchanged capture wrapper')
    require(candidate==checked_source('fastprotflex/memory_lifetime_v2/candidate/FUNCTION_GDV.r').resolve(),
            'Adopted candidate path differs from the audited derivation')
    cases = {c['sample_id']:c for c in provenance['cases']}
    require(set(cases)=={'2o70_F','A007V','itrace_1qrn'} and len(audit['results'])==3 and
            {r['sample_id'] for r in audit['results']}==set(cases), 'Actual three-source cases differ from the predefined cases')
    checks = ('passed','same_atom_count','same_columns','all_atom_identities_identical',
              'all_atom_coordinates_identical','all_atom_table_identical','all_predictions_finite',
              'all_atom_prediction_vector_identical_num_eq_FALSE')
    points = []
    for result in audit['results']:
        case = cases[result['sample_id']]
        comparison = result['all_atom_RDS_comparison']
        require(result['status']=='complete_exact_all_atom_match' and result['return_code']==0 and
                result['timed_out'] is False and result['R_and_package_versions_match_baseline'] is True and
                result['library_thread_caps']==1 and 0<float(result['sampling_interval_seconds'])<=.1 and
                all(comparison[k] is True for k in checks) and
                comparison['maximum_absolute_prediction_difference']==0 and
                comparison['reference_atom_count']==comparison['candidate_atom_count']==case['n_atoms'],
                f'Actual all-atom equality failed: {case["sample_id"]}')
        require(result['candidate_function_sha256']==protocol['candidate_function_sha256'] and
                result['input_sha256']==case['input_sha256'] and result['reference_RDS_sha256']==case['baseline_RDS_sha256'],
                f'Validation input/reference/source SHA differs: {case["sample_id"]}')
        directory = audit_path.parent/case['sample_id']
        require(load_json(directory/'validation.json')==result and load_json(directory/'comparison.json')==comparison,
                f'Actual case JSON no longer matches audit: {case["sample_id"]}')
        for name in ('validation.json','comparison.json','stdout.log','stderr.log','R_session_info.txt'):
            add((directory/name).relative_to(ROOT.resolve()).as_posix())
        for name,digest in (
                ('atom_predictions_exact.rds',result['candidate_RDS_sha256']),
                ('FUNCTION_GDV.r',protocol['candidate_function_sha256']),('predict.r',provenance['official_predict_sha256']),
                ('capture_full_precision.R',protocol['capture_wrapper_sha256']),
                ('capture_and_validate.R',audit['R_validation_wrapper_sha256']),('protein.pdb',case['input_sha256'])):
            bound(directory/name,digest,'Actual validated case artifact')
        baseline = bound(case['baseline_RDS'],case['baseline_RDS_sha256'],'Exact baseline RDS')
        bound(case['baseline_validation'],case['baseline_validation_sha256'],'Original baseline validation')
        add((baseline.parent/'R_session_info.txt').relative_to(ROOT).as_posix())
        peak = float(result['peak_R_process_tree_RSS_bytes'])/1024**2
        require(math.isfinite(peak) and peak>0 and abs(peak-float(result['peak_R_process_tree_RSS_MiB']))<1e-8 and
                float(result['observed_process_tree_CPU_seconds'])>0, 'Actual RAM/CPU measurements missing')
        points.append(dict(sample_id=case['sample_id'],n_atoms=case['n_atoms'],measured_peak_MiB=peak))
    calibration = protocol['resource_calibration']
    def resource_fit(value, expected_points):
        canonical = {k:v for k,v in value.items() if k!='calibration_sha256'}
        require(value['calibration_sha256']==hashlib.sha256(json.dumps(canonical,sort_keys=True).encode()).hexdigest(),
                'Resource calibration self-hash mismatch')
        require(value['measured_points']==expected_points and
                value['validation_audit_sha256']==protocol['lifecycle_validation_audit_sha256'] and
                Path(value['validation_audit_path']).resolve()==audit_path and
                value['maximum_measured_atoms']==max(p['n_atoms'] for p in expected_points) and
                value['minimum_measured_atoms']==min(p['n_atoms'] for p in expected_points) and
                value['safety_multiplier']==1.15 and value['additional_available_RAM_reserve_MiB']==256.0 and
                value['additional_available_commit_reserve_MiB']==1024.0 and
                value['prior_additional_available_RAM_reserve_MiB']==512.0 and
                value['resource_policy_version']==RESOURCE_POLICY and
                value['policy_change_reason']==protocol['policy_change_reason'],
                'Resource calibration is not bound to the actual measurements')
        slope = float(value['quadratic_slope_MiB_per_atom_squared'])
        intercept = float(value['quadratic_intercept_MiB'])
        require(math.isfinite(slope) and math.isfinite(intercept) and slope>=0 and intercept>=0 and
                all(intercept+slope*p['n_atoms']**2+1e-8>=p['measured_peak_MiB'] for p in expected_points) and
                value['minimum_peak_floor_MiB']==min(p['measured_peak_MiB'] for p in expected_points),
                'Resource calibration no longer upper-bounds the actual measured peaks')
        require(value['policy_change_thresholds_at8058']==resource_requirements(value,8058),
                'Resource policy change thresholds are not bound to this calibration')

    external_sources = [dict(path=str(official),sha256=protocol['official_algorithm_source_sha256'])]
    extension_evidence = None
    extension_keys = ('resource_extension_audit_path','resource_extension_audit_sha256',
                      'resource_extension_probe_path','resource_extension_probe_sha256')
    require(all(key in protocol for key in extension_keys),
            'This completed v2 stage requires its selected successful extreme_03 resource audit')
    if any(key in protocol or key in calibration for key in extension_keys):
        require(all(key in protocol and key in calibration and protocol[key]==calibration[key]
                    for key in extension_keys), 'Incomplete resource-extension protocol/calibration binding')
        extension_path = bound(protocol['resource_extension_audit_path'],protocol['resource_extension_audit_sha256'],
                               'Actual extreme resource audit')
        require(extension_path==checked_source(SELECTED_EXTREME_AUDIT).resolve() and
                sha(extension_path)==SELECTED_EXTREME_SHA,
                'Adopted resource evidence must be the selected successful extreme_03 audit')
        probe = bound(protocol['resource_extension_probe_path'],protocol['resource_extension_probe_sha256'],
                      'Actual extreme resource probe source')
        extension = load_json(extension_path)
        require(extension['status']=='extreme_resource_and_integrity_check_complete_no_reference_equivalence' and
                extension['sample_id']=='7zh9_A' and extension['n_atoms']==8058 and
                extension['original_full_precision_reference_available'] is False and
                extension['original_equivalence_compared'] is False and extension['resource_policy_version']==RESOURCE_POLICY,
                'Extreme probe is incomplete or incorrectly claims original-output equivalence')
        require(Path(extension['numerical_equivalence_audit_path']).resolve()==audit_path and
                extension['numerical_equivalence_audit_sha256']==protocol['lifecycle_validation_audit_sha256'] and
                extension['candidate_function_sha256']==protocol['candidate_function_sha256'] and
                extension['official_function_sha256']==protocol['official_algorithm_source_sha256'] and
                extension['official_predict_sha256']==provenance['official_predict_sha256'] and
                extension['capture_wrapper_sha256']==protocol['capture_wrapper_sha256'] and
                Path(extension['probe_script_path']).resolve()==probe and
                extension['probe_script_sha256']==protocol['resource_extension_probe_sha256'] and
                Path(extension['resume_entrypoint_path']).resolve()==entry and
                extension['resume_entrypoint_sha256']==protocol['execution_entrypoint_sha256'] and
                extension['cohort_manifest_sha256']==protocol['input_cohort_manifest_sha256'] and
                extension['schedule_sha256']==protocol['full_schedule_sha256'],
                'Extreme audit source/cohort/three-source bindings differ')
        base_calibration = calibration['base_three_source_resource_calibration']
        require(extension['base_three_source_resource_calibration']==base_calibration,
                'Extreme gate no longer preserves the actual original three-source calibration')
        resource_fit(base_calibration,points)
        gate = extension['gate']
        verify_launch_gate(gate,base_calibration,8058)
        verify_optional_resource_measurements(extension)
        require(gate['status']=='ready' and gate['source_calibration_sha256']==base_calibration['calibration_sha256'] and
                extension['return_code']==0 and extension['timed_out'] is False and
                extension['R_and_package_versions_match_baseline'] is True and extension['library_thread_caps']==1 and
                0<float(extension['sampling_interval_seconds'])<=.1 and extension['observed_process_tree_CPU_seconds']>0,
                'Extreme actual launch/resource/version proof is incomplete')
        integrity = extension['integrity']
        require(integrity['n_atoms']==8058 and integrity['n_CA']==1024 and
                all(integrity[key] is True for key in ('all_atom_identities_equal_input_in_order',
                    'all_atom_coordinates_verified','all_atom_predictions_finite','ca_identities_equal_frozen_truth_in_order',
                    'all_ca_predictions_finite','no_duplicate_CA_identity')) and
                0<=integrity['coordinate_maximum_absolute_difference']<integrity['coordinate_verification_tolerance']==1e-9 and
                integrity['absolute_RMSF_Angstrom'] is False and integrity['Hchi_supported'] is False,
                'Extreme atom/CA integrity or relative-output proof is incomplete')
        with (ROOT/'cohort/records.csv').open(encoding='utf-8-sig',newline='') as handle:
            records = [r for r in csv.DictReader(handle) if r['sample_id']=='7zh9_A']
        require(len(records)==1 and int(records[0]['n_residues'])==1024 and
                extension['input_sha256']==records[0]['structure_sha256']==sha(records[0]['structure_path']),
                'Extreme input is not the frozen 7zh9_A structure')
        external_sources.append(dict(path=records[0]['structure_path'],sha256=extension['input_sha256']))
        output_rows = extension['output_files']
        required_outputs = {'7zh9_A/'+name for name in ('FUNCTION_GDV.r','predict.r','capture_full_precision.R',
            'protein.pdb','prediction.npz','atom_predictions_exact.rds','atom_predictions_full_precision.csv',
            'proteinRMSF.pdb','R_session_info.txt','stdout.log','stderr.log','integrity.json')}
        output_keys = {row['relative_path'] for row in output_rows}
        require(len(output_rows)==len(output_keys) and
                required_outputs<=output_keys<=required_outputs|{'waiting_for_memory.json'},
                'Extreme probe output manifest lacks/duplicates a required real artifact')
        actual_outputs = {p.relative_to(extension_path.parent).as_posix() for p in extension_path.parent.rglob('*')
                          if p.is_file() and p.resolve()!=extension_path}
        require(actual_outputs==output_keys, 'Extreme probe directory contains unbound/missing real output files')
        output_lookup = {}
        for row in output_rows:
            output = extension_path.parent/relative_path(row['relative_path'])
            require(output.resolve().is_relative_to(extension_path.parent.resolve()), 'Extreme artifact escapes its audit directory')
            bound(output,row['sha256'],'Actual extreme probe artifact')
            require(output.stat().st_size==row['size_bytes'], 'Extreme probe artifact byte count differs')
            output_lookup[row['relative_path']] = row['sha256']
        require(load_json(extension_path.parent/'7zh9_A/integrity.json')==integrity and
                output_lookup['7zh9_A/FUNCTION_GDV.r']==protocol['candidate_function_sha256'] and
                output_lookup['7zh9_A/predict.r']==provenance['official_predict_sha256'] and
                output_lookup['7zh9_A/capture_full_precision.R']==protocol['capture_wrapper_sha256'] and
                output_lookup['7zh9_A/protein.pdb']==extension['input_sha256'], 'Extreme copied sources/integrity JSON differ')
        peak = float(extension['peak_R_process_tree_RSS_bytes'])/1024**2
        require(math.isfinite(peak) and peak>0 and abs(peak-float(extension['peak_R_process_tree_RSS_MiB']))<1e-8,
                'Actual extreme RAM measurement missing')
        points.append(dict(sample_id='7zh9_A',n_atoms=8058,measured_peak_MiB=peak))
        envelope_keys = ('version','measured_points','quadratic_intercept_MiB',
                         'quadratic_slope_MiB_per_atom_squared','maximum_measured_atoms','minimum_measured_atoms',
                         'minimum_peak_floor_MiB','safety_multiplier','additional_available_RAM_reserve_MiB',
                         'additional_available_commit_reserve_MiB','prior_additional_available_RAM_reserve_MiB',
                         'resource_policy_version','policy_change_reason')
        require(extension['combined_resource_envelope']=={key:calibration[key] for key in envelope_keys},
                'Actual probe combined resource envelope differs from the adopted four-point calibration')
        extension_evidence = dict(audit_path=str(extension_path),audit_sha256=sha(extension_path),
            original_equivalence_compared=False,resource_policy_version=RESOURCE_POLICY,
            actual_peak_R_process_tree_RSS_MiB=peak,
            actual_launch_gate=gate,
            resource_sampling={key:extension[key] for key in ('peak_R_process_tree_private_bytes',
                'private_bytes_counter_status','known_private_samples','unknown_private_samples',
                'minimum_system_available_commit_MiB','system_commit_sampling_status',
                'known_commit_samples','unknown_commit_samples','minimum_system_available_physical_MiB')},
            claim='8058 atoms: actual measured resources and atom/CA integrity only; no original precision-output equivalence claim.')
        require(all(full_audit[key]==protocol[key] for key in extension_keys), 'Full audit lacks adopted resource-extension bindings')
    require(calibration['version']==protocol['resource_gate_version'], 'Protocol/calibration gate versions disagree')
    resource_fit(calibration,points)
    require(full_audit['algorithm_lifecycle_version']=='lifetime_v2' and
            full_audit['lifecycle_validation_audit_sha256']==protocol['lifecycle_validation_audit_sha256'] and
            full_audit['resource_gate_version']==protocol['resource_gate_version'] and
            full_audit['resource_policy_version']==RESOURCE_POLICY and
            full_audit['cached_frozen_truth_revalidated_this_invocation'] is False and
            full_audit['final_scorer_frozen_truth_revalidation_required'] is True,
            'Full audit does not bind adopted v2 validation or cache truth-verification scope')
    add('fastprotflex/full_audit.lifetime_v2.json')
    require(load_json(ROOT/'fastprotflex/full_audit.lifetime_v2.json')==full_audit, 'v2 and published full audit differ')
    for result in full_audit['results']:
        record_version = result.get('cache_algorithm_lifecycle_version',result.get('algorithm_lifecycle_version','official_v1'))
        require(record_version in ('official_v1','lifetime_v2'), 'Unknown cached lifecycle version')
        if result.get('cache_reused_without_forward'):
            require(result['cache_frozen_truth_revalidated_this_invocation'] is False and
                    result['final_scorer_frozen_truth_revalidation_required'] is True and
                    result.get('cache_verification_scope'), 'Reused cache truth-verification scope is overstated')
        if record_version=='lifetime_v2':
            require(result['execution_entrypoint_sha256']==protocol['execution_entrypoint_sha256'] and
                    Path(result['execution_entrypoint_path']).resolve()==entry and
                    result['candidate_function_sha256']==protocol['candidate_function_sha256'] and
                    result['official_algorithm_source_sha256']==protocol['official_algorithm_source_sha256'] and
                    result['lifecycle_validation_audit_sha256']==protocol['lifecycle_validation_audit_sha256'] and
                    result['resource_gate_calibration_sha256']==calibration['calibration_sha256'],
                    f'v2 record source/calibration binding differs: {result["sample_id"]}')
            require(result['resource_policy_version']==RESOURCE_POLICY, 'New v2 case policy version differs')
            require(int(result['input_protein_heavy_atoms'])<=calibration['maximum_measured_atoms'],
                    'v2 full inference exceeds the actually measured resource range')
            verify_launch_gate(result,calibration,int(result['input_protein_heavy_atoms']))
            if extension_evidence:
                require(all(result[key]==protocol[key] for key in extension_keys),
                        f'v2 record lacks adopted resource-extension bindings: {result["sample_id"]}')
    return dict(algorithm_lifecycle_version=version,actual_validation_audits=[selected],
                execution_entrypoint_path=str(entry),execution_entrypoint_sha256=sha(entry),
                calibration_sha256=calibration['calibration_sha256'],
                resource_policy_version=RESOURCE_POLICY,
                cache_truth_verification_scope=dict(cached_frozen_truth_revalidated_this_invocation=False,
                    final_scorer_frozen_truth_revalidation_required=True,
                    final_scoring_completion_required=True,
                    basis='The final current-source scoring audit must independently cover all3063 records and all961321 frozen-truth positions.'),
                resource_extension=extension_evidence,external_sources=external_sources)


def verify_prerequisites(memory_v2_audits):
    SHA_OBSERVATIONS.clear()
    require(not DESTINATION.exists(), f'Refusing to overwrite existing snapshot: {DESTINATION}')
    base = load_json(BASE_SCOPE)
    require(base['completed_entries']==11 and base['records']==3063, 'Unexpected base scope')
    sources = {}
    def add(relative):
        relative = relative_path(relative)
        source = checked_source(relative)
        key = relative.casefold()
        if key in sources:
            require(sources[key]['source']==source, f'Case-colliding source paths: {relative}')
        else:
            sources[key] = dict(path=relative,source=source,initial_sha256=sha(source))
    for row in base['files']:
        add(row['path'])
    base_unique_count = len(sources)
    for relative in [
        'freeze_completed_light_stage.py',
        'controller_state.json','controller_phase_results.json','controller.log',
        'fastprotflex/full_audit.json','fastprotflex/run_protocol.json',
        'fastprotflex/full_results.csv','fastprotflex/full_schedule.csv',
        'fastprotflex/run_local.py','fastprotflex/capture_full_precision.R',
        'fastprotflex/copied_target_parse_audit.json','fastprotflex/README.md',
        'backflip/run_benchmark.py','backflip/validate_outputs.py',
        'backflip/summarize_results.py','backflip/dependency_patch_audit.json',
        'backflip/input_audit.json','backflip/manifest.json','backflip/README.md',
        'backflip/run_remaining_variants.ps1','backflip/remaining_variants_queue.log',
        'backflip/all_variants_summary.log',
        'dynaprot/predictions/summary.json','dynaprot/predictions/metadata.json',
        'cohort/records.csv','cohort/chains.csv','cohort/bootstrap_membership.csv',
    ]:
        add(relative)
    for name in ('adoption_preparation_event.json','working_set_trim_events.json',
                 'resource_policy_adoption_preparation.json'):
        relative = 'fastprotflex/memory_lifetime_v2/'+name
        if (ROOT/relative).exists():
            add(relative)
    for tag in ('backflip-2.1','backflip-1.0','backflip-1.0-seq'):
        for name in ('completion_audit.json','run_metadata.json','status.jsonl'):
            add(f'backflip/results/{tag}/{name}')

    scoring = load_json(ROOT/'evaluation/scoring_audit.json')
    exact_methods(scoring['complete_methods'], 'Scoring audit')
    require(scoring['records']==3063 and scoring['residues']==961321, 'Scoring cohort mismatch')
    require(scoring.get('scoring_code_sha256')==sha(ROOT/'evaluate_benchmark.py') and
            scoring.get('script_sha256')==sha(ROOT/'evaluate_benchmark.py'),
            'Scoring output is not bound to the current scorer source')
    require(scoring['methods']['FastProtFlex']['absolute_angstrom'] is False,
            'FastProtFlex must remain a relative score, without Angstrom MAE')

    with (ROOT/'evaluation/coverage_summary.csv').open(encoding='utf-8-sig',newline='') as handle:
        coverage = [r for r in csv.DictReader(handle) if r['method'] in EXPECTED]
    require(len(coverage)==42, 'Expected exactly 14 x 3 primary coverage rows')
    require(len({(r['method'],r['group']) for r in coverage})==42, 'Duplicate coverage keys')
    for row in coverage:
        require(row['group'] in GROUPS and int(row['records'])==GROUPS[row['group']] and
                int(row['complete_records'])==GROUPS[row['group']], f'Incomplete coverage: {row}')
        require(int(row['expected_residues'])==int(row['covered_residues']), f'Residue coverage gap: {row}')
    for method in EXPECTED:
        require(sum(int(r['covered_residues']) for r in coverage if r['method']==method)==961321,
                f'{method}: expected 961321 covered positions')

    plotting = load_json(ROOT/'figures/plot_audit.json')
    exact_methods(plotting['methods'], 'Plot audit')
    require(plotting['records']==3063 and plotting['primary_groups']==GROUPS, 'Plot cohort mismatch')
    bindings = {
        'plot_script_sha256':'plot_benchmark.py',
        'summary_sha256':'evaluation/group_summary.csv',
        'scoring_audit_sha256':'evaluation/scoring_audit.json',
        'training_exposure_audit_sha256':'references/training_exposure_protocol_sources_sha256.json',
    }
    for key, relative in bindings.items():
        require(plotting.get(key)==sha(ROOT/relative), f'Stale plot binding: {key}')
    require(utc(plotting['generated_utc'])>=utc(scoring['generated_utc']), 'Plot predates scoring')

    qa = load_json(ROOT/'QA_REPORT.json')
    require(qa['status']==QA_STATUS, f'Current QA is not accepted: {qa.get("status")}')
    exact_methods(qa['plot_methods'], 'QA plot methods')
    require(qa['qa_script_sha256']==sha(ROOT/'qa_current.py'), 'QA source changed since QA ran')
    require(qa['passed_custom_checks']==qa['total_custom_checks'], 'Custom QA has failures')
    require(qa['artifact_refresh']['source_hash_and_generation_bindings_valid'] is True,
            'QA reports stale render/source bindings')
    require(utc(qa['created_utc'])>=utc(plotting['generated_utc']), 'QA predates current render')
    for path, digest in qa['input_sha256'].items():
        source = Path(path)
        require(source.resolve().is_relative_to(ROOT.resolve()) and source.is_file(),
                f'QA input path is missing/outside task root: {path}')
        require(sha(source)==digest, f'QA input changed after inspection: {path}')
        add(source.relative_to(ROOT).as_posix())
    reviews = qa['visual_reviews']
    png_names = {'benchmark_current.png','benchmark_distributions_current.png'}
    require(len(reviews)==2 and {r['file'] for r in reviews}==png_names,
            'Both current PNGs require distinct actual visual reviews')
    for review in reviews:
        require(review['status']=='inspected_current_hash_with_notes' and review.get('observations'),
                f'Actual visual review is pending: {review["file"]}')
        require(review['sha256']==sha(ROOT/'figures'/review['file']), f'Visual SHA stale: {review["file"]}')
        require(utc(review['reviewed_utc'])>=utc(plotting['generated_utc']),
                f'Visual review predates render: {review["file"]}')

    sensitivity = load_json(ROOT/'evaluation/training_exposure_sensitivity_audit.json')
    require(sensitivity['output_rows']==60, 'Expected new 60-row secondary analysis, including BF1-seq')
    require(sensitivity['membership_sha256']==MEMBERSHIP_SHA and
            sha(ROOT/'references/training_exposure_membership_frozen.csv')==MEMBERSHIP_SHA,
            'Frozen training-exposure membership changed')
    require(sha(ROOT/'references/training_exposure_protocol_sources_sha256.json')==EXPOSURE_SHA,
            'Frozen training-exposure source protocol changed')
    for key, relative in {
        'source_scoring_audit_sha256':'evaluation/scoring_audit.json',
        'source_per_record_metrics_sha256':'evaluation/per_record_metrics.csv',
        'script_sha256':'evaluate_exposure_sensitivity.py',
        'protocol_sha256':'references/TRAINING_EXPOSURE_SENSITIVITY_PROTOCOL.json',
        'result_sha256':'evaluation/training_exposure_sensitivity.csv',
    }.items():
        require(sensitivity[key]==sha(ROOT/relative), f'Secondary analysis stale: {key}')
    require(len(sensitivity['eligibility'])==3 and all(x['eligible'] for x in sensitivity['eligibility']),
            'All three predeclared secondary comparators must now be eligible')
    expected_strata = {
        'DynaProt-M':{'all_atlas':1929,'train':1262,'val':39,'test':82,'unlisted':546},
        'backflip-1.0':{'all_atlas':1929,'train':1108,'val':137,'test':138,'unlisted':546},
        'backflip-1.0-seq':{'all_atlas':1929,'train':1108,'val':137,'test':138,'unlisted':546},
    }
    with (ROOT/'evaluation/training_exposure_sensitivity.csv').open(encoding='utf-8-sig',newline='') as handle:
        secondary_rows = list(csv.DictReader(handle))
    expected_keys = {(method, group, metric) for method, groups in expected_strata.items()
                     for group in groups for metric in ('spearman','top10_ap','mae','z_mae')}
    require(len(secondary_rows)==60 and
            {(r['identifier_method'],r['stratum'],r['metric']) for r in secondary_rows}==expected_keys,
            'Secondary analysis omits/duplicates a published category')
    for row in secondary_rows:
        require(row['dataset']=='atlas' and int(row['n_records'])==expected_strata[row['identifier_method']][row['stratum']],
                f'Secondary paired cohort mismatch: {row["identifier_method"]}/{row["stratum"]}')

    for tag in ('backflip-2.1','backflip-1.0','backflip-1.0-seq'):
        folder = ROOT/'backflip/results'/tag
        audit = load_json(folder/'completion_audit.json')
        require(audit['complete'] is True and audit['tag']==tag and
                audit['expected_records']==3063 and audit['expected_valid_residues']==961321 and
                audit['checked_output_hashes']==6126 and not audit['missing'] and not audit['failed'] and
                audit['singleton_shared_prediction_arrays_equal'] is True, f'{tag}: incomplete identity/hash audit')
        require(set(audit['contexts'])=={'whole_complex','independent_chain'}, f'{tag}: unexpected contexts')
        for context, item in audit['contexts'].items():
            require(item['completed_records']==3063 and item['valid_residues']==961321 and item['groups']==GROUPS,
                    f'{tag}/{context}: completion cohort mismatch')
        require(audit['contexts']['independent_chain']['shared_single_chain_records']==2193 and
                audit['contexts']['independent_chain']['actual_inference_records']==870,
                f'{tag}: incorrect single-chain shared-forward accounting')
        require(audit['audit_script_sha256']==sha(ROOT/'backflip/validate_outputs.py'),
                f'{tag}: identity validator source changed after completion')
        metadata = load_json(folder/'run_metadata.json')
        require(metadata['tag']==tag and metadata.get('completed_utc') and metadata['device']=='cpu',
                f'{tag}: completed local run metadata missing')
        require(metadata['runner_sha256']==sha(ROOT/'backflip/run_benchmark.py'), f'{tag}: runner SHA changed')

    fast = load_json(ROOT/'fastprotflex/full_audit.json')
    require(fast['mode']=='full' and fast['status']=='complete' and
            fast['records_requested']==3063 and fast['records_complete']==3063,
            'FastProtFlex full audit is incomplete')
    require(len(fast['results'])==3063 and len({r['sample_id'] for r in fast['results']})==3063 and
            all(r['status']=='complete' and r['absolute_RMSF_Angstrom'] is False and
                r['all_ca_predictions_finite'] is True and r['ca_identities_equal_frozen_truth_in_original_order'] is True
                for r in fast['results']), 'FastProtFlex per-record audit is incomplete or mislabels units')
    require(sum(int(r['ca_positions_verified']) for r in fast['results'])==961321,
            'FastProtFlex audited CA position count mismatch')
    fast_protocol = load_json(ROOT/'fastprotflex/run_protocol.json')
    fast_evidence = verify_fast_protocol(fast_protocol,fast,add,memory_v2_audits)
    if fast_evidence['algorithm_lifecycle_version']=='lifetime_v2':
        method = scoring['methods']['FastProtFlex']
        require(method['copied_target_atol']==1e-12 and
                0<=method['copied_target_max_abs_difference']<=method['copied_target_atol'],
                'Final independent scoring did not preserve the limited copied-target parser tolerance')
        fast_evidence['independent_final_scoring_audit_sha256'] = sha(ROOT/'evaluation/scoring_audit.json')

    config = load_json(ROOT/'long_queue_config.json')
    require(config.get('run_planet') is False and config.get('run_rmsfnet') is False,
            'Heavy models must remain deferred for this stage')
    controller = load_json(ROOT/'controller_state.json')
    require(controller['phase']=='finished_selected_cohort', 'Selected local queue has not finished successfully')
    exact_methods(controller['complete_entries'], 'Controller completed entries')
    require(not controller['missing_selected_entries'] and not controller['failed_phases'], 'Controller reports incomplete/failed phases')
    # The controller may truthfully finish before a person inspects the PNGs.
    # Preserve its then-current visual status; final visual acceptance is proved
    # independently by the current QA, exact PNG hashes and review timestamps.
    require(controller['deferred_heavy_models'].get('run_planet') is False and
            controller['deferred_heavy_models'].get('run_rmsfnet') is False, 'Controller heavy-model state differs')
    ledger = load_json(ROOT/'controller_phase_results.json')
    require(ledger and all(x['ok'] is True for x in ledger), 'Phase ledger reports a failure')
    needed = {'score_light_models','plot_light_models','report_light_models','exposure_light_models',
              'qa_light_models','inventory_light_models','publish_light_models'}
    require(needed <= {x['label'] for x in ledger}, 'Phase ledger lacks completed stage publication steps')
    require(not any(x['label'].startswith(('planet_','rmsfnet_full')) for x in ledger),
            'Heavy inference phases belong outside this selected light stage')
    for row in ledger:
        label = row['label']
        if label in needed-{'publish_light_models'}:
            add(f'controller_{label}.log')
    controller_completed_utc = utc(controller['updated_utc'])
    post_controller_visual_review = dict(
        controller_completed_utc=controller['updated_utc'],
        controller_visual_review_status_at_completion=controller.get('visual_review_status','not_recorded'),
        final_qa_created_utc=qa['created_utc'],final_qa_status=qa['status'],
        final_QA_after_controller_completion=utc(qa['created_utc'])>controller_completed_utc,
        actual_review_times={r['file']:r['reviewed_utc'] for r in reviews},
        controller_record_preserved_without_rewriting=True,
        basis='Final visual acceptance is the current hash-bound actual PNG review and QA, independent of the historical controller visual status.')

    exposure = load_json(ROOT/'references/training_exposure_protocol_sources_sha256.json')
    require(exposure['own_Full_Intra_OOF']['per_variant']['full']['records_with_same_family_in_train']==904,
            'Documented related-family exposure differs from expected 904 records')
    # Capture all initial hashes only after prerequisites, then verify them again
    # after copying. Stream hashes/copies never deserialize the aligned NPZ arrays.
    payload = []
    for row in sorted(sources.values(),key=lambda r:r['path']):
        payload.append(dict(path=row['path'],source=row['source'],sha256=sha(row['source']),
                            bytes=row['source'].stat().st_size))
    return payload, dict(base_scope_sha256=sha(BASE_SCOPE),base_files=len(base['files']),
        base_unique_files=base_unique_count,memory_v2_audits=fast_evidence['actual_validation_audits'],
        fastprotflex_protocol=fast_evidence,
        complete_methods=sorted(EXPECTED),qa_status=qa['status'],
        post_controller_visual_review=post_controller_visual_review,
        scoring_audit_sha256=sha(ROOT/'evaluation/scoring_audit.json'),
        visual_review_sha256={r['file']:r['sha256'] for r in reviews})


def freeze(payload, evidence):
    require(not DESTINATION.exists(), f'Refusing to overwrite existing snapshot: {DESTINATION}')
    staging = DESTINATION.parent/('.completed_14_entries.pending_'+uuid.uuid4().hex)
    staging.mkdir(parents=True,exist_ok=False)
    try:
        for row in payload:
            require(sha(row['source'])==row['sha256'], f'Source changed before copy: {row["path"]}')
            output = staging/row['path']
            output.parent.mkdir(parents=True,exist_ok=True)
            with row['source'].open('rb') as source, output.open('xb') as target:
                shutil.copyfileobj(source,target,length=1024**2)
            require(sha(output)==row['sha256'] and sha(row['source'])==row['sha256'],
                    f'Source/copy changed during copy: {row["path"]}')
        readme = '''# 已完成轻模型阶段：14项方法/配置

本快照在实际查看两张最终PNG、当前SHA绑定的QA通过、选中队列完成后冻结；全部14项均覆盖同一3,063条记录和961,321个有效RMSF位置。14项包括模型版本、输入配置与本文消融，不等于14个独立外部模型。

当前评价采用已发布权重与共同参考标签，训练条件混合，不能称为等训练数据竞赛。Full/Intra是记录/生物学分组OOF，并非序列家族留出；904条测试记录有相关family成员在训练中。DynaProt与BackFlip-1存在已文档化或官方候选清单的训练ID重叠，BF2没有可绑定清单。训练暴露的补充敏感性分析按预先冻结membership保留所有公开类别，不替代主队列。

BackFlip完整结构与逐链预测分别列出；单链共享同次前向不能记成额外独立计算。FastProtFlex、ANM/GNM为相对分数，不分配Å MAE。外部输出不充当同定义联合侧链熵Hχ。TCR–pMHC为单family，只作描述性点/分布，不给独立family置信区间。

PLANET、RMSF-net全量及其他重模型按用户选择后补，本快照不暗示已完成这些推理。QA保留原校验器的局限与当前人工视觉观察；不把自动校验当作未来图的视觉通过。论文主文件未由本冻结脚本改写。

SNAPSHOT_SCOPE.json说明范围，INDEPENDENT_MANIFEST.json保存全部载荷和范围文件SHA，INDEPENDENT_MANIFEST.sha256提供独立清单的校验值。原始模型目录与权重未整目录打包；aligned_predictions.npz仅作字节复制而未加载。冻结脚本不评分、不推理、不改变源文件。
'''
        if evidence['fastprotflex_protocol']['algorithm_lifecycle_version']=='lifetime_v2':
            readme += '\nFastProtFlex实际采用lifetime_v2：仅调整R临时对象释放与垃圾回收时机，保留官方数值程序；三来源预定样本的全原子身份、坐标与精确RDS预测逐值相同，最大差为0。实际入口、官方原件、候选、原v1协议、真实审计与内存校准均按SHA绑定。\n'
            readme += '资源策略明确版本化为物理可用内存至少1.15倍预测峰值加256MiB，并额外要求Windows系统可承诺内存至少预测峰值加1024MiB；commit计数未知时禁止启动。历史512MiB规则及原三源实测记录保留原样，不能回写为新规则。缓存复用时验证输入、来源、全部输出SHA与PDB到NPZ残基身份；本次运行未重新核对缓存target/valid/residue_index的全部真值，该项由当前源码绑定且全量覆盖的最终独立评分完成。准备事件、工作集页缓存整理事件与资源策略准备记录按实际存在的原件纳入，不能视为模型数值验证或绕过门槛。\n'
            if evidence['fastprotflex_protocol']['resource_extension']:
                readme += '8058原子7zh9_A极值探针只验证实际内存、输入输出身份与有限性，缺少原版精度输出，未进行原版逐值等价比较。四点资源校准保留原三点基准及两份独立审计，不能把资源极值探针称为第四个数值等价样本。\n'
        (staging/'STAGE_README_zh.md').write_text(readme,encoding='utf-8')
        public = [{k:v for k,v in r.items() if k!='source'} for r in payload]
        public.append(dict(path='STAGE_README_zh.md',sha256=sha(staging/'STAGE_README_zh.md'),
                           bytes=(staging/'STAGE_README_zh.md').stat().st_size))
        scope = dict(created_utc=datetime.now(timezone.utc).isoformat(),completed_entries=14,
            records=3063,valid_rmsf_positions=961321,original_root=str(ROOT),
            scope='Completed selected light stage; 14 method/configuration entries, mixed training exposure, record/biological-group OOF with 904 related-family-exposed test records. Heavy full queues deferred.',
            source_hashes_checked_before_and_after_copy=True,arrays_deserialized=False,
            evidence=evidence,files=public)
        (staging/'SNAPSHOT_SCOPE.json').write_text(json.dumps(scope,indent=2,ensure_ascii=False),encoding='utf-8')
        manifest_files = public+[dict(path='SNAPSHOT_SCOPE.json',sha256=sha(staging/'SNAPSHOT_SCOPE.json'),
                                     bytes=(staging/'SNAPSHOT_SCOPE.json').stat().st_size)]
        manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(),snapshot_name=DESTINATION.name,
                        algorithm='sha256',files=manifest_files,
                        excludes_self_and_checksum_to_avoid_circular_hashes=True)
        manifest_path = staging/'INDEPENDENT_MANIFEST.json'
        manifest_path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False),encoding='utf-8')
        manifest_sha = sha(manifest_path)
        (staging/'INDEPENDENT_MANIFEST.sha256').write_text(manifest_sha+'  INDEPENDENT_MANIFEST.json\n',encoding='ascii')
        for row in manifest_files:
            require(sha(staging/row['path'])==row['sha256'], f'Final payload mismatch: {row["path"]}')
        for row in payload:
            require(sha(row['source'])==row['sha256'], f'Source changed before final seal: {row["path"]}')
        for row in evidence['fastprotflex_protocol']['external_sources']:
            require(sha(row['path'])==row['sha256'], f'Official external source changed during freeze: {row["path"]}')
        require(sha(BASE_SCOPE)==evidence['base_scope_sha256'], 'Previous scope changed during freeze')
        require(not DESTINATION.exists(), f'Snapshot appeared concurrently: {DESTINATION}')
        staging.rename(DESTINATION)
    except Exception as error:
        (staging/'SNAPSHOT_INCOMPLETE.json').write_text(json.dumps(dict(status='incomplete_not_a_snapshot',
            error=str(error),source_files_untouched=True),indent=2),encoding='utf-8')
        raise RuntimeError(f'Freeze rejected; diagnostic staging preserved at {staging}: {error}') from error
    print(json.dumps(dict(status='frozen',snapshot=str(DESTINATION),payload_files=len(public),
                          independent_manifest_sha256=manifest_sha),indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--freeze',action='store_true',help='Copy only after all hard prerequisites pass')
    mode.add_argument('--check-only',action='store_true',help='Read-only prerequisite/hash check (default)')
    parser.add_argument('--memory-v2-audit',action='append',default=[],
                        help='Optional task-root-relative actual adopted v2 audit; must match the protocol exactly')
    args = parser.parse_args()
    payload,evidence = verify_prerequisites(args.memory_v2_audit)
    if args.freeze:
        freeze(payload,evidence)
    else:
        print(json.dumps(dict(status='prerequisites_pass_no_snapshot_written',files=len(payload),
                              destination=str(DESTINATION),evidence=evidence),indent=2))


if __name__=='__main__':
    try:
        main()
    except Exception as error:
        print(f'FREEZE REFUSED: {error}',file=sys.stderr)
        sys.exit(1)
