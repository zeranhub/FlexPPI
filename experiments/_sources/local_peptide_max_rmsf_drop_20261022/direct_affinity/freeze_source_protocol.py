"""Record the user-specified exploratory score before calculating its utility."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parents[1] / 'local_antigen_screening_20261021'
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
files = [OLD/'evaluation_review/fixed_analysis_protocol.json',
         OLD/'experimental_labels/label_rule_freeze.json',
         OLD/'structures/fixed_backbone_manifest.json']
for variant in ['reported_hybrid', 'all_constrained', 'all_fitted']:
    files += [OLD/f'experimental_labels/SPR_labels_candidates_{variant}.csv',
              OLD/f'experimental_labels/SPR_labels_observations_{variant}.csv']
protocol = {
  'protocol_id': 'user_peptide_max_rmsf_drop_direct_SPR_exploratory_20261022_v1',
  'created_at_utc': datetime.now(timezone.utc).isoformat(),
  'status': 'source_only_before_new_score_utility_calculation',
  'user_hypothesis': 'Prefer a mutant if at least one peptide position has a sufficiently large predicted RMSF decrease.',
  'knowledge_at_recording': {
    'previous_affinity_labels_and_previous_negative_mean_feature_results_known': True,
    'predictions_previously_generated_and_archived': True,
    'parent_prediction_header_and_parent_audit_read': True,
    'new_max_drop_candidate_scores_or_their_affinity_utility_not_yet_calculated': True,
    'not_blinded_confirmatory_or_prospective': True},
  'score': {
    'formula': 'max over all 9 peptide positions of (RMSF_parent - RMSF_mutant)',
    'unit': 'Angstrom, original predicted RMSF physical scale',
    'sign': 'Larger score is hypothesized to indicate stronger binding; retain signed values, including all-negative profiles.',
    'sites': 'All nine peptide residues, including the mutated position; geometry and sequence verified against archived metadata.',
    'no_Hchi_no_training_no_new_inference': True,
    'threshold': 'No absolute threshold calibrated or optimized; analyze continuous ranking and archive every candidate.'},
  'inherited_fixed_choices': {
    'primary_templates': {'1G4_NYESO':'2BNQ', 'A6_TAX':'1AO7'},
    'sensitivity_templates': {'1G4_NYESO':['2BNR'], 'A6_TAX':['1QSE','1QSF']},
    'primary_SPR_version': 'reported_hybrid',
    'SPR_versions': ['reported_hybrid','all_constrained','all_fitted'],
    'SPR_rule': 'Pettmann original Methods: constrained KD >20 uM uses constrained, otherwise Bmax-fitted. No score-based source selection.',
    'parents': {'1G4_NYESO':'SLLMWITQV', 'A6_TAX':'LLFGYPVYV'},
    'label': 'pKD=-log10(KD in M); candidate-minus-parent pKD within the same source study and assay at 37C, then median across studies for a duplicate candidate.',
    'batch': 'Published study/protocol matched means; individual run/batch pairing not reported.',
    'censoring': 'No bound is treated as exact. Inherited label censoring is checked; any censored item is retained separately, not silently dropped.'},
  'primary_analysis': {
    'unit': '26 unique single mutants: 13 per TCR panel; parents excluded from correlations.',
    'outcome': 'archived within-study parent-relative delta_pKD, larger means stronger binding',
    'statistic': 'Descriptive Spearman correlation of max RMSF drop versus delta_pKD on all 26 primary-template single mutants, and each 13-candidate panel separately.',
    'pooled_caution': 'Two TCR panels are not 26 independent systems; pooled value is descriptive and panels must remain visible.',
    'no_p_value_or_population_CI': True,
    'no_deletion_flipped_direction_or_best_template_selection': True},
  'sensitivity_analysis': {
    'scope': 'All five template panels x all three SPR versions; same single mutant candidate scope within each panel; per-study correlations also retained.',
    'multiple_mutants': 'Two A6 double mutants shown separately with raw scores and labels; never merge into primary n.',
    'not_independent_biological_repeats': True},
  'interpretation_limits': [
    'Small, retrospectively assembled affinity panels, not complete mutational libraries; cannot estimate screening hit rate.',
    'Only one primary candidate is slightly stronger than its parent and all A6 mutants are weaker: cannot establish high-affinity design success.',
    'Fixed parental backbone and historical 2431 no-ITrace joint checkpoint; not the manuscript 2430 checkpoint.',
    'Correlation cannot establish a causal RMSF-to-affinity mechanism, bound-to-unbound entropy change, or universal direction.'],
  'source_files_sha256': {str(p):sha(p) for p in files}}
ROOT.mkdir(parents=True, exist_ok=True)
target=ROOT/'source_only_protocol.json'
if target.exists():
    raise FileExistsError(target)
target.write_text(json.dumps(protocol,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'protocol':str(target),'sha256':sha(target)},ensure_ascii=False))
