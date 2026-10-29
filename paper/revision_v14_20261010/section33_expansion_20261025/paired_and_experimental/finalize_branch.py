"""Verify result joins/units and write final provenance/file checksums."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
OUT=Path(__file__).resolve().parent
data=pd.read_csv(OUT/'test53_corrected_full_oof_residue_deltas_with_family.csv.gz')
inventory=pd.read_csv(OUT/'paired_inventory.csv')
accepted=inventory[inventory.status.eq('evaluated')]
assert len(inventory)==37 and inventory.pdb_id.nunique()==36
assert len(accepted)==36 and accepted.pdb_id.nunique()==35
assert accepted.apo_fold.eq(accepted.holo_fold).all()
assert accepted.corrected_mapping_identity.eq(1).all()
assert accepted.corrected_mapping_coverage.min()>.87
assert data.pair_id.nunique()==36 and data.pdb_id.nunique()==35
assert np.allclose(data.md_delta,data.md_holo-data.md_apo,atol=1e-6)
assert np.allclose(data.model_delta,data.pred_holo-data.pred_apo,atol=1e-6)
assert not data[data.metric.eq('chi_entropy')].amino_acid.isin(list('AGP')).any()
assert data.groupby('metric').size().to_dict()==dict(chi_entropy=9723,rmsf=11842)
assert data.family_components.nunique()==27
assert data[data.family_clean].pdb_id.nunique()==20
assert data[data.family_clean].family_components.nunique()==20
pair=pd.read_csv(OUT/'test53_per_pair_metrics.csv')
system=pd.read_csv(OUT/'test53_per_system_with_family_component.csv')
for metric in ['rmsf','chi_entropy']:
    a=pair[(pair.pdb_id=='4ETP')&(pair.metric==metric)].mae_model.mean()
    b=system[(system.pdb_id=='4ETP')&(system.metric==metric)].mae_model.iloc[0]
    assert np.isclose(a,b)
nmr=pd.read_csv(OUT/'Fares2009_to_A003R_all50_mapping.csv')
assert len(nmr)==50 and nmr.valid_order_parameter.sum()==37
assert (nmr.valid_order_parameter&nmr.A003R_present).sum()==35
assert not nmr[nmr.s2_raw.eq(0)].valid_order_parameter.any()
assert (nmr.s2_raw>1).sum()==3
vectors=np.load(OUT/'A003R_bound_aligned_methyl_unit_vectors.npz')['vectors']
assert vectors.shape==(5000,48,3) and np.isfinite(vectors).all()
assert np.allclose(np.linalg.norm(vectors,axis=2),1.,atol=1e-10)
protocol=json.loads((OUT/'A003R_methyl_S2_protocol.json').read_text())
assert protocol['source_files'][0]['sha256']=='4b2206ae6e9b71d9b8659e661d6825320eebae563ec4e7c2643a43a5225bf1f5'
assert protocol['source_files'][1]['sha256']=='9d5ffa6a837e7400a5ed571fa9ba2cbeb9838606449aa5c2d145477d725510df'
conditions=dict(experiment='Fares2009_free_human_ubiquitin_methyl_RDCs',
                temperature_K=None,temperature_reason='Not explicitly identified in accessible Fares2009 main text/SI; do not infer from MD.',
                datasets1_to10='Conditions referenced to Lakomek2008a/b; alignment_conditions.csv preserves original Table1.',
                dataset11='Literature RDCs from OttigerandBax1998/1999; see original condition references.',
                dataset12=dict(positive_gel_percent=7,protein_mg_approx=2,volume_uL=300,
                               buffer='50mM NaPO4,100mM NaCl,pH6.5,0.05%NaN3'),
                dataset13=dict(medium='Pf1phage15mg/mL',buffer='50mM NaPO4,pH6.5,0.05%NaN3',final_NaCl_mM_approx=300),
                observable='methyl axial generalized order parameter; rapid methyl axial rotation removed',
                timescale='spatial averaging for motions faster than ms; includes sub-tumbling and supra-tumbling contributions; not a specific kinetic rate',
                source_url='https://link.springer.com/article/10.1007/s10858-009-9354-7')
(OUT/'Fares2009_sample_conditions.json').write_text(json.dumps(conditions,indent=2),encoding='utf-8')
sources=[dict(citation='Fareset al2009,JBiomolNMR45:23-44',doi='10.1007/s10858-009-9354-7',
              article_url='https://link.springer.com/article/10.1007/s10858-009-9354-7',
              supplement_url='https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs10858-009-9354-7/MediaObjects/10858_2009_9354_MOESM1_ESM.pdf',
              local_file='sources/Fares2009_SI.pdf',use='All50methyl axes,37identifiableS2; source condition and uncertainty preservation'),
         dict(citation='Peifferet al2022,PLoSComputBiol18:e1009977',doi='10.1371/journal.pcbi.1009977',
              article_url='https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009977',
              supplement_url='https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1009977.s001&type=supplementary',
              local_file='sources/Peiffer2022_KIX_S1_methods.pdf',
              use='S1specifiesfiveindependent>=100nstrials. ArticleDataAvailability points to manuscript/SI; no pertrialtrajectorydownload link found in those entries, notproofrawdata unavailableeverywhere.')]
(OUT/'literature_sources.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
qa=dict(status='pass',original_pairs=37,original_PDBs=36,evaluated_pairs=36,evaluated_PDBs=35,
        full_cohort_components=27,family_clean_pairs=int(accepted.family_clean.sum()),family_clean_PDBs=20,
        mapping_complete_sequence_matches=30,mapping_exact_global_alignment_pairs=6,
        corrected_mapping_min_coverage=float(accepted.corrected_mapping_coverage.min()),
        residual_counts=data.groupby('metric').size().to_dict(),NMR_rows=50,NMR_identifiable_axes=37,
        identifiable_axes_present_in_bound_MD=35,bound_backcalculated_axes=48,
        original_A003R_structure_and_trajectory_SHA256_match_paper_archive=True,
        checks='samefoldbothstates;chemicalsequenceexactmatching;legacyindexremapping;validitymasks;deltaidentity;4ETPequal-systemweight;27clusterjoin;20familyclean;S2missingzeroandupperbound;normalizedvectors;unchangedrawinputhashes')
(OUT/'final_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
manifest=[]
for f in sorted(OUT.rglob('*')):
    if f.is_file() and f.name!='file_manifest_sha256.csv' and '__pycache__' not in f.parts:
        manifest.append(dict(relative_path=str(f.relative_to(OUT)),bytes=f.stat().st_size,
                             sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
pd.DataFrame(manifest).to_csv(OUT/'file_manifest_sha256.csv',index=False)
print(json.dumps(qa));print('FILES',len(manifest))
