"""Preserve independent v54 checks as a read-only baseline for v55."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent.parent
OLD=ROOT.parent/'mel5_composite_v54_20261022/figure4'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
names=['source_CDR_absolute.csv','source_peptide_four_states.csv','source_CDR_three_contrasts.csv','source_interface_three_contrasts.csv','source_CDR_heatmap_all135.csv','source_interface_heatmap_all90.csv','source_heatmap_arrays.npz']
baseline={
 'status':'ready_for_new_figure','v54_reviewed_png_SHA256':'514bc36ab0ec79c7348f1c7ccc804212906ab927334527c0e380c4e45cf63d87',
 'authorized_visual_changes':['Remove the ranked23 Hchi bar panel.','Combine the old e/f heatmaps as one new panel d with upper and lower metric blocks.','Reflow existing a/b/c panels.'],
 'must_preserve':{'peptide_four_states_values':40,'absolute_CDR_two_state_values':90,'CDR_delta_values':135,'pMHC_interface_delta_slots':90,'interface_Hchi_valid_per_row':[23,23,22],'common_CDR_sites':45,'pMHC_union_sites':30},
 'plotting_ranked23_removed_but_source_retention_allowed':True,
 'sources':[{'filename':n,'path':str(OLD/n),'SHA256':sha(OLD/n)} for n in names],
 'minimal_prose_edits':{'P91':'Keep all23 Hchi/Q72-largest-decrease numerical statement; point to new Figure4d lower block.','P93':'Replace Figure4e/f with Figure4d upper/lower block.','caption':'Remove old rankedbar description and p-prefix explanation. Rename old e/f to d with upper CDR deltaRMSF and lower pMHC deltaHchi; preserve all contrast, unit and eligibility definitions.'}}
(ROOT/'source_audit/source_preservation_baseline.json').write_text(json.dumps(baseline,ensure_ascii=False,indent=2),encoding='utf-8')
print('Read-only v54 source baseline recorded.')
