import sys
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts')]
from run_apo_holo_exploration import subset, KEYS, align_states
from analyze_apo_holo_exploration import loops, corr


def test_subset_preserves_coordinates_and_sequence_positions():
    example={key:torch.arange(6) for key in KEYS}
    example['chain_index']=torch.tensor([0,0,1,1,2,2])
    example['backbone']=torch.arange(72).reshape(6,4,3)
    result=subset(example,np.array([2,3,4,5]))
    assert torch.equal(result['backbone'],example['backbone'][2:])
    assert result['chain_index'].tolist()==[0,0,1,1]
    assert result['residue_index'].tolist()==[2,3,4,5]
    result['backbone'][0,0,0]=-1
    assert example['backbone'][2,0,0]!=-1


def test_alignment_rejects_unrelated_chain():
    a=dict(chain=np.array(['A']*20),aa=np.array(list('A'*20)))
    b=dict(chain=np.array(['B']*20),aa=np.array(list('V'*20)))
    with pytest.raises(ValueError):
        align_states(a,b,['B'])


def test_alignment_returns_only_exact_residues():
    a=dict(chain=np.array(['A']*20),aa=np.array(list('ACDEFGHIKLMNPQRSTVWY')))
    b=dict(chain=np.array(['B']*20),aa=np.array(list('ACDEFGHIKLMNPQRSTVWA')))
    _,ia,ib,quality,_=align_states(a,b,['B'])
    assert len(ia)==19
    assert np.array_equal(a['aa'][ia],b['aa'][ib])
    assert quality[0]['identity']==.95


def test_star_loop_requires_complete_records():
    text='loop_\n_Order_param.Seq_ID\n_Order_param.Order_param_val\n1 0.8\n2 0.4\nstop_\n'
    assert list(loops(text,'_Order_param.'))[1]['Order_param_val']=='0.4'
    with pytest.raises(ValueError):
        list(loops(text.replace('2 0.4','2'),'_Order_param.'))


def test_constant_predictor_has_undefined_correlation():
    assert np.isnan(corr([0,0,0,0],[1,2,3,4]))
