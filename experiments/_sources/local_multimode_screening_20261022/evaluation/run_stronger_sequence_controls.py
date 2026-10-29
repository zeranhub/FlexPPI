"""Additional stronger baselines after initial results; no new dynamic method.

The primary study is unchanged. This adds sequence-only feedback and BLOSUM-only
comparators to test whether an apparent gain over sequence+static merely repairs
a weak static block. Same labels, budgets, policies and seeds as the primary.
"""
import numpy as np
from Bio.Align import substitution_matrices
import evaluate_multimode as base

original=base.construct_kernels
def stronger(sequences,static,R,H,Rscale,Hscale):
    kernels,settings=original(sequences,static,R,H,Rscale,Hscale)
    aa='ACDEFGHIKLMNPQRSTVWY';matrix=substitution_matrices.load('BLOSUM62')
    X=np.asarray([[matrix[a,b] for a in s for b in aa] for s in sequences])
    chem,setting=base.distance_kernel(base.pairwise_sq(X))
    kernels['sequence_only_control']=kernels['sequence']
    kernels['BLOSUM_only_control']=chem
    settings['BLOSUM_only_control']=setting
    return kernels,settings
base.construct_kernels=stronger
base.EXTRA_STRATEGIES=['sequence_only_control','BLOSUM_only_control']
if __name__=='__main__':base.main()
