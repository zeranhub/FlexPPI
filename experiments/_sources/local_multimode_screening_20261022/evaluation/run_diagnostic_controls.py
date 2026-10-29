"""Pre-evaluation controls: stronger sequence chemistry and shuffled profiles.

These complement the frozen primary analysis; they cannot replace a negative
primary finding. Sequence-to-profile correspondence is broken only in placebo
controls. The genuine dynamic profiles are never shuffled in the primary run.
"""
import numpy as np
from Bio.Align import substitution_matrices
import evaluate_multimode as base

original=base.construct_kernels
def with_controls(sequences,static,R,H,Rscale,Hscale):
    kernels,settings=original(sequences,static,R,H,Rscale,Hscale)
    n=len(sequences)
    matrix=substitution_matrices.load('BLOSUM62')
    aa='ACDEFGHIKLMNPQRSTVWY'
    X=np.asarray([[matrix[a,b] for a in s for b in aa] for s in sequences])
    kchem,setting=base.distance_kernel(base.pairwise_sq(X));settings['BLOSUM62_chemistry']=setting
    scale=np.sqrt(np.mean(np.asarray(static)**2,axis=0));scale[scale==0]=1.
    ks,_=base.distance_kernel(base.pairwise_sq(np.asarray(static)/scale))
    kr,_=base.distance_kernel(base.pairwise_sq(np.asarray(R)/Rscale))
    kh,_=base.distance_kernel(base.pairwise_sq(np.asarray(H)/Hscale))
    kernels['BLOSUM_static']=(kchem+ks)/2
    kernels['BLOSUM_static_RH']=(kchem+ks+kr+kh)/4
    seq_kernel=kernels['sequence']
    for k in range(10):
        permutation=np.random.default_rng(2026102200+k).permutation(n)
        name=f'baseline_RH_shuffled{k:02d}'
        kernels[name]=(seq_kernel+ks+kr[np.ix_(permutation,permutation)]+kh[np.ix_(permutation,permutation)])/4
        settings[name]=dict(permutation_seed=2026102200+k,R_H_shuffled_together=True,
            destroys_candidate_correspondence=True,hidden_labels_used=False)
    return kernels,settings

base.construct_kernels=with_controls
base.EXTRA_STRATEGIES=['BLOSUM_static','BLOSUM_static_RH']+[f'baseline_RH_shuffled{k:02d}' for k in range(10)]
if __name__=='__main__':base.main()
