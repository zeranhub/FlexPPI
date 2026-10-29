"""Numerical sanity tests for post-readout scale-free sensitivity metrics."""
import unittest
import numpy as np
from relative_pilot import transform, summarize


class RelativePilotTest(unittest.TestCase):
    def test_positive_affine_invariance(self):
        values = np.array([1., 2., 4., 9., 16.])
        np.testing.assert_allclose(transform(values,'percentile'),transform(3*values+17,'percentile'))
        np.testing.assert_allclose(transform(values,'standardized'),transform(3*values+17,'standardized'))

    def test_rank_ties_and_range(self):
        np.testing.assert_allclose(transform([1.,1.,2.,3.],'percentile'),[100/6,100/6,200/3,100])

    def test_global_shift_not_a_rank_response(self):
        before=np.array([1.,2.,3.,4.]);after=before+9
        self.assertTrue(np.all(after-before == 9))
        np.testing.assert_array_equal(transform(after,'percentile')-transform(before,'percentile'),np.zeros(4))

    def test_constant_correlation_is_undefined(self):
        result=summarize([0.,0.,0.],[1.,2.,3.])
        self.assertTrue(np.isnan(result['pearson']) and np.isnan(result['spearman']))
        self.assertEqual(result['mae'],2.)

    def test_unchanged_baseline(self):
        observed=np.array([-10.,0.,20.])
        self.assertEqual(summarize(observed,np.zeros(3))['mae'],10.)


if __name__=='__main__':
    unittest.main()
