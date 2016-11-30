"""Unit test for mod_leastsqr.py"""
from __future__ import division
import unittest
import numpy as np
from prime.postrefine import mod_leastsqr
from cctbx.array_family import flex
import cPickle as pickle

class KnownValues(unittest.TestCase):
  x = np.array([0,1,2,3])
  y = np.array([-1, 0.2, 0.9, 2.1])
  A = np.vstack([x, np.ones(len(x))]).T
  r_sqr = 0.990099009901
  se = 0.158113883008
  n_params = 2
  m, c = np.linalg.lstsq(A, y)[0]
  y_model = (m * x) + c
  pickle_filename = 'test_data/int_shot-s00-E1_0_00002_20140518051704255_masked.pickle'
  test_pickle = pickle.load(open(pickle_filename, 'rb'))

  def testCoefficientOfDetermination(self):
    r_sqr = mod_leastsqr.coefficient_of_determination(flex.double(self.y), flex.double(self.y_model))
    self.assertAlmostEqual(self.r_sqr, r_sqr)

  def testStandardErrorOfTheEstimate(self):
    se = mod_leastsqr.standard_error_of_the_estimate(flex.double(self.y), flex.double(self.y_model), self.n_params)
    self.assertAlmostEqual(self.se, se)

  def testLeastSqrGetFilteredData(self):


if __name__ == "__main__":
  unittest.main()
