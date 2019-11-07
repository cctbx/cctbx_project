#!/usr/bin/env python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function
from scitbx import matrix
from scitbx.array_family import flex
import random
from scitbx.math.lefebvre import matrix_inverse_error_propagation
from six.moves import range
from six.moves import map

"""Implementation of the propagation of errors formula for matrix inversion
given in Lefebvre et al. (1999) http://arxiv.org/abs/hep-ex/9909031. As in
that paper the analytical formula is tested against Monte Carlo simulation."""

def random_vector_close_to(vector):
  vec = matrix.col(vector)
  vec2 = vec.rotate_around_origin(matrix.col(
              (random.random(),
               random.random(),
               random.random())).normalize(),
               random.gauss(0, 5), deg = True)
  length_multiplier = max(0.5, random.gauss(1, 0.1))
  return vec2 * length_multiplier

def create_Bmat():
  """Create a random crystal model, in P1 for maximum generality. Return the
  B matrix of this crystal"""

  from dxtbx.model import Crystal

  vecs = list(map(random_vector_close_to,
             [(20, 0, 0),
              (0, 20, 0),
              (0, 0, 20)]))

  return matrix.sqr(Crystal(*vecs, space_group_symbol = "P 1").get_B())

def create_sig_mat(mat, fraction_sd=0.01):
  """Create a matrix of errors of the elements of mat, where each error is a
  normal deviate with a sigma some fraction of the absolute element value"""

  vals = [random.gauss(0.0, abs(fraction_sd*e)) for e in mat]
  return matrix.sqr(vals)

def perturb_mat(mat, fraction_sd=0.01):
  """Perturb the elements of mat by normal deviates with sigmas of some
  fraction of the absolute element values."""

  return mat + create_sig_mat(mat, fraction_sd)

def cov(a, b):
  """Return the sample covariance of vectors a and b"""
  a = flex.double(a)
  b = flex.double(b)
  n = len(a)
  assert n == len(b)
  resid_a = a - flex.mean(a)
  resid_b = b - flex.mean(b)
  return flex.sum(resid_a*resid_b) / (n - 1)

def calc_monte_carlo_covariances(mats):
  """Given the sequence of matrices mats, calculate the variance-covariance
  matrix of the elements"""

  # check input
  assert all([e.is_square() for e in mats])
  n = mats[0].n_rows()
  assert all([e.n_rows() == n for e in mats])

  # create an empty var-cov matrix
  covmat = flex.double(flex.grid(n**2,n**2), 0.0)

  for i in range(covmat.all()[0]):
    for j in range(covmat.all()[1]):
      a = [m[i] for m in mats]
      b = [m[j] for m in mats]
      covmat[i,j] = cov(a,b)

  return covmat

def calc_monte_carlo_population_covariances(mats, mean_matrix):
  """Given the sequence of matrices mats, calculate the variance-covariance
  matrix of the elements, using the known mean values for each elt in mat
  to avoid the approximation implied in taking the sample covariance"""

  # check input
  assert all([e.is_square() for e in mats])
  n = mats[0].n_rows()
  assert all([e.n_rows() == n for e in mats])

  # create an empty var-cov matrix
  covmat = flex.double(flex.grid(n**2,n**2), 0.0)

  for i in range(covmat.all()[0]):
    for j in range(covmat.all()[1]):
      resid_a = flex.double([m[i] - mean_matrix[i] for m in mats])
      resid_b = flex.double([m[j] - mean_matrix[j] for m in mats])
      covmat[i,j] = flex.mean(resid_a*resid_b)

  return covmat

def test_lefebvre():
  """Run the test presented in part 4 of the paper Lefebvre et al. (1999),
  http://arxiv.org/abs/hep-ex/9909031."""

  # Simple test following the paper. First MC sim
  mat = matrix.sqr((0.7,0.2,0.4,0.6))
  sig_mat = mat*0.01
  p_mats = [perturb_mat(mat) for i in range(10000)]
  inv_p_mats = [m.inverse() for m in p_mats]

  # Here use the more exact formula for population covariance as we know the
  # true means rather than sample covariance where we calculate the means
  # from the data
  cov_inv_mat_MC = calc_monte_carlo_population_covariances(inv_p_mats,
    mean_matrix=mat.inverse())

  # Now analytical formula
  cov_mat = matrix.diag([e**2 for e in sig_mat])

  cov_inv_mat = matrix_inverse_error_propagation(mat, cov_mat)
  cov_inv_mat = cov_inv_mat.as_flex_double_matrix()

  # Get fractional differences
  frac = (cov_inv_mat_MC - cov_inv_mat) / cov_inv_mat_MC

  # assert all fractional errors are within 5%. When the random seed is allowed
  # to vary this fails about 7 times in every 100 runs. This is mainly due to
  # the random variation in the MC method itself, as with just 10000 samples
  # one MC estimate compared to another MC estimate regularly has some
  # fractional errors greater than 5%.
  assert all([abs(e) < 0.05 for e in frac])

  return

def test_B_matrix():
  """Test errors in an inverted B matrix, when the errors in B are known (and
  are independent)"""

  Bmat = create_Bmat()
  inv_Bmat = Bmat.inverse()

  # Note that elements in the strict upper triangle of B are zero
  # so their perturbed values will be zero too and the errors in these elements
  # will also be zero (this is what we want)

  # calculate the covariance of elements of inv_B by Monte Carlo simulation,
  # assuming that each element of B has an independent normal error given by a
  # sigma of 1% of the element value
  perturbed_Bmats = [perturb_mat(Bmat, fraction_sd=0.01) for i in range(10000)]
  invBs = [m.inverse() for m in perturbed_Bmats]
  cov_invB_MC = calc_monte_carlo_population_covariances(invBs, inv_Bmat)

  # Now calculate using the analytical formula. First need the covariance
  # matrix of B itself. This is the diagonal matrix of errors applied in the
  # simulation.
  n = Bmat.n_rows()
  sig_B = Bmat * 0.01
  cov_B = matrix.diag([e**2 for e in sig_B])

  # Now can use the analytical formula
  cov_invB = matrix_inverse_error_propagation(Bmat, cov_B)
  cov_invB = cov_invB.as_flex_double_matrix()

  # Get fractional differences
  frac = flex.double(flex.grid(n**2, n**2), 0.0)
  for i in range(frac.all()[0]):
    for j in range(frac.all()[1]):
      e1 = cov_invB_MC[i,j]
      e2 = cov_invB[i,j]
      if e1 < 1e-10: continue # avoid divide-by-zero errors below
      if e2 < 1e-10: continue # avoid elements that are supposed to be zero
      frac[i,j] = (e1 - e2) / e1

  # assert all fractional errors are within 5%
  assert all([abs(e) < 0.05 for e in frac])

  return

if __name__ == '__main__':

  # set random seet to avoid rare test failures
  random.seed(42)

  # run the Lefebvre paper test
  test_lefebvre()
  print("OK")

  # run a similar test for a B matrix
  try:
    import dxtbx # import dependency
  except ImportError:
    pass
  else:
    test_B_matrix()
    print("OK")
