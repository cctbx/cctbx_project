#
# lefebvre.py
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from six.moves import range

def matrix_inverse_error_propagation(mat, cov_mat):
  """Implement analytical formula of Lefebvre et al. (1999)
  http://arxiv.org/abs/hep-ex/9909031 to calculate the covariances of elements
  of mat^-1, given the covariances of mat itself. The input covariance matrix,
  and the return value of the function, have elements ordered by treating mat as
  a 1D vector using row-major ordering."""

  # initialise covariance matrix
  assert mat.is_square()
  assert cov_mat.is_square()
  n = mat.n_rows()
  assert cov_mat.n_rows() == n**2

  # use flex for nice 2D indexing
  inv_mat = flex.double(mat.inverse())
  inv_mat.reshape(flex.grid(n, n))
  cov_mat = cov_mat.as_flex_double_matrix()

  inv_cov_mat = flex.double(flex.grid(n**2,n**2), 0.0)
  for alpha in range(n):
    for beta in range(n):
      for a in range(n):
        for b in range(n):

          # index into inv_cov_mat after flattening inv_mat
          u = alpha * n + beta
          v = a * n + b
          # skip elements in the lower triangle
          if v < u: continue

          # The element u,v of the result is the calculation
          # cov(m^-1[alpha, beta], m^-1[a, b])
          elt = 0.0
          for i in range(n):
            for j in range(n):
              for k in range(n):
                for l in range(n):

                  # index into cov_mat after flattening mat
                  x = i * n + j
                  y = k * n + l
                  elt += inv_mat[alpha, i] * inv_mat[j, beta] * \
                       inv_mat[a, k] * inv_mat[l, b] * \
                       cov_mat[x, y]
          inv_cov_mat[u, v] = elt

  inv_cov_mat.matrix_copy_upper_to_lower_triangle_in_place()
  return inv_cov_mat.as_scitbx_matrix()

