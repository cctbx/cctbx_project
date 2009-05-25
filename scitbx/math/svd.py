from __future__ import division
import libtbx
import scitbx.math
from scitbx.array_family import flex

class real(object):
  """ SVD decomposition of a real matrix.

      References:
        [1] Tony F. Chan,
            An improved algorithm for computing the singular value
            decomposition,
            ACM Trans. Math. Soft. 8 (1982), no. 1, 72-83.
  """

  crossover = 5/3
  """ crossover from [1, table II] (note: LAPACK uses m > 1.6*n) """

  def __init__(self, a, accumulate_u, accumulate_v, epsilon=None):
    """ Compute the decomposition of a (A = U S V^T).
        The matrix a is overwritten in the process.
        The singular values are available in self.sigma.
        The matrices U and V are assembled only if requested, in which case
        they are avalaible in self.u and self.v.
        The argument epsilon is the "precision" to compute the singular
        values with, in the sense of the epsilon machine.
    """
    from scitbx.math import bidiagonal_matrix_kind
    m, n = a.focus()
    u = flex.double()
    u.reshape(flex.grid(0,0))
    v = u

    common_svd_args = libtbx.group_args(accumulate_u=accumulate_u,
                                        accumulate_v=accumulate_v)
    if epsilon is not None:
      common_svd_args.epsilon = epsilon

    if m > self.crossover * n:
      # many more rows than columns
      qr = scitbx.math.householder_qr_decomposition(
        a, may_accumulate_q=accumulate_u)
      r = a.matrix_copy_upper_triangle()
      if accumulate_u:
        qr.accumulate_q_in_place()
        q = a
      bidiag = scitbx.math.householder_bidiagonalisation(r)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      d, f = r.matrix_upper_bidiagonal()
      svd = scitbx.math.svd_decomposition_of_bidiagonal_matrix(
        diagonal=d, off_diagonal=f,
        kind=bidiagonal_matrix_kind.upper_diagonal,
        u=u, v=v,
        **common_svd_args.__dict__)
      svd.compute()
      assert svd.has_converged
      svd.sort()
      if accumulate_u:
        u = q.matrix_multiply(u)
    elif n > self.crossover * m:
      # many more columns than rows
      lq = scitbx.math.householder_lq_decomposition(
        a, may_accumulate_q=accumulate_u)
      l = a.matrix_copy_lower_triangle()
      if accumulate_v:
        lq.accumulate_q_in_place()
        q = a
      bidiag = scitbx.math.householder_bidiagonalisation(l)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      d, f = l.matrix_upper_bidiagonal()
      svd = scitbx.math.svd_decomposition_of_bidiagonal_matrix(
        diagonal=d, off_diagonal=f,
        kind=bidiagonal_matrix_kind.upper_diagonal,
        u=u, v=v,
        **common_svd_args.__dict__)
      svd.compute()
      assert svd.has_converged
      svd.sort()
      if accumulate_v:
        v = q.matrix_transpose().matrix_multiply(v)
    else:
      # balanced number of rows and columns
      bidiag = scitbx.math.householder_bidiagonalisation(a)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      if m >= n:
        d, f = a.matrix_upper_bidiagonal()
        kind=bidiagonal_matrix_kind.upper_diagonal
      else:
        d, f = a.matrix_lower_bidiagonal()
        kind=bidiagonal_matrix_kind.lower_diagonal
      svd = scitbx.math.svd_decomposition_of_bidiagonal_matrix(
        diagonal=d, off_diagonal=f, kind=kind,
        u=u, v=v,
        **common_svd_args.__dict__)
      svd.compute()
      assert svd.has_converged
      svd.sort()

    self.u, self.v, self.sigma = u, v, d

  def reconstruct(self):
    """ Reconstruct the matrix A from its singular values and vectors """
    assert self.u and self.v, "Missing singular vectors"
    return scitbx.math.reconstruct_svd(self.u, self.v, self.sigma)
