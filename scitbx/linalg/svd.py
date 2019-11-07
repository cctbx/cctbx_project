from __future__ import absolute_import, division, print_function
import libtbx
import scitbx.linalg
from scitbx.array_family import flex

class real_bidiagonal(scitbx.linalg.svd_decomposition_of_bidiagonal_matrix):

  def compute(self):
    super(real_bidiagonal, self).compute()
    assert self.has_converged
    self.sort()


class real(object):
  """ SVD decomposition of a real matrix.

      TODO: fix crashes for n x 1 or 1 x n matrices

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
    from scitbx.linalg import bidiagonal_matrix_kind
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
      qr = scitbx.linalg.householder_qr_decomposition(
        a, may_accumulate_q=accumulate_u)
      r = a.matrix_copy_upper_triangle()
      if accumulate_u:
        qr.accumulate_q_in_place()
        q = a
      bidiag = scitbx.linalg.householder_bidiagonalisation(r)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      d, f = r.matrix_upper_bidiagonal()
      svd = real_bidiagonal(diagonal=d, off_diagonal=f,
                            kind=bidiagonal_matrix_kind.upper_diagonal,
                            u=u, v=v,
                            **common_svd_args.__dict__)
      svd.compute()
      if accumulate_u:
        u = q.matrix_multiply(u)
    elif n > self.crossover * m:
      # many more columns than rows
      lq = scitbx.linalg.householder_lq_decomposition(
        a, may_accumulate_q=accumulate_u)
      l = a.matrix_copy_lower_triangle()
      if accumulate_v:
        lq.accumulate_q_in_place()
        q = a
      bidiag = scitbx.linalg.householder_bidiagonalisation(l)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      d, f = l.matrix_upper_bidiagonal()
      svd = real_bidiagonal(diagonal=d, off_diagonal=f,
                            kind=bidiagonal_matrix_kind.upper_diagonal,
                            u=u, v=v,
                            **common_svd_args.__dict__)
      svd.compute()
      if accumulate_v:
        v = q.matrix_transpose().matrix_multiply(v)
    else:
      # balanced number of rows and columns
      bidiag = scitbx.linalg.householder_bidiagonalisation(a)
      if accumulate_u: u = bidiag.u()
      if accumulate_v: v = bidiag.v()
      if m >= n:
        d, f = a.matrix_upper_bidiagonal()
        kind=bidiagonal_matrix_kind.upper_diagonal
      else:
        d, f = a.matrix_lower_bidiagonal()
        kind=bidiagonal_matrix_kind.lower_diagonal
      svd = real_bidiagonal(diagonal=d, off_diagonal=f,
                            kind=kind,
                            u=u, v=v,
                            **common_svd_args.__dict__)
      svd.compute()

    self.u, self.v, self.sigma = u, v, d
    self._svd = svd

  def reconstruct(self):
    """ Reconstruct the matrix A from its singular values and vectors """
    assert self.u and self.v, "Missing singular vectors"
    return scitbx.linalg.reconstruct_svd(self.u, self.v, self.sigma)

  def numerical_rank(self, delta):
    return self._svd.numerical_rank(delta)

def inverse_via_svd(a,is_singular_threshold=1e-16):
  svd_obj = real(a.deep_copy(),True,True,1e-16)
  u = svd_obj.u
  v = svd_obj.v
  sigma = svd_obj.sigma
  selection = flex.bool( sigma < flex.max( sigma )*is_singular_threshold  ).iselection()
  inv_sigma = sigma.deep_copy()
  inv_sigma.set_selected(selection,1.0)
  inv_sigma =  1.0/inv_sigma
  inv_sigma.set_selected(selection,0.0)
  ia = scitbx.linalg.reconstruct_svd(v,u,inv_sigma) #.matrix_transpose(), inv_sigma)
  return(ia, sigma)
