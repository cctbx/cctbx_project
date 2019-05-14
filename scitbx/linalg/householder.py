from __future__ import absolute_import, division, print_function
from scitbx.linalg import ext

class householder_qr_decomposition(ext.householder_qr_decomposition):

  def __init__(self, matrix, may_accumulate_q):
    super(householder_qr_decomposition, self).__init__(
      matrix, may_accumulate_q)
    self.r = matrix

class householder_lq_decomposition(ext.householder_lq_decomposition):

  def __init__(self, matrix, may_accumulate_q):
    super(householder_lq_decomposition, self).__init__(
      matrix, may_accumulate_q)
    self.l = matrix

class householder_bidiagonalisation(ext.householder_bidiagonalisation):

  def __init__(self, matrix):
    super(householder_bidiagonalisation, self).__init__(matrix)
    self.bidiagonal = matrix
