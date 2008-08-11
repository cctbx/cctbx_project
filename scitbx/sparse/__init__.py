from __future__ import generators
import boost.python
from scitbx_sparse_ext import *
from scitbx.array_family import flex

class _matrix(boost.python.injector, matrix):

  def cols(self):
    for j in xrange(self.n_cols): yield self.col(j)

  def as_dense_matrix(self):
    result = flex.double(flex.grid(self.n_rows, self.n_cols))
    for j,c in enumerate(self.cols()):
      for i,val in c:
        result[i,j] = val
    return result
