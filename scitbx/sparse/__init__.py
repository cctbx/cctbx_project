from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
import boost_adaptbx.boost.optional # import dependency
from six.moves import range
# Indeed Boost.Python binding needs boost::optional
# through scitbx/random/boost_python/random.h
ext = bp.import_ext("scitbx_sparse_ext")
from scitbx_sparse_ext import *
from scitbx.array_family import flex
import scitbx.random
scitbx.random.variate.register_module(ext)

@bp.inject_into(matrix)
class _():

  def cols(self):
    for j in range(self.n_cols): yield self.col(j)

  def as_dense_matrix(self):
    result = flex.double(flex.grid(self.n_rows, self.n_cols))
    for j,c in enumerate(self.cols()):
      for i,val in c:
        result[i,j] = val
    return result
