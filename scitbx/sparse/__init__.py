import boost.python
import boost.optional # Boost.Python binding needs boost::optional
                      # through scitbx/random/boost_python/random.h
ext = boost.python.import_ext("scitbx_sparse_ext")
from scitbx_sparse_ext import *
from scitbx.array_family import flex
import scitbx.random
scitbx.random.variate.register_module(ext)

class _matrix(boost.python.injector, matrix):

  def cols(self):
    for j in xrange(self.n_cols): yield self.col(j)

  def as_dense_matrix(self):
    result = flex.double(flex.grid(self.n_rows, self.n_cols))
    for j,c in enumerate(self.cols()):
      for i,val in c:
        result[i,j] = val
    return result


class _flex_double(boost.python.injector, flex.double):
  """ Inject method to flex.double to allow
  a = flex.double(...)
  a.reshape(flex.grid(...))
  a += sparse.matrix(...)
  a -= sparse.matrix(...)
  """

  def __iadd__(self, sparse_matrix):
    sparse_matrix._add_to(self)
    return self

  def __isub__(self, sparse_matrix):
    sparse_matrix._substract_from(self)
    return self
