import boost.python
boost.python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
import scitbx_array_family_flex_ext as ext

import sys

def export_to(target_module_name):
  export_list = ["to_list", "select", "linear_regression"]
  target_module = sys.modules[target_module_name]
  g = globals()
  for attr in export_list:
    setattr(target_module, attr, g[attr])

def to_list(array):
  """Workaround for C++ exception handling bugs
     (list(array) involves C++ exceptions)"""
  result = []
  for i in xrange(array.size()):
    result.append(array[i])
  return result

def select(sequence, permutation):
  result = []
  for i in permutation:
    result.append(sequence[i])
  return result

class linear_regression(ext.linear_regression):

  def __init__(self, x, y, epsilon=1.e-15):
    ext.linear_regression.__init__(self, x, y, epsilon)

class _linear_regression(boost.python.injector, ext.linear_regression):

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "y_intercept:", self.y_intercept()
    print >> f, "slope:", self.slope()

def exercise_triple(flex_triple, flex_order=None, as_double=00000):
  from libtbx.test_utils import approx_equal
  import pickle
  a = flex_triple()
  a = flex_triple(((1,2,3), (2,3,4), (3,4,5)))
  assert a.size() == 3
  assert tuple(a) == ((1,2,3), (2,3,4), (3,4,5))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert tuple(a) == tuple(b)
  if (flex_order is not None):
    assert flex_order(a, b) == 0
  if (as_double):
    assert approx_equal(tuple(a.as_double()), (1,2,3,2,3,4,3,4,5))
    b = flex_triple().from_double(a.as_double())
    assert tuple(a) == tuple(b)
