import boost.python
boost.python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
import scitbx_array_family_flex_ext as ext

import time
import sys, os

builtin_int = __builtins__["int"]
builtin_long = __builtins__["long"]

class grid_(boost.python.injector, grid):

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "origin:", self.origin()
    print >> f, "last:", self.last()
    print >> f, "focus:", self.focus()
    print >> f, "all:", self.all()
    return self

def export_to(target_module_name):
  export_list = [
    "to_list",
    "select",
    "get_random_seed",
    "random_generator",
    "set_random_seed",
    "random_size_t",
    "random_double",
    "random_permutation",
    "py_object",
    "linear_regression"]
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

def get_random_seed():
  try:
    result = builtin_long(os.getpid() * (2**16)) \
           + builtin_long(time.time() * (2**8))
  except:
    result = time.time()
  return builtin_int(result % (2**31-1))

random_generator = ext.mersenne_twister(seed=get_random_seed())

def set_random_seed(value):
  random_generator.seed(value=value)

random_size_t = random_generator.random_size_t
random_double = random_generator.random_double

def random_permutation(size):
  return sort_permutation(random_double(size=size))

class py_object:

  def __init__(self, accessor, value=None, values=None, value_factory=None):
    assert [value, values, value_factory].count(None) >= 2
    self._accessor = accessor
    if (value_factory is not None):
      self._data = [value_factory() for i in xrange(accessor.size_1d())]
    elif (values is not None):
      assert len(values) == accessor.size_1d()
      self._data = values[:]
    else:
      self._data = [value for i in xrange(accessor.size_1d())]

  def accessor(self):
    return self._accessor

  def data(self):
    return self._data

  def __getitem__(self, index):
    return self._data[self._accessor(index)]

  def __setitem__(self, index, value):
    self._data[self._accessor(index)] = value

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
    b = flex_triple(a.as_double())
    assert tuple(a) == tuple(b)
