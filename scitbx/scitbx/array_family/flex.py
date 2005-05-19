from __future__ import generators
import scitbx.boost_python.slice

import boost.python
boost.python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *
import scitbx_array_family_flex_ext as ext
from libtbx import introspection

import scitbx.stl.map
import md5
import time
import sys, os

builtin_int = __builtins__["int"]
builtin_long = __builtins__["long"]

def bool_md5(self):
  result = md5.new()
  result.update(self.__getstate__()[1])
  return result
bool.md5 = bool_md5

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
    "min_default",
    "max_default",
    "mean_default",
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

def min_default(values, default):
  if (values.size() == 0): return default
  return min(values)

def max_default(values, default):
  if (values.size() == 0): return default
  return max(values)

def mean_default(values, default):
  if (values.size() == 0): return default
  return mean(values)

def select(sequence, permutation=None, flags=None):
  result = []
  if (permutation is not None):
    assert flags is None
    for i in permutation:
      result.append(sequence[i])
  else:
    assert flags is not None
    for s,f in zip(sequence, flags):
      if (f): result.append(s)
  return result

def get_random_seed():
  try:
    result = builtin_long(os.getpid() * (2**16)) \
           + builtin_long(time.time() * (2**8))
  except KeyboardInterrupt: raise
  except:
    result = time.time()
  return builtin_int(result % (2**31-1))

random_generator = ext.mersenne_twister(seed=get_random_seed())

def set_random_seed(value):
  random_generator.seed(value=value)

random_size_t = random_generator.random_size_t
random_double = random_generator.random_double
random_permutation = random_generator.random_permutation
random_integer = random_generator.random_integer

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

class _linear_regression_core(boost.python.injector,
                              ext.linear_regression_core):

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "is_well_defined:", self.is_well_defined()
    print >> f, "y_intercept:", self.y_intercept()
    print >> f, "slope:", self.slope()

class _linear_correlation(boost.python.injector, ext.linear_correlation):

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "is_well_defined:", self.is_well_defined()
    print >> f, "mean_x:", self.mean_x()
    print >> f, "mean_y:", self.mean_y()
    print >> f, "coefficient:", self.coefficient()

class histogram_slot_info:

  def __init__(self, low_cutoff, high_cutoff, n):
    introspection.adopt_init_args()

class _histogram(boost.python.injector, ext.histogram):

  def slot_infos(self):
    low_cutoff = self.data_min()
    for i,n in enumerate(self.slots()):
      high_cutoff = self.data_min() + self.slot_width() * (i+1)
      yield histogram_slot_info(low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    for info in self.slot_infos():
      print >> f, "%s%.8g - %.8g: %d" % (
        prefix, info.low_cutoff, info.high_cutoff, info.n)

def exercise_triple(flex_triple, flex_order=None, as_double=False):
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
