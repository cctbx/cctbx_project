import boost.python
ext = boost.python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *

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
