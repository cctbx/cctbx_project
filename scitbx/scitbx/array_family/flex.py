import libtbx.boost_python
ext = libtbx.boost_python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import *

import sys

def to_list(array):
  """Workaround for C++ exception handling bugs
     (list(array) involves C++ exceptions)"""
  result = []
  for i in xrange(array.size()):
    result.append(array[i])
  return result

def linear_regression(x, y, epsilon=1.e-15):
  return ext.linear_regression(x, y, epsilon)

def linear_regression_show_summary(self, f=None):
  if (f is None): f = sys.stdout
  print >> f, "y_intercept:", self.y_intercept()
  print >> f, "slope:", self.slope()

ext.linear_regression.show_summary = linear_regression_show_summary
