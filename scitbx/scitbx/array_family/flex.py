from scitbx.python_utils import misc
ext = misc.import_ext("scitbx_boost.array_family.flex_scitbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

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

def linear_regression_show_summary(self, f=sys.stdout):
  print >> f, "y_intercept:", self.y_intercept()
  print >> f, "slope:", self.slope()

ext.linear_regression.show_summary = linear_regression_show_summary
