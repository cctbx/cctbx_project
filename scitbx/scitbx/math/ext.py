import boost.python
ext = boost.python.import_ext("scitbx_math_ext")
from scitbx_math_ext import *

import sys

class _basic_statistics(boost.python.injector, ext.basic_statistics):

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix+"n:", self.n
    if (self.n > 0):
      print >> f, prefix+"min:", self.min
      print >> f, prefix+"max:", self.max
      print >> f, prefix+"max_absolute:", self.max_absolute
      print >> f, prefix+"sum:", self.sum
      print >> f, prefix+"mean:", self.mean
      print >> f, prefix+"mean_absolute_deviation_from_mean:", \
                     self.mean_absolute_deviation_from_mean
      print >> f, prefix+"biased_variance:", self.biased_variance
      print >> f, prefix+"biased_standard_deviation:", \
                     self.biased_standard_deviation
      if (self.n > 1):
        print >> f, prefix+"bias_corrected_variance:", \
                       self.bias_corrected_variance
        print >> f, prefix+"bias_corrected_standard_deviation:", \
                       self.bias_corrected_standard_deviation
        print >> f, prefix+"skew:", self.skew
        print >> f, prefix+"kurtosis:", self.kurtosis
        print >> f, prefix+"kurtosis_excess:", self.kurtosis_excess
