import scitbx.array_family.flex # import dependency
import scitbx.array_family.shared # import dependency

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

class _line_search_more_thuente_1994(
        boost.python.injector, ext.line_search_more_thuente_1994):

  def show_status(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix+"xtol:", self.xtol
    print >> f, prefix+"ftol:", self.ftol
    print >> f, prefix+"gtol:", self.gtol
    print >> f, prefix+"stpmin:", self.stpmin
    print >> f, prefix+"stpmax:", self.stpmax
    print >> f, prefix+"maxfev:", self.maxfev
    print >> f, prefix+"info_code:", self.info_code
    print >> f, prefix+"info_meaning:", self.info_meaning
    print >> f, prefix+"stp:", self.stp
    print >> f, prefix+"nfev:", self.nfev

class _unimodular_generator(boost.python.injector, ext.unimodular_generator):

  def all(self):
    while (not self.at_end()): yield self.next()
