from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex # import dependency
import scitbx.array_family.shared # import dependency

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_math_ext")
from scitbx_math_ext import *

import sys

@bp.inject_into(ext.basic_statistics)
class _():

  def show(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print(prefix+"n:", self.n, file=f)
    if (self.n > 0):
      print(prefix+"min:", self.min, file=f)
      print(prefix+"max:", self.max, file=f)
      print(prefix+"max_absolute:", self.max_absolute, file=f)
      print(prefix+"sum:", self.sum, file=f)
      print(prefix+"mean:", self.mean, file=f)
      print(prefix+"mean_absolute_deviation_from_mean:", \
                     self.mean_absolute_deviation_from_mean, file=f)
      print(prefix+"biased_variance:", self.biased_variance, file=f)
      print(prefix+"biased_standard_deviation:", \
                     self.biased_standard_deviation, file=f)
      if (self.n > 1):
        print(prefix+"bias_corrected_variance:", \
                       self.bias_corrected_variance, file=f)
        print(prefix+"bias_corrected_standard_deviation:", \
                       self.bias_corrected_standard_deviation, file=f)
        print(prefix+"skew:", self.skew, file=f)
        print(prefix+"kurtosis:", self.kurtosis, file=f)
        print(prefix+"kurtosis_excess:", self.kurtosis_excess, file=f)

@bp.inject_into(ext.line_search_more_thuente_1994)
class _():

  def show_status(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print(prefix+"xtol:", self.xtol, file=f)
    print(prefix+"ftol:", self.ftol, file=f)
    print(prefix+"gtol:", self.gtol, file=f)
    print(prefix+"stpmin:", self.stpmin, file=f)
    print(prefix+"stpmax:", self.stpmax, file=f)
    print(prefix+"maxfev:", self.maxfev, file=f)
    print(prefix+"info_code:", self.info_code, file=f)
    print(prefix+"info_meaning:", self.info_meaning, file=f)
    print(prefix+"stp:", self.stp, file=f)
    print(prefix+"nfev:", self.nfev, file=f)

@bp.inject_into(ext.unimodular_generator)
class _():

  def all(self):
    while (not self.at_end()): yield next(self)
