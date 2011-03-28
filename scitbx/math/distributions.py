from scitbx.math import ext
from ext import normal_distribution
if (hasattr(ext, "students_t_distribution")):
  from ext import students_t_distribution
else:
  class students_t_distribution(object):
    def __init__(self, *args, **keyword_args):
      raise RuntimeError("Implementation not available in this build.")
import boost.python

class __distribution_mixin(object):

  def mean(self):
    return ext.mean(self)

  def median(self):
    return ext.median(self)

  def mode(self):
    return ext.mode(self)

  def standard_deviation(self):
    return ext.standard_deviation(self)

  def variance(self):
    return ext.variance(self)

  def kurtosis(self):
    return ext.kurtosis(self)

  def skewness(self):
    return ext.skewness(self)

  def pdf(self, x):
    return ext.pdf(self, x)

  def cdf(self, x):
    return ext.cdf(self, x)

  def quantile(self, p):
    return ext.quantile(self, p)

  def quantiles(self, n):
    return ext.quantiles(self, n)

class _(boost.python.injector, normal_distribution, __distribution_mixin):
  pass

class _(boost.python.injector, students_t_distribution, __distribution_mixin):
  pass
