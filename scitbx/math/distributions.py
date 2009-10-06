from scitbx.math import ext
from ext import normal_distribution
from ext import students_t_distribution
import boost.python

class _injector(boost.python.injector,
                normal_distribution,
                students_t_distribution):

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
