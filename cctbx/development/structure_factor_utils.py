from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import scitbx.math
from libtbx.complex_math import abs_arg
from libtbx import adopt_init_args
import math

def check_regression(x, y, label, min_correlation=0, verbose=0):
  xy_regr = flex.linear_regression(x, y)
  xy_corr = flex.linear_correlation(x, y)
  assert xy_regr.is_well_defined()
  if (0 or verbose):
    print(label, "correlation: %.4f slope: %.3f" % (
      xy_corr.coefficient(), xy_regr.slope()))
  assert min_correlation == 0 or xy_corr.coefficient() >= min_correlation

class collector(object):

  def __init__(self, label, min_corr_ampl=0, max_mean_w_phase_error=0,
               verbose=0):
    adopt_init_args(self, locals())
    self.amp1 = flex.double()
    self.amp2 = flex.double()
    self.sum_amp1_minus_amp2_sq = 0
    self.sum_amp1_sq = 0
    self.sum_w_phase_error = 0
    self.sum_w = 0

  def add(self, h, f1, f2):
    a1, p1 = abs_arg(f1, deg=True)
    a2, p2 = abs_arg(f2, deg=True)
    if (0 or self.verbose):
      print(h)
      print(" ", a1, p1)
      print(" ", a2, p2)
      print(" " * 20, "%.2f" % scitbx.math.phase_error(
        phi1=p1, phi2=p2, deg=True))
    self.amp1.append(a1)
    self.amp2.append(a2)
    self.sum_amp1_minus_amp2_sq += (a1 - a2)**2
    self.sum_amp1_sq += a1**2
    self.sum_w_phase_error += (a1 + a2) * scitbx.math.phase_error(
      phi1=p1, phi2=p2, deg=True)
    self.sum_w += (a1 + a2)

  def report(self):
    if (self.sum_amp1_sq):
      r = self.sum_amp1_minus_amp2_sq / self.sum_amp1_sq
      if (0 or self.verbose):
        print(self.label, "R-factor: %.3f" % (r,))
    if (self.sum_w):
      self.mean_w_phase_error = self.sum_w_phase_error / self.sum_w
    check_regression(
      self.amp1, self.amp2, self.label + " ampl", self.min_corr_ampl,
      self.verbose)
    if (0 or self.verbose):
      print(self.label + (" mean weighted phase error: %.2f" % (
        self.mean_w_phase_error,)))
    if (self.max_mean_w_phase_error):
      assert self.mean_w_phase_error <= self.max_mean_w_phase_error

def check_correlation(label, h1, match, f1, f2,
                      min_corr_ampl=0,
                      max_mean_w_phase_error=0,
                      verbose=0):
  coll = collector(
    label, min_corr_ampl, max_mean_w_phase_error, verbose)
  if (match == 0):
    assert f1.size() == f2.size()
    for i,f1_i in enumerate(f1):
      coll.add(h1[i], f1_i, f2[i])
  else:
    for i,j in match.pairs():
      coll.add(h1[i], f1[i], f2[j])
  coll.report()

def check_phase_restrictions(miller_array, epsilon=1.e-10, verbose=0):
  space_group = miller_array.space_group()
  phases = flex.arg(miller_array.data())
  for i,h in enumerate(miller_array.indices()):
    f = miller_array.data()[i]
    if (verbose):
      print(h, f, abs(f), phases[i]*180/math.pi)
    if (abs(f.real) > epsilon and abs(f.imag) > epsilon):
      assert space_group.is_valid_phase(h, phases[i])
