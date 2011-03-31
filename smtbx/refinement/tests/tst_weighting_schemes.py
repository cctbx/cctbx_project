from __future__ import division
from cctbx.array_family import flex
from libtbx.test_utils import show_diff
import libtbx.utils
import scitbx.random
from smtbx.refinement import least_squares
import smtbx.development

import random
import math

def exercise_optimise_shelxl_weights():
  def calc_goof(fo2, fc, w, k, n_params):
    fc2 = fc.as_intensity_array()
    w = w(fo2.data(), fo2.sigmas(), fc2.data(), k)
    return math.sqrt(flex.sum(
      w * flex.pow2(fo2.data() - k*fc2.data()))/(fo2.size() - n_params))
  xs = smtbx.development.sucrose()
  k = 0.05 + 10 * flex.random_double()
  fc = xs.structure_factors(anomalous_flag=False, d_min=0.7).f_calc()
  fo = fc.as_amplitude_array()
  fo = fo.customized_copy(data=fo.data()*math.sqrt(k))
  fo = fo.customized_copy(sigmas=0.03*fo.data())
  sigmas = fo.sigmas()
  for i in range(fo.size()):
    fo.data()[i] += 2 * scitbx.random.variate(
      scitbx.random.normal_distribution(sigma=sigmas[i]))() \
      + 0.5*random.random()
  fo2 = fo.as_intensity_array()
  fc2 = fc.as_intensity_array()
  w = least_squares.mainstream_shelx_weighting(a=0.1)
  s = calc_goof(fo2, fc, w, k, xs.n_parameters())
  w2 = w.optimise_parameters(fo2, fc2, k, xs.n_parameters())
  s2 = calc_goof(fo2, fc, w2, k, xs.n_parameters())
  # sort data and setup binning by fc/fc_max
  fc_sq = fc.as_intensity_array()
  fc_sq_over_fc_sq_max = fc_sq.data()/flex.max(fc_sq.data())
  permutation = flex.sort_permutation(fc_sq_over_fc_sq_max)
  fc_sq_over_fc_sq_max = fc_sq.customized_copy(
    data=fc_sq_over_fc_sq_max).select(permutation)
  fc_sq = fc_sq.select(permutation)
  fo_sq = fo2.select(permutation)
  n_bins = 10
  bin_max = 0
  bin_limits = flex.size_t(1, 0)
  bin_count = flex.size_t()
  for i in range(n_bins):
    bin_limits.append(int(math.ceil((i+1) * fc_sq.size()/n_bins)))
    bin_count.append(bin_limits[i+1] - bin_limits[i])
  goofs_w = flex.double()
  goofs_w2 = flex.double()
  for i_bin in range(n_bins):
    sel = flex.size_t_range(bin_limits[i_bin], bin_limits[i_bin+1])
    goofs_w2.append(calc_goof(fo_sq.select(sel),
                              fc_sq.select(sel),
                              w2, k, xs.n_parameters()))
    goofs_w.append(calc_goof(fo_sq.select(sel),
                              fc_sq.select(sel),
                              w, k, xs.n_parameters()))
  a = flex.mean_and_variance(goofs_w).unweighted_sample_variance()
  b = flex.mean_and_variance(goofs_w2).unweighted_sample_variance()
  assert a > b or abs(1-s) > abs(1-s2)
  assert a > b # flat analysis of variance
  assert abs(1-s) > abs(1-s2) # GooF close to 1

def exercise_weighting_schemes():
  unit_weighting = least_squares.unit_weighting()
  assert unit_weighting.type() == "unit"
  assert str(unit_weighting) == "w=1"
  shelx_weighting = least_squares.mainstream_shelx_weighting(0.1234, 0.5678)
  assert shelx_weighting.type() == "calc"
  assert not show_diff(
    str(shelx_weighting),
    "w=1/[\s^2^(Fo^2^)+(0.1234P)^2^+0.5678P] where P=(Fo^2^+2Fc^2^)/3")

def run(args):
  if "--fix_random_seeds" in args:
    random.seed(1)
    flex.set_random_seed(1)
    scitbx.random.set_random_seed(1)
  libtbx.utils.show_times_at_exit()
  for i in range(10):
    exercise_optimise_shelxl_weights()
  exercise_weighting_schemes()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
