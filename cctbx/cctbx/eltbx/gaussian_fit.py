from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from scitbx import lbfgs
from scitbx.python_utils.misc import adopt_init_args
import random
import math
import sys

def d_star_sq_points(d_min, n_points):
  d_star_max = 1. / d_min
  d_step = d_star_max / n_points
  result = flex.double()
  for i in xrange(n_points):
    result.append((i * d_step)**2)
  return result

class minimize:

  def __init__(self, diff_gaussian, d_star_sq, weights, b_min,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=10)):
    adopt_init_args(self, locals())
    self.n = diff_gaussian.n_ab() * 2
    if (weights is None):
      self.weights = flex.double(d_star_sq.size(), 1)
    self.x = flex.double(self.n, 0)
    self.first_target_value = None
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=00000)
    self.final_target_value = self.f
    self.final_diff_gaussian = self.diff_gaussian_shifted

  def apply_shifts(self):
    self.diff_gaussian_shifted = self.diff_gaussian.apply_shifts(
      self.x, self.b_min)

  def compute_target(self, compute_gradients):
    dg = self.diff_gaussian_shifted
    tt = dg.target_terms_at_points(self.d_star_sq)
    self.f = flex.sum(self.weights * flex.pow2(tt))
    if (compute_gradients):
      self.g = dg.sum_of_gradients_at_points(
        self.d_star_sq, self.weights, tt, 00000)
    else:
      self.g = None

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
      self.diff_gaussian_shifted = self.diff_gaussian
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.x, self.f, self.g

def make_start_approximation(reference_gaussian, n_terms,
                             min_partitioning_factor=0.1,
                             d_max_random_points=10.,
                             d_min_random_points=1/12.):
  assert n_terms <= reference_gaussian.n_ab()
  assert 0 < min_partitioning_factor < 1
  assert d_max_random_points > 0
  assert 0 < d_min_random_points < d_max_random_points
  partitioning = flex.double()
  for i_term in xrange(n_terms):
    partitioning.append(0.1 + random.random())
  min_partitioning = 1. / n_terms * min_partitioning_factor
  while 1:
    partitioning /= flex.sum(partitioning)
    i = flex.min_index(partitioning)
    if (partitioning[i] >= min_partitioning):
      break
    partitioning[i] = min_partitioning
  assert abs(flex.sum(partitioning)-1) < 1.e-6
  f0_reference = reference_gaussian.at_d_star_sq(0)
  a = partitioning * f0_reference
  b = flex.double()
  d_star_min = 1 / d_max_random_points
  d_star_max = 1 / d_min_random_points
  d_star_range = d_star_max - d_star_min
  for i_term in xrange(n_terms):
    d_star_sq = d_star_min + d_star_range * random.random()
    stol_sq = d_star_sq / 4.
    f0_reference = reference_gaussian.at_stol_sq(stol_sq)
    f0_part = partitioning[i_term] * f0_reference
    b.append(-math.log(f0_part / a[i_term]) / stol_sq)
    assert abs(a[i_term] * math.exp(-b[i_term] * stol_sq) - f0_part) < 1.e-6
  return xray_scattering.difference_gaussian(
    reference_gaussian,
    xray_scattering.gaussian(iter(a), iter(b)))

def multi_trial(reference_gaussian, n_terms, d_star_sq, n_trials, verbose=0):
  best_fit = None
  for i in xrange(n_trials):
    gdiff = make_start_approximation(
      reference_gaussian=reference_gaussian,
      n_terms=n_terms)
    assert (gdiff.at_d_star_sq(0) - reference_gaussian.at_d_star_sq(0)) < 1.e-4
    minimized = minimize(
      gdiff, d_star_sq=d_star_sq, weights=d_star_sq, b_min=-1)
    if (0 or verbose):
      print minimized.final_target_value
      minimized.final_diff_gaussian.show()
      print
    if (   best_fit is None
        or best_fit.final_target_value > minimized.final_target_value):
      best_fit = minimized
  return best_fit

def find_d_min(reference_gaussian, n_terms, n_trials_per_d_min,
               start_interval=[10,1/12.],
               min_interval=0.01,
               max_target_value=5.e-4,
               verbose=0):
  interval = start_interval[:]
  fits = []
  for d_min in interval:
    d_star_sq = d_star_sq_points(d_min=d_min, n_points=100)
    fits.append(multi_trial(
      reference_gaussian=reference_gaussian,
      n_terms=n_terms,
      d_star_sq=d_star_sq,
      n_trials=n_trials_per_d_min,
      verbose=verbose))
  while 1:
    if (0 or verbose):
      print interval, [fit.final_target_value for fit in fits]
    if (fits[1].final_target_value < max_target_value):
      break
    d_min = (interval[0] + interval[1]) / 2.
    d_star_sq = d_star_sq_points(d_min=d_min, n_points=100)
    center_fit = multi_trial(
        reference_gaussian=reference_gaussian,
        n_terms=n_terms,
        d_star_sq=d_star_sq,
        n_trials=n_trials_per_d_min,
        verbose=verbose)
    if (center_fit.final_target_value < max_target_value):
      interval[0] = d_min
      fits[0] = center_fit
    else:
      interval[1] = d_min
      fits[1] = center_fit
    if (interval[0] - interval[1] < min_interval):
      break
  if (fits[1].final_target_value < max_target_value):
    best_fit = fits[1]
    d_min = interval[1]
  else:
    assert fits[0].final_target_value < max_target_value
    best_fit = fits[0]
    d_min = interval[0]
  return best_fit, d_min

def show_fit_summary(source, label, gaussian, d_min, q, q_other=None):
  n_terms = str(gaussian.n_ab())
  if (gaussian.c() != 0): n_terms += "+c"
  n_terms += ","
  print "%24s: %s n_terms=%-4s d_min=%.2f, q=%.6g" % (
    source, label, n_terms, d_min, q),
  if (q_other is not None and q_other > q):
    print "Better"
  print

def show_literature_fits(label, d_min, reference_gaussian, n_terms, d_star_sq,
                         q_other=None):
  for lib in [xray_scattering.it1992,
              xray_scattering.two_gaussian_agarwal_isaacs,
              xray_scattering.two_gaussian_agarwal_1978,
              xray_scattering.one_gaussian_agarwal_1978]:
    if (lib == xray_scattering.it1992):
      lib_gaussian = xray_scattering.it1992(label, 1).fetch()
      lib_source = "IT1992"
    elif (lib.table.has_key(label)):
      lib_gaussian = lib.table[label]
      lib_source = lib.source_short
    else:
      lib_gaussian = None
    if (lib_gaussian is not None and lib_gaussian.n_ab() == n_terms):
      diff_gaussian = xray_scattering.difference_gaussian(
        reference_gaussian, lib_gaussian)
      tt = diff_gaussian.target_terms_at_points(d_star_sq)
      q = flex.sum(flex.pow2(tt))
      show_fit_summary(lib_source, label, lib_gaussian, d_min, q, q_other)

def run(args=[], verbose=0):
  for g5c in xray_scattering.wk1995_iterator():
    if (len(args) > 0 and g5c.label() not in args): continue
    for n_terms in [5,4,3,2,1]:
      best_fit, d_min = find_d_min(
        reference_gaussian=g5c.fetch(),
        n_terms=n_terms,
        n_trials_per_d_min=10,
        verbose=verbose)
      show_fit_summary(
        "Best fit", g5c.label(), best_fit.final_diff_gaussian, d_min,
        best_fit.final_target_value)
      show_literature_fits(
        label=g5c.label(),
        reference_gaussian=g5c.fetch(),
        n_terms=n_terms,
        d_min=d_min,
        d_star_sq=best_fit.d_star_sq,
        q_other=best_fit.final_target_value)
      best_fit.final_diff_gaussian.show()
      print
      sys.stdout.flush()
      if (min(best_fit.final_diff_gaussian.b()) <= best_fit.b_min):
        print "Error: Bad fit."
        print
        sys.stdout.flush()

if (__name__ == "__main__"):
  run(sys.argv[1:])
