from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from scitbx import lbfgs
from scitbx.python_utils import easy_pickle
from scitbx.python_utils.misc import adopt_init_args
import random
import math
import sys, os

def d_star_sq_points(d_min, n_points):
  d_star_max = 1. / d_min
  d_step = d_star_max / n_points
  result = flex.double()
  for i in xrange(n_points):
    result.append((i * d_step)**2)
  return result

class minimize:

  def __init__(self, diff_gaussian, d_star_sq,
                     weights=None,
                     b_min=-1,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=7)):
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

class d_interval:

  def __init__(self, d_max, d_min):
    assert 0 < d_min < d_max
    self.d_max = d_max
    self.d_min = d_min
    self.d_star_min = 1./d_max
    self.d_star_max = 1./d_min
    self.d_star_range = self.d_star_max - self.d_star_min

  def d_star(self, f_range):
    return self.d_star_min + f_range * self.d_star_range

  def d_star_sq(self, f_range):
    return self.d_star(f_range)**2

  def stol_sq(self, f_range):
    return self.d_star_sq(f_range) / 4

  def d(self, f_range):
    return 1./self.d_star(f_range)

def make_start_approximation(reference_gaussian, n_terms, start_interval,
                             f_trial_point=None,
                             min_partitioning_factor=0.01):
  assert n_terms <= reference_gaussian.n_ab()
  assert f_trial_point is None or n_terms == 1
  assert 0 < min_partitioning_factor < 1
  partitioning = flex.double()
  if (n_terms == 1):
    partitioning.append(1)
  else:
    for i_term in xrange(n_terms):
      partitioning.append(min_partitioning_factor + random.random())
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
  f_range = f_trial_point
  for i_term in xrange(n_terms):
    if (f_trial_point is None):
      f_range = random.random()
    stol_sq = start_interval.stol_sq(f_range)
    f0_reference = reference_gaussian.at_stol_sq(stol_sq)
    f0_part = partitioning[i_term] * f0_reference
    b.append(-math.log(f0_part / a[i_term]) / stol_sq)
    assert abs(a[i_term] * math.exp(-b[i_term] * stol_sq) - f0_part) < 1.e-6
  result = xray_scattering.difference_gaussian(
    reference_gaussian,
    xray_scattering.gaussian(iter(a), iter(b)))
  assert abs(result.at_d_star_sq(0)-reference_gaussian.at_d_star_sq(0)) < 1.e-4
  if (f_trial_point is not None):
    stol_sq = start_interval.stol_sq(f_trial_point)
    assert abs(  result.at_stol_sq(stol_sq)
               - reference_gaussian.at_stol_sq(stol_sq)) < 1.e-4
  return result

class find_d_min:

  def __init__(self, diff_gaussian, max_target_value,
                     n_points,
                     start_interval,
                     dead_end_d_min,
                     min_interval):
    adopt_init_args(self, locals())
    interval = [start_interval.d_max, start_interval.d_min-min_interval]
    good_min = None
    while 1:
      self.d_min = (interval[0] + interval[1]) / 2.
      d_star_sq = d_star_sq_points(d_min=self.d_min, n_points=n_points)
      self.min = minimize(diff_gaussian, d_star_sq=d_star_sq)
      if (    self.min.final_target_value <= max_target_value
          and min(self.min.final_diff_gaussian.b()) > self.min.b_min):
        if (interval[0] - interval[1] < min_interval):
          break
        diff_gaussian = self.min.final_diff_gaussian
        interval[0] = self.d_min
        good_min = self.min
      else:
        if (dead_end_d_min is not None and self.d_min > dead_end_d_min):
          self.d_min = None
          self.min = None
          break
        if (interval[0] - interval[1] < min_interval):
          if (good_min is not None):
            self.d_min = interval[0]
            self.min = good_min
          else:
            self.d_min = None
            self.min = None
          break
        interval[1] = self.d_min

class find_one_term_d_min:

  def __init__(self, reference_gaussian,
                     start_interval,
                     n_trial_points,
                     n_points,
                     max_target_value,
                     min_interval):
    best_fit = None
    for i_trial_point in xrange(n_trial_points+1):
      f_range = i_trial_point / float(n_trial_points)
      d_min = start_interval.d(f_range)
      d_star_sq = d_star_sq_points(d_min=d_min, n_points=n_points)
      diff_gaussian = make_start_approximation(
        reference_gaussian=reference_gaussian,
        n_terms=1,
        start_interval=start_interval,
        f_trial_point=f_range/2)
      fit = find_d_min(
        diff_gaussian=diff_gaussian,
        max_target_value=max_target_value,
        start_interval=start_interval,
        n_points=n_points,
          dead_end_d_min=None,
          min_interval=min_interval)
      if (fit.d_min is not None):
        if (best_fit is None or best_fit.d_min > fit.d_min):
          best_fit = fit
    if (best_fit is None):
      self.d_min = None
      self.min = None
    else:
      self.d_min = fit.d_min
      self.min = fit.min

def show_fit_summary(source, label, gaussian, d_min, q, q_other=None):
  n_terms = str(gaussian.n_ab())
  if (gaussian.c() != 0): n_terms += "+c"
  n_terms += ","
  print "%24s: %s n_terms=%-4s d_min=%.2f, q=%.6g" % (
    source, label, n_terms, d_min, q),
  if (q_other is not None and q_other > q):
    print "Better",
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

def write_plot(f, gaussian, d_min, n_points):
  assert d_min > 0
  d_star_max = 1. / d_min
  for i in xrange(n_points+1):
    d_star = d_star_max * i / n_points
    print >> f, d_star, gaussian.at_d_star_sq(d_star**2)
  print >> f, "&"

def write_plots(plots_dir, label, reference_gaussian, gaussian,
                d_min, n_points):
  label = label.replace("'", "prime") + "_" + str(gaussian.n_ab())
  file_name = os.path.join(plots_dir, label+".xy")
  f = open(file_name, "w")
  write_plot(f, reference_gaussian, d_min, n_points)
  write_plot(f, gaussian, d_min, n_points)
  f.close()

class fit_parameters:

  def __init__(self, n_points=50,
                     max_relative_target_value=1.e-4,
                     find_d_min_min_interval=0.01,
                     n_similar_fits=10,
                     max_d_min_spread_similar_fits=0.02,
                     max_calls_find_d_min_per_term=1000):
    adopt_init_args(self, locals())

def run(args=[], params=fit_parameters(), verbose=0):
  chunk_n = 1
  chunk_i = 0
  if (len(args) > 0 and len(args[0].split(",")) == 2):
    chunk_n, chunk_i = [int(i) for i in args[0].split(",")]
    args = args[1:]
  plots_dir = "plots"
  if (not os.path.isdir(plots_dir)):
    print "***************************************************************"
    print "No plots because target directory does not exist (mkdir %s)." % \
      plots_dir
    print "***************************************************************"
    print
    plots_dir = None
  results = {}
  results["fit_parameters"] = params
  start_interval = d_interval(d_max=15, d_min=1/12.)
  i_chunk = 0
  for g5c in xray_scattering.wk1995_iterator():
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    if (len(args) > 0 and g5c.label() not in args): continue
    reference_gaussian = g5c.fetch()
    f0 = reference_gaussian.at_d_star_sq(0)
    one_term_fit = find_one_term_d_min(
      reference_gaussian=reference_gaussian,
      start_interval=start_interval,
      n_trial_points=3,
      n_points=params.n_points,
      max_target_value=params.max_relative_target_value*f0,
      min_interval=params.find_d_min_min_interval)
    assert one_term_fit.d_min is not None
    show_fit_summary(
      "Best fit", g5c.label(),
      one_term_fit.min.final_diff_gaussian,
      one_term_fit.d_min,
      one_term_fit.min.final_target_value)
    sys.stdout.flush()
    results[g5c.label()] = []
    previous_d_min = None
    for n_terms in [2,3,4,5][:1]:
      fits = []
      fits_d_min = flex.double()
      fits_min_d_min = 15
      n_calls_find_d_min = 0
      while 1:
        diff_gaussian = make_start_approximation(
          reference_gaussian=reference_gaussian,
          n_terms=n_terms,
          start_interval=start_interval)
        n_calls_find_d_min += 1
        fit = find_d_min(
          diff_gaussian=diff_gaussian,
          max_target_value=params.max_relative_target_value*f0,
          start_interval=start_interval,
          n_points=params.n_points,
          dead_end_d_min=fits_min_d_min+params.max_d_min_spread_similar_fits,
          min_interval=params.find_d_min_min_interval)
        if (fit.d_min is not None):
          fits.append(fit)
          fits_d_min.append(fit.d_min)
          fits_min_d_min = min(fits_min_d_min, fit.d_min)
          if (fits_d_min.size() >= params.n_similar_fits):
            perm = flex.sort_permutation(fits_d_min)
            best_d_min = fits_d_min.select(perm[:params.n_similar_fits])
            if (  best_d_min[-1] - best_d_min[0]
                < params.max_d_min_spread_similar_fits):
              fit = fits[perm[0]]
              break
        if (n_calls_find_d_min==params.max_calls_find_d_min_per_term*n_terms):
          if (fits_d_min.size() == 0):
            fit = None
          else:
            perm = flex.sort_permutation(fits_d_min)
            best_d_min = fits_d_min.select(perm[:params.n_similar_fits])
            print "best_d_min:", ", ".join(["%.2f" % d for d in best_d_min])
            fit = fits[perm[0]]
          break
      print "n_calls_find_d_min:", n_calls_find_d_min
      print "successful fits:", fits_d_min.size()
      if (0 or verbose):
        print "fits_d_min:", ", ".join(["%.2f" % d for d in fits_d_min])
      if (fit is None):
        print "No fit: %s n_terms=%d" % (g5c.label(), n_terms)
        print
      else:
        if (previous_d_min is not None and previous_d_min > fit.d_min):
          print "Warning: previous d_min was larger."
        previous_d_min = fit.d_min
        show_fit_summary(
          "Best fit", g5c.label(), fit.min.final_diff_gaussian, fit.d_min,
          fit.min.final_target_value)
        show_literature_fits(
          label=g5c.label(),
          reference_gaussian=reference_gaussian,
          n_terms=n_terms,
          d_min=fit.d_min,
          d_star_sq=fit.min.d_star_sq,
          q_other=fit.min.final_target_value)
        fit.min.final_diff_gaussian.show()
        print
        if (plots_dir):
          write_plots(
            plots_dir=plots_dir,
            label=g5c.label(),
            reference_gaussian=reference_gaussian,
            gaussian=fit.min.final_diff_gaussian,
            d_min=fit.d_min,
            n_points=100)
        g = fit.min.final_diff_gaussian
        results[g5c.label()].append(xray_scattering.fitted_gaussian(
          d_min=fit.d_min, a=g.a(),b=g.b()))
      sys.stdout.flush()
      easy_pickle.dump("fits_%02d.pickle" % chunk_i, results)

if (__name__ == "__main__"):
  from cctbx.eltbx.gaussian_fit import run
  run(sys.argv[1:])
