from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from scitbx import lbfgs
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils.random_utils import random_selection

def d_star_sq_points(d_min, n_points):
  d_star_max = 1. / d_min
  d_step = d_star_max / n_points
  result = flex.double()
  for i in xrange(n_points):
    result.append((i * d_step)**2)
  return result

class minimize:

  def __init__(self, diff_gaussian, d_star_sq, b_min,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=10)):
    adopt_init_args(self, locals())
    self.n = diff_gaussian.n_ab() * 2
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
    self.f = flex.sum(flex.pow2(tt))
    if (compute_gradients):
      self.g = dg.sum_of_gradients_at_points(self.d_star_sq, tt)
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

def make_difference_gaussian(reference_gaussian,
                             n_terms=None,
                             random_select=00000):
  assert n_terms <= reference_gaussian.n_ab()
  a = flex.double(reference_gaussian.a())
  b = flex.double(reference_gaussian.b())
  if (not random_select):
    permutation = flex.sort_permutation(flex.abs(a), 1)[:n_terms]
  else:
    permutation = random_selection(
      n_candidates=reference_gaussian.n_ab(),
      n_keep=n_terms)
  return xray_scattering.difference_gaussian(
    reference_gaussian,
    xray_scattering.gaussian(
      iter(a.select(permutation)),
      iter(b.select(permutation))))

def multi_trial(reference_gaussian, n_terms, d_star_sq, n_trials, verbose=0):
  best_fit = None
  random_select = 00000
  for i in xrange(n_trials):
    gdiff = make_difference_gaussian(
      reference_gaussian=reference_gaussian,
      n_terms=n_terms,
      random_select=random_select)
    minimized = minimize(gdiff, d_star_sq=d_star_sq, b_min=-1)
    if (0 or verbose):
      print minimized.final_target_value
      minimized.final_diff_gaussian.show()
      print
    random_select = 0001
    if (   best_fit is None
        or best_fit.final_target_value > minimized.final_target_value):
      best_fit = minimized
  return best_fit

def show_agarwal_fits(label, d_min, reference_gaussian, d_star_sq):
  for lib in [xray_scattering.two_gaussian_agarwal_isaacs,
              xray_scattering.two_gaussian_agarwal_1978,
              xray_scattering.one_gaussian_agarwal_1978]:
    if (lib.table.has_key(label)):
      lib_gaussian = lib.table[label]
      diff_gaussian = xray_scattering.difference_gaussian(
        reference_gaussian, lib_gaussian)
      tt = diff_gaussian.target_terms_at_points(d_star_sq)
      q = flex.sum(flex.pow2(tt))
      print "%33s: %s n_terms=%d, d_min=%.1f, q=%.6g" % (
        lib.source_short, label, lib_gaussian.n_ab(), d_min, q)

def run(verbose=0):
  for g5c in xray_scattering.wk1995_iterator():
    if (g5c.label() not in ["C", "N", "O", "S"]): continue
    for d_min in [6,4,3,2,1.5,1]:
      d_star_sq = d_star_sq_points(d_min=d_min, n_points=100)
      reference_gaussian = g5c.fetch()
      show_agarwal_fits(
        label=g5c.label(),
        d_min=d_min,
        reference_gaussian=reference_gaussian,
        d_star_sq=d_star_sq)
      for n_terms in [4,3,2,1]:
        best_fit = multi_trial(
          reference_gaussian=reference_gaussian,
          n_terms=n_terms,
          d_star_sq=d_star_sq,
          n_trials=10,
          verbose=verbose)
        print "%33s: %s n_terms=%d, d_min=%.1f, q=%.6g" % (
          "Best fit", g5c.label(), n_terms, d_min,
          best_fit.final_target_value)
        if (0 or verbose):
          best_fit.final_diff_gaussian.show()
          print
          print "=" * 78
          print
        assert min(best_fit.final_diff_gaussian.b()) > best_fit.b_min
        reference_gaussian = best_fit.final_diff_gaussian
      if (not verbose): print

if (__name__ == "__main__"):
  run()
