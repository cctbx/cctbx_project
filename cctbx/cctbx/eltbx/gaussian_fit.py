from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
from scitbx import lbfgs
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils.random_utils import random_selection
from scitbx.test_utils import approx_equal

class minimize:

  def __init__(self, diff_gaussian, d_min, n_points=100,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=None):
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
    self.diff_gaussian_shifted = self.diff_gaussian.apply_shifts(self.x)

  def compute_target(self, compute_gradients):
    d_star_max = 1. / self.d_min
    d_step = d_star_max / self.n_points
    dg = self.diff_gaussian_shifted
    self.f = 0
    if (compute_gradients):
      self.g = flex.double(self.n)
      for i in xrange(self.n_points):
        d_star_sq = (i * d_step)**2
        f, g = dg.target_and_target_gradients_at_d_star_sq(d_star_sq)
        self.f += f
        self.g = self.g + g
    else:
      self.g = None
      for i in xrange(self.n_points):
        d_star_sq = (i * d_step)**2
        f = dg.target_at_d_star_sq(d_star_sq)
        self.f += f

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

class difference_gaussian(xray_scattering.gaussian):

  def __init__(self, reference_gaussian,
                     n_terms=None, random_select=00000,
                     gaussian=None):
    assert [n_terms, gaussian].count(None) == 1
    self.reference_gaussian = reference_gaussian
    if (n_terms is None):
      assert gaussian.c() == 0
      xray_scattering.gaussian.__init__(self,
        gaussian.a(), gaussian.b())
    else:
      assert n_terms <= reference_gaussian.n_ab()
      a = flex.double(reference_gaussian.a())
      b = flex.double(reference_gaussian.b())
      if (not random_select):
        permutation = flex.sort_permutation(flex.abs(a), 1)[:n_terms]
      else:
        permutation = random_selection(
          n_candidates=reference_gaussian.n_ab(),
          n_keep=n_terms)
      xray_scattering.gaussian.__init__(self,
        iter(a.select(permutation)),
        iter(b.select(permutation)))

  def new(self, gaussian):
    return difference_gaussian(
      reference_gaussian=self.reference_gaussian,
      gaussian=gaussian)

  def apply_shifts(self, shifts, b_min=-1):
    assert shifts.size() == self.n_ab() * 2
    a = flex.double(self.a()) + shifts[:self.n_ab()]
    b = flex.double(self.b()) + shifts[self.n_ab():]
    selection = b < b_min
    b.set_selected(selection, flex.double(selection.count(0001), b_min))
    return self.new(xray_scattering.gaussian(iter(a), iter(b)))

  def target_term_at_d_star_sq(self, d_star_sq):
    return (self.at_d_star_sq(d_star_sq)
            -self.reference_gaussian.at_d_star_sq(d_star_sq))

  def target_at_d_star_sq(self, d_star_sq):
    return self.target_term_at_d_star_sq(d_star_sq)**2

  def target_and_target_gradients_at_d_star_sq(self, d_star_sq):
    tt = self.target_term_at_d_star_sq(d_star_sq)
    gr = self.gradients_at_d_star_sq(d_star_sq)
    return tt**2, gr*(2*tt)

  def gradients_at_d_star_sq(self, d_star_sq):
    grg = xray_scattering.gaussian.gradients_at_d_star_sq(self, d_star_sq)
    return flex.double(grg.a() + grg.b())

  def finite_diff_gradients_at_d_star_sq(self, d_star_sq, eps=1.e-2):
    gr = flex.double()
    for i in xrange(self.n_ab()):
      t = []
      for seps in (eps, -eps):
        a = list(self.a())
        a[i] += seps
        t.append(
          xray_scattering.gaussian(a, self.b()).at_d_star_sq(d_star_sq))
      gr.append((t[0]-t[1])/(2*eps))
    for i in xrange(self.n_ab()):
      t = []
      for seps in (eps, -eps):
        b = list(self.b())
        b[i] += seps
        t.append(
          xray_scattering.gaussian(self.a(), b).at_d_star_sq(d_star_sq))
      gr.append((t[0]-t[1])/(2*eps))
    return gr

  def finite_diff_target_gradients_at_d_star_sq(self, d_star_sq, eps=1.e-2):
    gr = flex.double()
    for i in xrange(self.n_ab()):
      t = []
      for seps in (eps, -eps):
        a = list(self.a())
        a[i] += seps
        t.append(self.new(xray_scattering.gaussian(
          a, self.b())).target_at_d_star_sq(d_star_sq))
      gr.append((t[0]-t[1])/(2*eps))
    for i in xrange(self.n_ab()):
      t = []
      for seps in (eps, -eps):
        b = list(self.b())
        b[i] += seps
        t.append(self.new(xray_scattering.gaussian(
          self.a(), b)).target_at_d_star_sq(d_star_sq))
      gr.append((t[0]-t[1])/(2*eps))
    return gr

def exercise():
  g5c = xray_scattering.wk1995("C")
  gdiff = difference_gaussian(reference_gaussian=g5c.fetch(), n_terms=4)
  assert approx_equal(gdiff.target_at_d_star_sq(0), 25.1178804546)
  gshifted = gdiff.apply_shifts(flex.double(8,-1))
  assert approx_equal(gshifted.a(),
                      [-5.2410698, 1.657506, 0.49090898, 0.078078985])
  assert approx_equal(gshifted.b(),
                      [-1, 13.780758, 41.086842, -0.223225])
  assert approx_equal(gshifted.gradients_at_d_star_sq(0), [1,1,1,1,0,0,0,0])
  assert approx_equal(gshifted.finite_diff_gradients_at_d_star_sq(0),
                      [1,1,1,1,0,0,0,0], eps=1.e-4)
  for i in xrange(10):
    d_star_sq = (i / 10.)**2
    assert approx_equal(
      gshifted.gradients_at_d_star_sq(d_star_sq),
      gshifted.finite_diff_gradients_at_d_star_sq(d_star_sq), eps=1.e-4)
    assert approx_equal(
       gshifted.target_and_target_gradients_at_d_star_sq(d_star_sq)[1],
       gshifted.finite_diff_target_gradients_at_d_star_sq(d_star_sq),eps=1.e-3)

def multi_trial(reference_gaussian, n_terms, d_min, n_trials):
  best_fit = None
  random_select = 00000
  for i in xrange(n_trials):
    gdiff = difference_gaussian(
      reference_gaussian=reference_gaussian,
      n_terms=n_terms,
      random_select=random_select)
    minimized = minimize(gdiff, d_min=d_min)
    print minimized.final_target_value
    minimized.final_diff_gaussian.show()
    print
    random_select = 0001
    if (   best_fit is None
        or best_fit.final_target_value > minimized.final_target_value):
      best_fit = minimized
  return best_fit

def run():
  exercise()
  for g5c in xray_scattering.wk1995_iterator():
    if (g5c.label() not in ["C", "N", "O", "S"]): continue
    best_fit = multi_trial(
      reference_gaussian=g5c.fetch(),
      n_terms=2,
      d_min=1.5,
      n_trials=10)
    print "Best fit:", g5c.label()
    print best_fit.final_target_value
    best_fit.final_diff_gaussian.show()
    print
    print "=" * 78
    print

if (__name__ == "__main__"):
  run()
