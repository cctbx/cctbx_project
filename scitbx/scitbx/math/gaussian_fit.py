import scitbx.math.gaussian
from scitbx import lbfgs
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
import math

def n_less_than(sorted_array, cutoff, eps=1.e-6):
  selection = sorted_array < cutoff + eps
  result = selection.count(0001)
  assert selection[:result].all_eq(0001)
  return result

class minimize:

  def __init__(self, gaussian_fit, target_power,
                     use_sigmas=00000,
                     enforce_positive_b=0001,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=7)):
    adopt_init_args(self, locals())
    assert target_power in [2,4]
    self.n = gaussian_fit.n_terms() * 2
    self.x = flex.double(self.n, 0)
    self.first_target_value = None
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.apply_shifts()
    self.compute_target(compute_gradients=00000)
    self.final_target_value = self.f
    self.final_gaussian_fit = self.gaussian_fit_shifted

  def apply_shifts(self):
    self.gaussian_fit_shifted = self.gaussian_fit.apply_shifts(
      self.x, self.enforce_positive_b)

  def compute_target(self, compute_gradients):
    differences = self.gaussian_fit_shifted.differences()
    self.f = self.gaussian_fit_shifted.target_function(
      self.target_power, self.use_sigmas, differences)
    if (compute_gradients):
      self.g = self.gaussian_fit_shifted.gradients_d_abc(
        self.target_power, self.use_sigmas, differences)
      if (self.enforce_positive_b):
        self.g = self.gaussian_fit.gradients_d_shifts(self.x, self.g)
    else:
      self.g = None

  def __call__(self):
    if (self.first_target_value is None):
      assert self.x.all_eq(0)
      self.gaussian_fit_shifted = self.gaussian_fit
    else:
      self.apply_shifts()
    self.compute_target(compute_gradients=0001)
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.x, self.f, self.g

def make_start_gaussian(null_fit,
                        existing_gaussian,
                        i_x,
                        start_fraction):
  x_sq = null_fit.table_x()[i_x]**2
  y0_table = null_fit.table_y()[0]
  yx_table = null_fit.table_y()[i_x]
  y0_existing = existing_gaussian.at_x_sq(0)
  yx_existing = existing_gaussian.at_x_sq(x_sq)
  n_terms = existing_gaussian.n_terms() + 1
  if (n_terms == 1):
    a = flex.double([y0_table])
    b = flex.double()
    yx_part = yx_table
  else:
    scale_old = 1 - start_fraction
    a = flex.double(existing_gaussian.array_of_a()) * scale_old
    a.append(y0_table - flex.sum(a))
    b = flex.double(existing_gaussian.array_of_b())
    yx_part = yx_table - yx_existing * scale_old
  addl_b = 0
  if (a[-1] != 0):
    r = yx_part / a[-1]
    if (0 < r <= 1):
      addl_b = -math.log(r) / x_sq
  b.append(addl_b)
  if (addl_b != 0):
    assert abs(a[-1] * math.exp(-b[-1] * x_sq) - yx_part) < 1.e-6
  result = scitbx.math.gaussian.fit(
    null_fit.table_x(),
    null_fit.table_y(),
    null_fit.table_sigmas(),
    scitbx.math.gaussian.sum(iter(a), iter(b)))
  if (addl_b != 0):
    assert abs(result.at_x_sq(0) - y0_table) < 1.e-4
  if (n_terms == 1):
    assert abs(result.at_x_sq(x_sq) - yx_table) < 1.e-4
  return result

class find_max_x:

  def __init__(self, gaussian_fit,
                     target_power,
                     minimize_using_sigmas,
                     n_repeats_minimization,
                     enforce_positive_b_mod_n,
                     b_min,
                     max_max_error):
    self.min = None
    self.max_error = None
    table_x = gaussian_fit.table_x()
    table_y = gaussian_fit.table_y()
    sigmas = gaussian_fit.table_sigmas()
    prev_n_points = 0
    good_n_points = 0
    i_x_high = table_x.size() - 1
    while 1:
      if (good_n_points == 0):
        x = (table_x[0] + table_x[i_x_high]) / 2
        n_points = n_less_than(sorted_array=table_x, cutoff=x)
        if (n_points == prev_n_points):
          n_points -= 1
          if (n_points < gaussian_fit.n_terms()*2):
            break
        prev_n_points = n_points
      else:
        n_points = good_n_points + 1
      min_gaussian_fit = scitbx.math.gaussian.fit(
        table_x[:n_points],
        table_y[:n_points],
        sigmas[:n_points],
        gaussian_fit)
      best_minimized = None
      best_max_error = None
      for i in xrange(n_repeats_minimization):
        enforce_positive_b_this_time = (i % enforce_positive_b_mod_n == 0)
        try:
          minimized = minimize(
            gaussian_fit=min_gaussian_fit,
            target_power=target_power,
            use_sigmas=minimize_using_sigmas,
            enforce_positive_b=enforce_positive_b_this_time)
        except RuntimeError, e:
          if (str(e).find("lbfgs error: ") < 0): raise
          if (enforce_positive_b_mod_n == 1): raise
          minimized = None
          max_error = None
          break
        if (min(minimized.final_gaussian_fit.array_of_b()) < b_min):
          break
        min_gaussian_fit = minimized.final_gaussian_fit
        max_error = flex.max(
          minimized.final_gaussian_fit.significant_relative_errors())
        if (    (best_max_error > max_error or best_max_error is None)
            and min(minimized.final_gaussian_fit.array_of_b()) >= b_min):
          best_minimized = minimized
          best_max_error = max_error
      if (best_minimized is not None):
        minimized = best_minimized
        max_error = best_max_error
      if (    max_error > max_max_error
          or min(minimized.final_gaussian_fit.array_of_b()) < b_min):
        if (good_n_points != 0):
          break
        i_x_high = n_points - 1
      else:
        good_n_points = n_points
        good_min = minimized
        good_max_error = max_error
        gaussian_fit = minimized.final_gaussian_fit
        if (good_n_points == table_x.size()):
          break
    if (good_n_points != 0):
      self.min = good_min
      self.max_error = good_max_error

class find_max_x_multi:

  def __init__(self, null_fit,
                     existing_gaussian,
                     target_powers,
                     minimize_using_sigmas,
                     enforce_positive_b_mod_n,
                     b_min,
                     max_max_error,
                     n_start_fractions,
                     n_repeats_minimization,
                     factor_y_x_begin=0.9,
                     factor_y_x_end=0.1):
    i_x_begin = None
    i_x_end = None
    y0 = null_fit.table_y()[0]
    for i,target_value in null_fit.table_y().items():
      if (i_x_begin is None and target_value < y0 * factor_y_x_begin):
        i_x_begin = i
      if (i_x_end is None and target_value < y0 * factor_y_x_end):
        i_x_end = i
        break
    assert i_x_begin is not None
    assert i_x_end is not None
    n_terms = existing_gaussian.n_terms() + 1
    if (n_terms == 1): n_start_fractions = 2
    best_fit = None
    for i_x in xrange(i_x_begin, i_x_end):
      for i_start_fraction in xrange(0,n_start_fractions):
        gaussian_fit = make_start_gaussian(
          null_fit=null_fit,
          existing_gaussian=existing_gaussian,
          i_x=i_x,
          start_fraction=i_start_fraction/float(n_start_fractions))
        for target_power in target_powers:
          fit = find_max_x(
            gaussian_fit=gaussian_fit,
            target_power=target_power,
            minimize_using_sigmas=minimize_using_sigmas,
            n_repeats_minimization=n_repeats_minimization,
            enforce_positive_b_mod_n=enforce_positive_b_mod_n,
            b_min=b_min,
            max_max_error=max_max_error)
          if (fit.min is not None):
            if (best_fit is None
                or best_fit.min.final_gaussian_fit.table_x().size()
                      < fit.min.final_gaussian_fit.table_x().size()
                or best_fit.min.final_gaussian_fit.table_x().size()
                     == fit.min.final_gaussian_fit.table_x().size()
                  and best_fit.max_error > fit.max_error):
              best_fit = fit
    if (best_fit is None):
      self.min = None
      self.max_error = None
    else:
      self.min = best_fit.min
      self.max_error = best_fit.max_error
