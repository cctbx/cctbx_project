import scitbx.math.gaussian
from scitbx.math import golay_24_12_generator
from scitbx import lbfgs
from scitbx import lbfgsb
from cctbx.array_family import flex
from scitbx.python_utils.math_utils import ifloor
from scitbx.python_utils.misc import adopt_init_args
from scitbx.python_utils import dicts
import math
import sys

def n_less_than(sorted_array, cutoff, eps=1.e-6):
  selection = sorted_array < cutoff + eps
  result = selection.count(0001)
  assert selection[:result].all_eq(0001)
  return result

class LargeNegativeB(RuntimeError): pass

class minimize_mixin:

  def apply_shifts(self):
    self.gaussian_fit_shifted = self.gaussian_fit.apply_shifts(
      self.x, self.shift_sqrt_b)
    min_b = min(self.gaussian_fit_shifted.array_of_b())
    if (min_b < self.hard_b_min):
      raise LargeNegativeB("gaussian_fit error: large negative b")

  def compute_fg(self):
    self.apply_shifts()
    differences = self.gaussian_fit_shifted.differences()
    self.f = self.gaussian_fit_shifted.target_function(
      self.target_power, self.use_sigmas, differences)
    self.g = self.gaussian_fit_shifted.gradients_d_abc(
      self.target_power, self.use_sigmas, differences)
    if (self.shift_sqrt_b):
      self.g = self.gaussian_fit.gradients_d_shifts(self.x, self.g)

  def finalize(self):
    self.compute_fg()
    self.final_target_value = self.f
    self.final_gaussian_fit = self.gaussian_fit_shifted
    self.max_error = flex.max(
      self.final_gaussian_fit.significant_relative_errors())

  def show_minimization_parameters(self, f=None):
    if (f is None): f = sys.stdout
    s = str(self)
    if (hasattr(self, "shift_sqrt_b_mod_n")):
      s += " shift_sqrt_b_mod_n:%d" % self.shift_sqrt_b_mod_n
      s += " i_repeat:%d" % self.i_repeat
    print >> f, s

class minimize_lbfgs(minimize_mixin):

  def __init__(self, gaussian_fit, target_power,
                     use_sigmas,
                     shift_sqrt_b,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=20),
                     hard_b_min=-1):
    adopt_init_args(self, locals())
    assert target_power in [2,4]
    self.n = gaussian_fit.n_terms() * 2
    self.x = flex.double(self.n, 0)
    self.first_target_value = None
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.finalize()

  def __str__(self):
    return "lbfgs:m=%d:tp=%d:us=%d:sq=%d" % (
      self.lbfgs_core_params.m,
      self.target_power,
      int(self.use_sigmas),
      int(self.shift_sqrt_b))

  def __call__(self):
    self.compute_fg()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.x, self.f, self.g

class minimize_lbfgsb(minimize_mixin):

  def __init__(self, gaussian_fit, target_power,
                     use_sigmas,
                     shift_sqrt_b,
                     apply_lower_bounds_on_b,
                     lbfgsb_m=20,
                     hard_b_min=-1):
    adopt_init_args(self, locals())
    assert target_power in [2,4]
    self.n = gaussian_fit.n_terms() * 2
    self.run()
    self.finalize()

  def __str__(self):
    return "lbfgsb:m=%d:tp=%d:us=%d:sq=%d:lb=%d" % (
      self.lbfgsb_m,
      self.target_power,
      int(self.use_sigmas),
      int(self.shift_sqrt_b),
      int(self.apply_lower_bounds_on_b))

  def run(self):
    l = flex.double(self.n, 0)
    u = flex.double(self.n, 0)
    nbd = flex.int(self.n, 0)
    if (self.apply_lower_bounds_on_b):
      bound_flags = self.gaussian_fit.bound_flags(00000, 0001)
      nbd.set_selected(bound_flags, flex.int(bound_flags.count(0001), 1))
    self.minimizer = lbfgsb.minimizer(
      n=self.n,
      m=self.lbfgsb_m,
      l=l,
      u=u,
      nbd=nbd)
    self.x = flex.double(self.n, 0)
    self.f = 0
    self.g = flex.double(self.n, 0)
    while 0001:
      task = self.minimizer.process(self.x, self.f, self.g)
      if (task[:2] == "FG"):
        self.compute_fg()
      elif (task[:5] != "NEW_X"):
        break

def minimize_multi_lbfgs(start_fit,
                         target_powers,
                         minimize_using_sigmas,
                         shift_sqrt_b_mod_n,
                         b_min,
                         n_repeats_minimization):
  best_min = None
  for target_power in target_powers:
    min_gaussian_fit = start_fit
    for i in xrange(n_repeats_minimization):
      if (shift_sqrt_b_mod_n > 0):
        shift_sqrt_b_this_time = (i % shift_sqrt_b_mod_n == 0)
      else:
        shift_sqrt_b_this_time = 0
      try:
        minimized = minimize_lbfgs(
          gaussian_fit=min_gaussian_fit,
          target_power=target_power,
          use_sigmas=minimize_using_sigmas,
          shift_sqrt_b=shift_sqrt_b_this_time)
        minimized.shift_sqrt_b_mod_n = shift_sqrt_b_mod_n
        minimized.i_repeat = i
      except LargeNegativeB:
        minimized = None
        break
      except RuntimeError, e:
        if (str(e).find("lbfgs error: ") < 0): raise
        print e
        print "Aborting this minimization."
        print
        sys.stdout.flush()
        minimized = None
        break
      if (min(minimized.final_gaussian_fit.array_of_b()) < b_min):
        minimized = None
        break
      min_gaussian_fit = minimized.final_gaussian_fit
      if (best_min is None or best_min.max_error > minimized.max_error):
        best_min = minimized
  return best_min

def minimize_multi_lbfgsb(start_fit,
                          target_powers,
                          minimize_using_sigmas,
                          shift_sqrt_b_mod_n,
                          b_min,
                          n_repeats_minimization):
  best_min = None
  for target_power in target_powers:
    for apply_lower_bounds_on_b in [00000, 0001]:
      min_gaussian_fit = start_fit
      for i in xrange(n_repeats_minimization):
        if (shift_sqrt_b_mod_n > 0):
          shift_sqrt_b_this_time = (i % shift_sqrt_b_mod_n == 0)
        else:
          shift_sqrt_b_this_time = 0
        try:
          minimized = minimize_lbfgsb(
            gaussian_fit=min_gaussian_fit,
            target_power=target_power,
            use_sigmas=minimize_using_sigmas,
            shift_sqrt_b=shift_sqrt_b_this_time,
            apply_lower_bounds_on_b=apply_lower_bounds_on_b)
          minimized.shift_sqrt_b_mod_n = shift_sqrt_b_mod_n
          minimized.i_repeat = i
        except LargeNegativeB:
          minimized = None
          break
        if (min(minimized.final_gaussian_fit.array_of_b()) < b_min):
          minimized = None
          break
        min_gaussian_fit = minimized.final_gaussian_fit
        if (best_min is None or best_min.max_error > minimized.max_error):
          best_min = minimized
  return best_min

minimize_multi_histogram = dicts.with_default_value(0)

def show_minimize_multi_histogram(f=None, reset=0001):
  global minimize_multi_histogram
  minimizer_types = minimize_multi_histogram.keys()
  counts = flex.double(minimize_multi_histogram.values())
  perm = flex.sort_permutation(counts, 0001)
  minimizer_types = flex.select(minimizer_types, perm)
  counts = counts.select(perm)
  n_total = flex.sum(counts)
  for m,c in zip(minimizer_types, counts):
    print >> f, "%-32s  %5.3f" % (m, c/n_total)
  print >> f
  if (reset):
    minimize_multi_histogram = dicts.with_default_value(0)

def minimize_multi(start_fit,
                   target_powers,
                   minimize_using_sigmas,
                   shift_sqrt_b_mod_n,
                   b_min,
                   n_repeats_minimization):
  best_min_list = []
  for current_shift_sqrt_b_mod_n in shift_sqrt_b_mod_n:
    for minimize_multi_type in [minimize_multi_lbfgs, minimize_multi_lbfgsb]:
      best_min_list.append(minimize_multi_type(
        start_fit=start_fit,
        target_powers=target_powers,
        minimize_using_sigmas=minimize_using_sigmas,
        shift_sqrt_b_mod_n=current_shift_sqrt_b_mod_n,
        b_min=b_min,
        n_repeats_minimization=n_repeats_minimization))
  best_best_min = None
  for best_min in best_min_list:
    if (best_min is None): continue
    if (best_best_min is None or best_best_min.max_error > best_min.max_error):
      best_best_min = best_min
  minimize_multi_histogram[str(best_best_min)] += 1
  return best_best_min

def find_max_x(gaussian_fit,
               target_powers,
               minimize_using_sigmas,
               n_repeats_minimization,
               shift_sqrt_b_mod_n,
               b_min,
               max_max_error):
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
    start_fit = scitbx.math.gaussian.fit(
      table_x[:n_points],
      table_y[:n_points],
      sigmas[:n_points],
      gaussian_fit)
    best_min = minimize_multi(
      start_fit=start_fit,
      target_powers=target_powers,
      minimize_using_sigmas=minimize_using_sigmas,
      shift_sqrt_b_mod_n=shift_sqrt_b_mod_n,
      b_min=b_min,
      n_repeats_minimization=n_repeats_minimization)
    if (best_min is None or best_min.max_error > max_max_error):
      if (good_n_points != 0):
        break
      i_x_high = n_points - 1
    else:
      good_n_points = n_points
      good_min = best_min
      gaussian_fit = best_min.final_gaussian_fit
      if (good_n_points == table_x.size()):
        break
  if (good_n_points != 0):
    return good_min
  return None

def make_start_gaussian(null_fit,
                        existing_gaussian,
                        i_split,
                        i_x,
                        start_fraction,
                        b_range=1.e-3):
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
    b = flex.double(existing_gaussian.array_of_b())
    b_max = flex.max(flex.abs(b))
    b_min = b_max * b_range
    sel = b < b_min
    b.set_selected(sel, flex.double(sel.count(0001), b_min))
    if (i_split < 0):
      a = flex.double(existing_gaussian.array_of_a()) * scale_old
      a.append(y0_table - flex.sum(a))
      yx_part = yx_table - yx_existing * scale_old
    else:
      t_split = scitbx.math.gaussian.term(
        existing_gaussian.array_of_a()[i_split],
        existing_gaussian.array_of_b()[i_split])
      a = flex.double(existing_gaussian.array_of_a())
      a.append(a[i_split] * start_fraction)
      a[i_split] *= scale_old
      yx_part = t_split.at_x_sq(x_sq) * start_fraction
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
  if (addl_b != 0 and i_split < 0):
    assert abs(result.at_x_sq(0) - y0_table) < 1.e-4
  if (n_terms == 1):
    assert abs(result.at_x_sq(x_sq) - yx_table) < 1.e-4
  return result

def find_max_x_multi(null_fit,
                     existing_gaussian,
                     target_powers,
                     minimize_using_sigmas,
                     shift_sqrt_b_mod_n,
                     b_min,
                     max_max_error,
                     n_start_fractions,
                     n_repeats_minimization,
                     factor_y_x_begin=0.9,
                     factor_y_x_end=0.1,
                     factor_x_step=2.):
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
  i_x_step = max(1, ifloor((i_x_end-i_x_begin) / (factor_x_step*n_terms)))
  if (n_terms == 1): n_start_fractions = 2
  best_min = None
  for i_x in xrange(i_x_begin, i_x_end, i_x_step):
    for i_split in xrange(-1, existing_gaussian.n_terms()):
      for i_start_fraction in xrange(0,n_start_fractions):
        gaussian_fit = make_start_gaussian(
          null_fit=null_fit,
          existing_gaussian=existing_gaussian,
          i_split=i_split,
          i_x=i_x,
          start_fraction=i_start_fraction/float(n_start_fractions))
        for target_power in target_powers:
          good_min = find_max_x(
            gaussian_fit=gaussian_fit,
            target_powers=[target_power],
            minimize_using_sigmas=minimize_using_sigmas,
            n_repeats_minimization=n_repeats_minimization,
            shift_sqrt_b_mod_n=shift_sqrt_b_mod_n,
            b_min=b_min,
            max_max_error=max_max_error)
          if (good_min is not None):
            if (best_min is None
                or best_min.final_gaussian_fit.table_x().size()
                 < good_min.final_gaussian_fit.table_x().size()
                or best_min.final_gaussian_fit.table_x().size()
                == good_min.final_gaussian_fit.table_x().size()
                  and best_min.max_error > good_min.max_error):
              best_min = good_min
  return best_min

def make_golay_based_start_gaussian(null_fit, code):
  assert len(code) == 24
  a_starts = [1,4,16,32]
  b_starts = [1,4,16,32]
  a = flex.double()
  b = flex.double()
  for i_term in xrange(6):
    i_bits = i_term * 2
    bits_a = code[i_bits], code[i_bits+12]
    bits_b = code[i_bits+1], code[i_bits+12+1]
    a.append(a_starts[bits_a[0]*2+bits_a[1]])
    b.append(b_starts[bits_b[0]*2+bits_b[1]])
  a = a * null_fit.table_y()[0] / flex.sum(a)
  return scitbx.math.gaussian.fit(
    null_fit.table_x(),
    null_fit.table_y(),
    null_fit.table_sigmas(),
    scitbx.math.gaussian.sum(iter(a), iter(b)))

def fit_with_golay_starts(label,
                          null_fit,
                          null_fit_more,
                          n_terms,
                          target_powers,
                          minimize_using_sigmas,
                          shift_sqrt_b_mod_n,
                          b_min,
                          n_repeats_minimization,
                          negligible_max_error=0.001,
                          print_to=None):
  assert n_terms == 6
  if (label is not None and print_to is None):
    print_to = sys.stdout
  good_min = None
  for golay_code in golay_24_12_generator():
    start_fit = make_golay_based_start_gaussian(
      null_fit=null_fit,
      code=golay_code)
    best_min = minimize_multi(
      start_fit=start_fit,
      target_powers=target_powers,
      minimize_using_sigmas=minimize_using_sigmas,
      shift_sqrt_b_mod_n=shift_sqrt_b_mod_n,
      b_min=b_min,
      n_repeats_minimization=n_repeats_minimization)
    if (best_min is not None):
      if (good_min is None or good_min.max_error > best_min.max_error):
        good_min = best_min
        fit_more = scitbx.math.gaussian.fit(
          null_fit_more.table_x(),
          null_fit_more.table_y(),
          null_fit_more.table_sigmas(),
          good_min.final_gaussian_fit)
        if (print_to is not None):
          print >> print_to, label, "max_error fitted=%.4f, more=%.4f" % (
            good_min.max_error,
            flex.max(fit_more.significant_relative_errors()))
          fit_more.show(f=print_to)
          fit_more.show_table(f=print_to)
          print >> print_to
          print_to.flush()
        if (good_min.max_error <= negligible_max_error):
          break
  if (print_to is not None):
    if (good_min is None):
      print >> print_to, "Final: %s: No successful minimization." % label
    else:
      print >> print_to, "Final:", label, "max_error fitted=%.4f, more=%.4f" % (
        good_min.max_error,
        flex.max(fit_more.significant_relative_errors()))
    good_min.show_minimization_parameters(f=print_to)
    print >> print_to
    show_minimize_multi_histogram(f=print_to)
  return good_min
