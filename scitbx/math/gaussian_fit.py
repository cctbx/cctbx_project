import scitbx.math.gaussian
from scitbx.math import golay_24_12_generator
from scitbx import lbfgs
from scitbx import lbfgsb
from scitbx.examples import immoptibox_ports
from scitbx.array_family import flex
from libtbx.math_utils import ifloor
from libtbx import adopt_init_args
from libtbx import easy_pickle
import time
import math
import sys

minimize_multi_histogram = {"None": 0}

def n_less_than(sorted_array, cutoff, eps=1.e-6):
  selection = sorted_array < cutoff + eps
  result = selection.count(True)
  assert selection[:result].all_eq(True)
  return result

class LargeNegativeB(RuntimeError): pass

class minimize_mixin(object):

  def apply_shifts(self):
    self.gaussian_fit_shifted = self.gaussian_fit.apply_shifts(
      self.x, self.shift_sqrt_b!=0)
    min_b = min(self.gaussian_fit_shifted.array_of_b())
    if (min_b < self.hard_b_min):
      raise LargeNegativeB("gaussian_fit error: large negative b")

  def finalize(self):
    self.final_gaussian_fit = self.gaussian_fit_shifted
    self.max_x = self.final_gaussian_fit.table_x()[-1]
    self.max_error = flex.max(
      self.final_gaussian_fit.significant_relative_errors())

  def is_better_than(self, other):
    if (other is None): return True
    if (self.max_x > other.max_x): return True
    if (self.max_x < other.max_x): return False
    if (self.max_error < other.max_error): return True
    return False

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "n_terms:", self.final_gaussian_fit.n_terms(),
    print >> f, "max_x: %.2f" % self.final_gaussian_fit.table_x()[-1],
    print >> f, "max_error: %.4f" % self.max_error
    f.flush()
    return self

  def show_minimization_parameters(self, f=None):
    if (f is None): f = sys.stdout
    s = str(self)
    if (hasattr(self, "shift_sqrt_b_mod_n")):
      s += " shift_sqrt_b_mod_n:%d" % self.shift_sqrt_b_mod_n
      s += " i_repeat:%d" % self.i_repeat
    print >> f, s
    f.flush()
    return self

class minimize_lbfgs_mixin(minimize_mixin):

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
    minimize_mixin.finalize(self)

class minimize_lbfgs(minimize_lbfgs_mixin):

  def __init__(self, gaussian_fit, target_power,
                     use_sigmas,
                     shift_sqrt_b,
                     lbfgs_termination_params=None,
                     lbfgs_core_params=lbfgs.core_parameters(m=20),
                     hard_b_min=-1):
    adopt_init_args(self, locals())
    minimize_multi_histogram.setdefault(str(self), 0)
    assert target_power in [2,4]
    self.x = flex.double(gaussian_fit.n_terms() * 2, 0)
    self.first_target_value = None
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      core_params=lbfgs_core_params)
    self.finalize()

  def __str__(self):
    return "lbfgs:np=%d:m=%d:tp=%d:us=%d:sq=%d" % (
      self.gaussian_fit.n_parameters(),
      self.lbfgs_core_params.m,
      self.target_power,
      int(self.use_sigmas),
      int(self.shift_sqrt_b))

  def compute_functional_and_gradients(self):
    self.compute_fg()
    if (self.first_target_value is None):
      self.first_target_value = self.f
    return self.f, self.g

class minimize_lbfgsb(minimize_lbfgs_mixin):

  def __init__(self, gaussian_fit, target_power,
                     use_sigmas,
                     shift_sqrt_b,
                     apply_lower_bounds_on_b,
                     lbfgsb_m=20,
                     hard_b_min=-1,
                     iprint=-1,
                     use_fortran_library=False):
    adopt_init_args(self, locals())
    minimize_multi_histogram.setdefault(str(self), 0)
    assert target_power in [2,4]
    self.n = gaussian_fit.n_terms() * 2
    try:
      self.run()
    except LargeNegativeB:
      raise
    except:
      easy_pickle.dump(
        "lbfgsb_exception_%.0f.pickle" % time.time(),
        gaussian_fit)
      raise
    self.finalize()

  def __str__(self):
    return "lbfgsb:np=%d:m=%d:tp=%d:us=%d:sq=%d:lb=%d" % (
      self.gaussian_fit.n_parameters(),
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
      bound_flags = self.gaussian_fit.bound_flags(False, True)
      nbd.set_selected(bound_flags, 1)
    self.minimizer = lbfgsb.minimizer(
      n=self.n,
      m=self.lbfgsb_m,
      l=l,
      u=u,
      nbd=nbd,
      iprint=self.iprint)
    self.x = flex.double(self.n, 0)
    self.f = 0
    self.g = flex.double(self.n, 0)
    while True:
      if (self.minimizer.process(self.x, self.f, self.g,
                                 self.use_fortran_library)):
        self.compute_fg()
      elif (self.minimizer.is_terminated()):
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
    for apply_lower_bounds_on_b in [False, True]:
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

class immoptibox_mixin(minimize_mixin):

  def f(self, x):
    self.x = x
    self.apply_shifts()
    self.x = None
    return self.gaussian_fit_shifted.differences()

  def jacobian(self, x):
    self.x = x
    self.apply_shifts()
    self.x = None
    return self.gaussian_fit_shifted.least_squares_jacobian_abc()

  def gradients(self, x, f_x=None):
    if (f_x is None): f_x = self.f(x=x)
    return self.jacobian(x=x).matrix_transpose().matrix_multiply(f_x)

  def hessian(self, x):
    self.x = x
    self.apply_shifts()
    self.x = None
    return self.gaussian_fit_shifted.least_squares_hessian_abc_as_packed_u() \
      .matrix_packed_u_as_symmetric()

class minimize_levenberg_marquardt(immoptibox_mixin):

  def __init__(self, gaussian_fit, k_max=500, hard_b_min=-1):
    adopt_init_args(self, locals())
    minimize_multi_histogram.setdefault(str(self), 0)
    self.shift_sqrt_b = 0
    self.minimizer = immoptibox_ports.levenberg_marquardt(
      function=self,
      x0=flex.double(self.gaussian_fit.n_parameters(), 0),
      tau=1.e-8,
      eps_1=1.e-16,
      eps_2=1.e-16,
      k_max=k_max)
    self.final_target_value = self.minimizer.f_x_star.norm()**2
    self.finalize()

  def __str__(self):
    return "levenberg_marquardt:np=%d:k_max=%d" % (
      self.gaussian_fit.n_parameters(),
      self.k_max)

class minimize_damped_newton(immoptibox_mixin):

  def __init__(self, gaussian_fit, k_max=500, hard_b_min=-1):
    adopt_init_args(self, locals())
    minimize_multi_histogram.setdefault(str(self), 0)
    self.shift_sqrt_b = 0
    self.minimizer = immoptibox_ports.damped_newton(
      function=self,
      x0=flex.double(self.gaussian_fit.n_parameters(), 0),
      tau=1.e-8,
      eps_1=1.e-16,
      eps_2=1.e-16,
      k_max=k_max)
    self.final_target_value = self.minimizer.f_x_star.norm()**2
    self.finalize()

  def __str__(self):
    return "damped_newton:np=%d:k_max=%d" % (
      self.gaussian_fit.n_parameters(),
      self.k_max)

def show_minimize_multi_histogram(f=None, reset=True):
  global minimize_multi_histogram
  minimizer_types = minimize_multi_histogram.keys()
  counts = flex.double(minimize_multi_histogram.values())
  perm = flex.sort_permutation(data=counts, reverse=True)
  minimizer_types = flex.select(minimizer_types, perm)
  counts = counts.select(perm)
  n_total = flex.sum(counts)
  for m,c in zip(minimizer_types, counts):
    print >> f, "%-39s  %5.3f %6d" % (m, c/max(1,n_total), c)
  print >> f
  if (reset):
    minimize_multi_histogram = {"None": 0}

def minimize_multi(start_fit,
                   target_powers,
                   minimize_using_sigmas,
                   shift_sqrt_b_mod_n,
                   b_min,
                   n_repeats_minimization):
  best_min_list = []
  for immoptibox_minimizer in [minimize_levenberg_marquardt,
                               minimize_damped_newton]:
    try: minimized = immoptibox_minimizer(gaussian_fit=start_fit)
    except LargeNegativeB: pass
    else: best_min_list.append(minimized)
  for current_shift_sqrt_b_mod_n in shift_sqrt_b_mod_n:
    for minimize_multi_type in [minimize_multi_lbfgs, minimize_multi_lbfgsb]:
      try:
        minimized = minimize_multi_type(
          start_fit=start_fit,
          target_powers=target_powers,
          minimize_using_sigmas=minimize_using_sigmas,
          shift_sqrt_b_mod_n=current_shift_sqrt_b_mod_n,
          b_min=b_min,
          n_repeats_minimization=n_repeats_minimization)
      except RuntimeError, e:
        if (str(e).find("SCITBX_ASSERT(b >= 0)") < 0):
          raise
      else: best_min_list.append(minimized)
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
  while True:
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
    b.set_selected(sel, flex.double(sel.count(True), b_min))
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
  if (a[-1] != 0 and x_sq != 0):
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
  for i,target_value in enumerate(null_fit.table_y()):
    if (i_x_begin is None and target_value < y0 * factor_y_x_begin):
      i_x_begin = i
    if (i_x_end is None and target_value < y0 * factor_y_x_end):
      i_x_end = i+1
      break
  if (i_x_end is None):
    i_x_end = null_fit.table_y().size()
  if (i_x_begin is None):
    i_x_begin = min(existing_gaussian.n_parameters(), i_x_end-1)
  assert i_x_end > 0
  assert i_x_begin < i_x_end
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
          if (good_min is not None and good_min.is_better_than(best_min)):
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
                          params,
                          print_to=None):
  if (label is not None and print_to is None):
    print_to = sys.stdout
  good_min = None
  if (print_to is not None):
    print >> print_to, "label:", label
    print_to.flush()
  n_golay_codes_total = 2**12
  n_golay_codes_processed = 0
  for golay_code in golay_24_12_generator():
    n_golay_codes_processed += 1
    start_fit = make_golay_based_start_gaussian(
      null_fit=null_fit,
      code=golay_code)
    best_min = minimize_multi(
      start_fit=start_fit,
      target_powers=params.target_powers,
      minimize_using_sigmas=params.minimize_using_sigmas,
      shift_sqrt_b_mod_n=params.shift_sqrt_b_mod_n,
      b_min=params.b_min,
      n_repeats_minimization=params.n_repeats_minimization)
    if (best_min is not None):
      if (good_min is None or good_min.max_error > best_min.max_error):
        good_min = best_min
        good_min.n_golay_codes_processed = n_golay_codes_processed
        fit_more = scitbx.math.gaussian.fit(
          null_fit_more.table_x(),
          null_fit_more.table_y(),
          null_fit_more.table_sigmas(),
          good_min.final_gaussian_fit)
        if (print_to is not None):
          print >> print_to, label, "max_error fitted=%.4f" % (
            good_min.max_error),
          if (null_fit_more.table_x().size() > null_fit.table_x().size()):
            print >> print_to, "more=%.4f" % (
              flex.max(fit_more.significant_relative_errors())),
          print >> print_to, "Golay code #%d (%.2f%% of all)" % (
            good_min.n_golay_codes_processed,
            100.*good_min.n_golay_codes_processed/n_golay_codes_total)
          fit_more.show(f=print_to)
          fit_more.show_table(f=print_to)
          print >> print_to
        print_to.flush()
        if (good_min.max_error <= params.negligible_max_error):
          break
  if (print_to is not None):
    print >> print_to, "Total number of Golay codes processed:", \
      n_golay_codes_processed
    print >> print_to
    if (good_min is None):
      print >> print_to, "Final: %s: No successful minimization." % label
    else:
      print >> print_to, "Final:", label, "max_error fitted=%.4f" %(
        good_min.max_error),
      if (null_fit_more.table_x().size() > null_fit.table_x().size()):
        print >> print_to, "more=%.4f" % (
          flex.max(fit_more.significant_relative_errors())),
      print >> print_to, "Golay code #%d (%.2f%% of all)" % (
        good_min.n_golay_codes_processed,
        100.*good_min.n_golay_codes_processed/n_golay_codes_total)
    good_min.show_minimization_parameters(f=print_to)
    print >> print_to
    show_minimize_multi_histogram(f=print_to)
  return good_min

def decremental_fit(existing_gaussian, params):
  n_terms = existing_gaussian.n_terms() - 1
  good_min = None
  last_a = list(existing_gaussian.array_of_a())
  last_b = list(existing_gaussian.array_of_b())
  for i_del in xrange(existing_gaussian.n_terms()):
    a_del = last_a[i_del]
    sel_a = last_a[:i_del] + last_a[i_del+1:]
    sel_b = last_b[:i_del] + last_b[i_del+1:]
    for i_add in xrange(n_terms):
      a = sel_a[:]
      a[i_add] += a_del
      start_fit = scitbx.math.gaussian.fit(
        existing_gaussian.table_x(),
        existing_gaussian.table_y(),
        existing_gaussian.table_sigmas(),
        scitbx.math.gaussian.sum(a, sel_b))
      while 1:
        best_min = scitbx.math.gaussian_fit.minimize_multi(
          start_fit=start_fit,
          target_powers=params.target_powers,
          minimize_using_sigmas=params.minimize_using_sigmas,
          shift_sqrt_b_mod_n=params.shift_sqrt_b_mod_n,
          b_min=params.b_min,
          n_repeats_minimization=params.n_repeats_minimization)
        if (best_min is None):
          break
        if (best_min.max_error <= params.max_max_error):
          if (best_min.is_better_than(good_min)):
            good_min = best_min
          break
        start_fit = scitbx.math.gaussian.fit(
          start_fit.table_x()[:-1],
          start_fit.table_y()[:-1],
          start_fit.table_sigmas()[:-1],
          best_min.final_gaussian_fit)
  return good_min
