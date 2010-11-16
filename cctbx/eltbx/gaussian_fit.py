from cctbx.eltbx import xray_scattering
from cctbx.array_family import flex
import scitbx.math.gaussian
from scitbx.math.gaussian_fit import find_max_x_multi
from scitbx.math.gaussian_fit import show_minimize_multi_histogram
from libtbx import adopt_init_args
import sys, os

# d = 1/(2*stol)
# stol = 1/(2*d)
international_tables_stols = flex.double(
  [0.00, 0.01, 0.02, 0.03, 0.04, 0.05,
   0.06, 0.07, 0.08, 0.09, 0.10,
   0.11, 0.12, 0.13, 0.14, 0.15,
   0.16, 0.17, 0.18, 0.19, 0.20,
   0.22, 0.24, 0.25, 0.26, 0.28, 0.30,
   0.32, 0.34, 0.35, 0.36, 0.38, 0.40,
   0.42, 0.44, 0.45, 0.46, 0.48, 0.50,
   0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00,
   1.10, 1.20, 1.30, 1.40, 1.50,
   1.60, 1.70, 1.80, 1.90, 2.00,
   2.50, 3.00, 3.50, 4.00, 5.00, 6.00])

def show_fit_summary(source, label, gaussian_fit, e,
                     e_other=None, n_terms_other=None):
  n_terms = str(gaussian_fit.n_terms())
  if (gaussian_fit.c() != 0): n_terms += "+c"
  n_terms += ","
  stol = gaussian_fit.table_x()[-1]
  d_min = 1/(2*stol)
  print "%24s: %s n_terms=%-4s stol=%.2f, d_min=%.2f, e=%.4f" % (
    source, label, n_terms, stol, d_min, e),
  if (e_other is not None and e_other > e and n_terms_other == n_terms):
    print "Better",
  print

def show_literature_fits(label, n_terms, null_fit, n_points, e_other=None):
  for lib in [xray_scattering.wk1995,
              xray_scattering.it1992,
              xray_scattering.two_gaussian_agarwal_isaacs,
              xray_scattering.two_gaussian_agarwal_1978,
              xray_scattering.one_gaussian_agarwal_1978]:
    if (lib == xray_scattering.wk1995):
      try:
        lib_gaussian = xray_scattering.wk1995(label, True).fetch()
        lib_source = "WK1995"
      except:
        lib_gaussian = None
    elif (lib == xray_scattering.it1992):
      try:
        lib_gaussian = xray_scattering.it1992(label, True).fetch()
        lib_source = "IT1992"
      except:
        lib_gaussian = None
    elif (lib.table.has_key(label)):
      lib_gaussian = lib.table[label]
      lib_source = lib.source_short
    else:
      lib_gaussian = None
    if (lib_gaussian is not None):
      gaussian_fit = scitbx.math.gaussian.fit(
        null_fit.table_x()[:n_points],
        null_fit.table_y()[:n_points],
        null_fit.table_sigmas()[:n_points],
        lib_gaussian)
      e = flex.max(gaussian_fit.significant_relative_errors())
      show_fit_summary(lib_source, label, gaussian_fit, e,
                       e_other, lib_gaussian.n_terms())

def write_plot(f, xs, ys):
  for x,y in zip(xs,ys):
    print >> f, x, y
  print >> f, "&"

def write_plots(plots_dir, label, gaussian_fit):
  label = label.replace("'", "prime")
  file_name = os.path.join(plots_dir, label+".xy")
  f = open(file_name, "w")
  write_plot(f, gaussian_fit.table_x(), gaussian_fit.table_y())
  write_plot(f, gaussian_fit.table_x(), gaussian_fit.fitted_values())
  f.close()
  file_name = os.path.join(plots_dir, label+".cmp")
  f = open(file_name, "w")
  for x,y,a,e in zip(gaussian_fit.table_x(),
                     gaussian_fit.table_y(),
                     gaussian_fit.fitted_values(),
                     gaussian_fit.significant_relative_errors()):
    print >> f, "%4.2f %7.4f %7.4f %8.5f %7.4f" % (x, y, a, a-y, e)
  f.close()
  return file_name

class fit_parameters(object):

  def __init__(self, max_n_terms=5,
                     target_powers=[2,4],
                     minimize_using_sigmas=False,
                     n_repeats_minimization=5,
                     shift_sqrt_b_mod_n=[0,1,2],
                     b_min=1.e-6,
                     max_max_error=0.01,
                     n_start_fractions=5,
                     negligible_max_error=0.001):
    adopt_init_args(self, locals())

  def quick(self):
    return fit_parameters(
      max_n_terms=2,
      target_powers=[2],
      minimize_using_sigmas=self.minimize_using_sigmas,
      n_repeats_minimization=1,
      shift_sqrt_b_mod_n=[1],
      b_min=self.b_min,
      max_max_error=0.02,
      n_start_fractions=2,
      negligible_max_error=0.02)

def incremental_fits(label, null_fit, params=None, plots_dir=None, verbose=0):
  if (params is None): params = fit_parameters()
  f0 = null_fit.table_y()[0]
  results = []
  previous_n_points = 0
  existing_gaussian = xray_scattering.gaussian([],[])
  while (existing_gaussian.n_terms() < params.max_n_terms):
    if (previous_n_points == null_fit.table_x().size()):
      print "%s: Full fit with %d terms. Search stopped." % (
        label, existing_gaussian.n_terms())
      print
      break
    n_terms = existing_gaussian.n_terms() + 1
    best_min = find_max_x_multi(
      null_fit=null_fit,
      existing_gaussian=existing_gaussian,
      target_powers=params.target_powers,
      minimize_using_sigmas=params.minimize_using_sigmas,
      n_repeats_minimization=params.n_repeats_minimization,
      shift_sqrt_b_mod_n=params.shift_sqrt_b_mod_n,
      b_min=params.b_min,
      max_max_error=params.max_max_error,
      n_start_fractions=params.n_start_fractions)
    if (best_min is None):
      print "Warning: No fit: %s n_terms=%d" % (label, n_terms)
      print
      break
    if (previous_n_points > best_min.final_gaussian_fit.table_x().size()):
      print "Warning: previous fit included more sampling points."
    previous_n_points = best_min.final_gaussian_fit.table_x().size()
    show_fit_summary(
      "Best fit", label, best_min.final_gaussian_fit, best_min.max_error)
    show_literature_fits(
      label=label,
      n_terms=n_terms,
      null_fit=null_fit,
      n_points=best_min.final_gaussian_fit.table_x().size(),
      e_other=best_min.max_error)
    best_min.final_gaussian_fit.show()
    best_min.show_minimization_parameters()
    existing_gaussian = best_min.final_gaussian_fit
    print
    show_minimize_multi_histogram()
    sys.stdout.flush()
    if (plots_dir):
      write_plots(
        plots_dir=plots_dir,
        label=label+"_%d"%n_terms,
        gaussian_fit=best_min.final_gaussian_fit)
    g = best_min.final_gaussian_fit
    results.append(xray_scattering.fitted_gaussian(
      stol=g.table_x()[-1], gaussian_sum=g))
  return results

def decremental_fits(label, null_fit, full_fit=None, params=None,
                     plots_dir=None, verbose=0):
  if (params is None): params = fit_parameters()
  results = []
  last_fit = scitbx.math.gaussian.fit(
    null_fit.table_x(),
    null_fit.table_y(),
    null_fit.table_sigmas(),
    full_fit)
  while (last_fit.n_terms() > 1):
    good_min = scitbx.math.gaussian_fit.decremental_fit(
      existing_gaussian=last_fit,
      params=params)
    if (good_min is None):
      print "%s n_terms=%d: No successful minimization. Aborting." % (
        label, last_fit.n_terms()-1)
      break
    show_fit_summary(
      "Best fit", label, good_min.final_gaussian_fit, good_min.max_error)
    show_literature_fits(
      label=label,
      n_terms=good_min.final_gaussian_fit.n_terms(),
      null_fit=null_fit,
      n_points=good_min.final_gaussian_fit.table_x().size(),
      e_other=good_min.max_error)
    good_min.final_gaussian_fit.show()
    good_min.show_minimization_parameters()
    last_fit = good_min.final_gaussian_fit
    print
    show_minimize_multi_histogram()
    sys.stdout.flush()
    if (plots_dir):
      write_plots(
        plots_dir=plots_dir,
        label=label+"_%d"%good_min.final_gaussian_fit.n_terms(),
        gaussian_fit=good_min.final_gaussian_fit)
    g = good_min.final_gaussian_fit
    results.append(xray_scattering.fitted_gaussian(
      stol=g.table_x()[-1], gaussian_sum=g))
  return results

def zig_zag_fits(label, null_fit, null_fit_more, params):
  six_term_best_min = scitbx.math.gaussian_fit.fit_with_golay_starts(
    label=label,
    null_fit=null_fit,
    null_fit_more=null_fit_more,
    params=params)
  results = []
  n_term_best_min = six_term_best_min
  have_all_points_in_previous = True
  while 1:
    while 1:
      if (n_term_best_min.final_gaussian_fit.n_terms() == 1):
        existing_gaussian = null_fit
      else:
        decr_best_min = scitbx.math.gaussian_fit.decremental_fit(
          existing_gaussian=n_term_best_min.final_gaussian_fit,
          params=params)
        assert decr_best_min is not None
        print "Decremental:",
        decr_best_min.show_summary()
        existing_gaussian = decr_best_min.final_gaussian_fit
        if (n_term_best_min.max_error <= params.negligible_max_error
            and have_all_points_in_previous):
          break
      incr_best_min = find_max_x_multi(
        null_fit=null_fit,
        existing_gaussian=existing_gaussian,
        target_powers=params.target_powers,
        minimize_using_sigmas=params.minimize_using_sigmas,
        n_repeats_minimization=params.n_repeats_minimization,
        shift_sqrt_b_mod_n=params.shift_sqrt_b_mod_n,
        b_min=params.b_min,
        max_max_error=params.max_max_error,
        n_start_fractions=params.n_start_fractions)
      assert incr_best_min is not None
      print "Incremental:",
      incr_best_min.show_summary()
      if (existing_gaussian is null_fit):
        break
      if (not incr_best_min.is_better_than(n_term_best_min)):
        break
      n_term_best_min = incr_best_min
    print " Settled on:",
    n_term_best_min.show_summary()
    results.append(xray_scattering.fitted_gaussian(
      stol=n_term_best_min.final_gaussian_fit.table_x()[-1],
      gaussian_sum=n_term_best_min.final_gaussian_fit,
      max_error=n_term_best_min.max_error))
    if (existing_gaussian is null_fit):
      break
    have_all_points_in_previous = (
         n_term_best_min.final_gaussian_fit.table_x().size()
      == decr_best_min.final_gaussian_fit.table_x().size())
    n_term_best_min = decr_best_min
  print
  return results
