from __future__ import absolute_import, division, print_function
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.development.format_gaussian_fits import expected_labels
from cctbx.eltbx.gaussian_fit import international_tables_stols
import scitbx.math.gaussian
from cctbx.array_family import flex
from libtbx.option_parser import OptionParser
from libtbx import easy_pickle
from six.moves import zip

def pick_nicest_fit(fit_0, fit_1):
  if (fit_0.max_error < fit_1.max_error): return fit_0
  if (fit_0.max_error > fit_1.max_error): return fit_1
  a0 = flex.double(fit_0.array_of_a())
  a1 = flex.double(fit_1.array_of_a())
  l0 = (flex.min(a0) > 0)
  l1 = (flex.min(a0) > 0)
  if (l0 and not l1): return fit_0
  if (l1 and not l0): return fit_1
  b0 = flex.double(fit_0.array_of_b())
  b1 = flex.double(fit_1.array_of_b())
  if (flex.min(b0) > flex.min(b1)): return fit_0
  if (flex.min(b0) < flex.min(b1)): return fit_1
  return fit_0

def reset_max_error(itvc_entry, fit):
  sel = international_tables_stols <= fit.stol + 1.e-6
  gaussian_fit = scitbx.math.gaussian.fit(
    international_tables_stols.select(sel),
    itvc_entry.table_y.select(sel),
    itvc_entry.table_sigmas.select(sel),
    fit)
  fit.max_error = flex.max(gaussian_fit.significant_relative_errors())

def main():
  parser = OptionParser(
    usage="usage: python %prog [options]"
         +" itvc_table all_fits_six_terms.pickle"
         +" all_fits_decremental.pickle all_fits_incremental.pickle")
  (options, args) = parser.parse_args()
  if (len(args) != 4):
    parser.print_help()
    return
  itvc_tab = itvc_section61_io.read_table6111(args[0])
  six_term_fits = easy_pickle.load(args[1])
  fits = []
  for file_name in args[2:]:
    fits.append(easy_pickle.load(file_name))
    for label,fit_group in fits[-1].all.items():
      for fit in fit_group:
        reset_max_error(itvc_tab.entries[label], fit)
  best_fits = {}
  n_less = 0
  n_greater = 0
  n_equal = 0
  n_less_list = [0] * 10
  n_greater_list = [0] * 10
  n_equal_list = [0] * 10
  for label in expected_labels(kissel_dir=None):
    fit_group_0 = fits[0].all.get(label, None)
    fit_group_1 = fits[1].all.get(label, None)
    if (fit_group_0 is None and fit_group_1 is None):
      best_fits[label] = None
      continue
    if (fit_group_0 is None):
      best_fits[label] = fit_group_0
      continue
    if (fit_group_1 is None):
      best_fits[label] = fit_group_1
      continue
    best_group = []
    all_n_terms = {}
    n_terms_dicts = []
    for fit_group in [fit_group_0, fit_group_1]:
      n_terms_dict = {}
      for fit in fit_group:
        n_terms_dict[fit.n_terms()] = fit
      n_terms_dicts.append(n_terms_dict)
      all_n_terms.update(n_terms_dicts[-1])
    all_n_terms = list(all_n_terms.keys())
    all_n_terms.sort()
    for n_terms in all_n_terms:
      fit_0 = n_terms_dicts[0].get(n_terms, None)
      fit_1 = n_terms_dicts[1].get(n_terms, None)
      if (fit_0 is None):
        best_group.append(fit_1)
        continue
      if (fit_1 is None):
        best_group.append(fit_0)
        continue
      if   (fit_0.stol < fit_1.stol):
        best_group.append(fit_1)
        status = "less"
        n_less += 1
        n_less_list[n_terms] += 1
      elif (fit_0.stol > fit_1.stol):
        best_group.append(fit_0)
        status = "greater"
        n_greater += 1
        n_greater_list[n_terms] += 1
      else:
        best_group.append(pick_nicest_fit(fit_0, fit_1))
        status = "equal"
        n_equal += 1
        n_equal_list[n_terms] += 1
      print("%-4s n_terms=%d %4.2f %4.2f %s" % (
        label, n_terms, fit_0.stol, fit_1.stol, status))
    best_fits[label] = best_group
  print("n_less:", n_less)
  print("n_greater:", n_greater)
  print("n_equal:", n_equal)
  print("total:", n_less + n_greater + n_equal)
  n_terms = -1
  for n_less,n_greater,n_equal in zip(n_less_list,n_greater_list,n_equal_list):
    n_terms += 1
    if (n_less == 0 and n_greater == 0 and n_equal == 0): continue
    print("n_terms:", n_terms)
    print("  n_less:", n_less)
    print("  n_greater:", n_greater)
    print("  n_equal:", n_equal)
    print("  total:", n_less + n_greater + n_equal)
  print()
  for label in expected_labels(kissel_dir=None):
    if (best_fits[label] is None):
      print("# Warning: Missing scattering_type:", label)
  print()
  print("Best fits:")
  print()
  for label in expected_labels(kissel_dir=None):
    fit_group = best_fits[label]
    if (fit_group is None): continue
    print("scattering_type:", label)
    assert len(six_term_fits.all[label]) == 1
    assert six_term_fits.all[label][0].n_terms() == 6
    fit_group.append(six_term_fits.all[label][0])
    reset_max_error(itvc_tab.entries[label], fit_group[-1])
    trimmed_fit_group =[]
    prev_fit = None
    for fit in fit_group:
      if (prev_fit is None
          or fit.stol > prev_fit.stol
          or (fit.stol == prev_fit.stol
              and fit.max_error < prev_fit.max_error)):
        trimmed_fit_group.append(fit)
        fit.show()
        prev_fit = fit
      else:
        print("# skipped: %s, n_terms: %d, stol: %.2f, max_error: %.4f" % (
          label, fit.n_terms(), fit.stol, fit.max_error))
    best_fits[label] = trimmed_fit_group
  print()
  easy_pickle.dump("best_fits.pickle", best_fits)

if (__name__ == "__main__"):
  main()
