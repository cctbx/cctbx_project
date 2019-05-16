from __future__ import absolute_import, division, print_function
from cctbx.eltbx import xray_scattering
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.development import kissel_io
from cctbx.eltbx.gaussian_fit import international_tables_stols
from cctbx.eltbx import tiny_pse
from cctbx.array_family import flex
import scitbx.math.gaussian
from libtbx.option_parser import OptionParser
from libtbx import adopt_init_args
from libtbx import easy_pickle
import os
from six.moves import range
from six.moves import zip

class labeled_fit(object):

  def __init__(self, label, gaussian_fit):
    adopt_init_args(self, locals())

class read_pickled_fits(object):

  def __init__(self, gaussian_fit_pickle_file_names):
    self.parameters = None
    self.all = {}
    for file_name in gaussian_fit_pickle_file_names:
      fits = easy_pickle.load(file_name)
      fp = fits["fit_parameters"].__dict__
      if (self.parameters is None):
        self.parameters = fp
      else:
        for k,v in fp.items():
          assert str(self.parameters[k]) == str(v)
      del fits["fit_parameters"]
      size_before = len(self.all)
      self.all.update(fits)
      assert len(self.all) == size_before + len(fits)

def expected_labels(kissel_dir):
  result = []
  if (kissel_dir is None):
    for wk in xray_scattering.wk1995_iterator():
      result.append(wk.label())
  else:
    for atomic_number in range(1,100):
      result.append(tiny_pse.table(atomic_number).symbol())
  return result

def run(gaussian_fit_pickle_file_names, itvc_file_name, kissel_dir):
  itvc_tab = None
  if (itvc_file_name is not None):
    itvc_tab = itvc_section61_io.read_table6111(itvc_file_name)
  fits = read_pickled_fits(gaussian_fit_pickle_file_names)
  #easy_pickle.dump("all_fits.pickle", fits)
  for k,v in fits.parameters.items():
    print("# %s:" % k, v)
  print()
  max_errors = flex.double()
  labeled_fits = []
  n_processed = 0
  for label in expected_labels(kissel_dir):
    try:
      fit_group = fits.all[label]
    except Exception:
      print("# Warning: Missing scattering_type:", label)
    else:
      print("scattering_type:", label)
      prev_fit = None
      for fit in fit_group:
        if (prev_fit is not None):
          if (fit.stol > prev_fit.stol):
            print("# Warning: decreasing stol")
          elif (fit.stol == prev_fit.stol):
            if (fit.max_error < prev_fit.max_error):
              print("# Warning: same stol but previous has larger error")
        prev_fit = fit
        fit.sort().show()
        gaussian_fit = None
        if (itvc_tab is not None and label != "O2-"):
          entry = itvc_tab.entries[label]
          sel = international_tables_stols <= fit.stol + 1.e-6
          gaussian_fit = scitbx.math.gaussian.fit(
            international_tables_stols.select(sel),
            entry.table_y.select(sel),
            entry.table_sigmas.select(sel),
            fit)
        elif (kissel_dir is not None):
          file_name = os.path.join(kissel_dir, "%02d_%s_rf" % (
            tiny_pse.table(label).atomic_number(), label))
          tab = kissel_io.read_table(file_name)
          sel = tab.itvc_sampling_selection() & (tab.x <= fit.stol + 1.e-6)
          gaussian_fit = scitbx.math.gaussian.fit(
            tab.x.select(sel),
            tab.y.select(sel),
            tab.sigmas.select(sel),
            fit)
        if (gaussian_fit is not None):
          max_errors.append(
            flex.max(gaussian_fit.significant_relative_errors()))
          labeled_fits.append(labeled_fit(label, gaussian_fit))
      n_processed += 1
  print()
  if (n_processed != len(fits.all)):
    print("# Warning: %d fits were not processed." % (
      len(fits.all) - n_processed))
    print()
  if (max_errors.size() > 0):
    print("Summary:")
    perm = flex.sort_permutation(data=max_errors, reverse=True)
    max_errors = max_errors.select(perm)
    labeled_fits = flex.select(labeled_fits, perm)
    quick_summary = {}
    for me,lf in zip(max_errors, labeled_fits):
      print(lf.label, "n_terms=%d max_error: %.4f" % (
        lf.gaussian_fit.n_terms(), me))
      quick_summary[lf.label + "_" + str(lf.gaussian_fit.n_terms())] = me
      if (me > 0.01):
        fit = lf.gaussian_fit
        re = fit.significant_relative_errors()
        for s,y,a,r in zip(fit.table_x(),fit.table_y(),fit.fitted_values(),re):
          comment = ""
          if (r > 0.01): comment = " large error"
          print("%4.2f %7.4f %7.4f %7.4f %7.4f%s" % (s,y,a,a-y,r,comment))
        print()
    print()
    #easy_pickle.dump("quick_summary.pickle", quick_summary)

def cross_check(args):
  quick_summaries = []
  for file_name in args:
    quick_summaries.append(easy_pickle.load(file_name))
  assert len(quick_summaries) == 2
  lines = []
  max_of_errors = flex.double()
  atomic_numbers = flex.double()
  n_less = 0
  n_greater = 0
  n_equal = 0
  for label_1,error_1 in quick_summaries[0].items():
    error_2 = quick_summaries[1].get(label_1, None)
    if (error_2 is not None):
      line = "%-10s %7.4f %7.4f" % (label_1, error_1, error_2)
      if   (error_1 < error_2):
        line += " less    %7.4f" % (error_2/error_1)
        n_less += 1
      elif (error_1 > error_2):
        line += " greater %7.4f" % (error_1/error_2)
        n_greater += 1
      else:
        line += " equal"
        n_equal += 1
      lines.append(line)
      max_of_errors.append(max(error_1, error_2))
      atomic_numbers.append(
        tiny_pse.table(label_1.split("_")[0]).atomic_number())
  for sort_key,reverse in [(max_of_errors,True), (atomic_numbers,False)]:
    perm = flex.sort_permutation(data=sort_key, reverse=reverse)
    perm_lines = flex.select(lines, perm)
    for line in perm_lines:
      print(line)
    print()
  print("n_less:", n_less)
  print("n_greater:", n_greater)
  print("n_equal:", n_equal)
  print("total:", n_less + n_greater + n_equal)

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name ...")
  parser.add_option("-t", "--itvc",
    action="store", metavar="FILE",
    help="file name for international tables data")
  parser.add_option("-k", "--kissel",
    action="store", metavar="DIRECTORY",
    help="directory name for Kissel data")
  parser.add_option("-c", "--cross_check",
    action="store_true", default=0,
    help="compare two quick_summary.pickle files")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  assert not (options.itvc and options.kissel)
  if (options.cross): assert len(args) == 2
  if (not options.cross):
    run(
      gaussian_fit_pickle_file_names=args,
      itvc_file_name=options.itvc,
      kissel_dir=options.kissel)
  else:
    cross_check(args)

def run():
  from cctbx.eltbx.development import format_gaussian_fits
  format_gaussian_fits.main()

if (__name__ == "__main__"):
  run()
