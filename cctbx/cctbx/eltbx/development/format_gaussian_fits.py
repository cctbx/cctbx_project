from cctbx.eltbx import xray_scattering
from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.gaussian_fit import international_tables_stols
from cctbx.array_family import flex
import scitbx.math.gaussian
from scitbx.python_utils import easy_pickle
from scitbx.python_utils.misc import adopt_init_args
from libtbx.optparse_wrapper import OptionParser

class labeled_fit:

  def __init__(self, label, gaussian_fit):
    adopt_init_args(self, locals())

class read_pickled_fits:

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

def run(gaussian_fit_pickle_file_names, itvc_file_name):
  itvc_tab = None
  if (itvc_file_name is not None):
    itvc_tab = itvc_section61_io.read_table6111(itvc_file_name)
  fits = read_pickled_fits(gaussian_fit_pickle_file_names)
  for k,v in fits.parameters.items():
    print "# %s:" % k, v
  print
  max_errors = flex.double()
  labeled_fits = []
  n_processed = 0
  for wk in xray_scattering.wk1995_iterator():
    if (wk.label() in ["H'", "D", "Siv"]): continue
    try:
      fit_group = fits.all[wk.label()]
    except:
      print "# Warning: Missing scattering_type:", wk.label()
    else:
      print "scattering_type:", wk.label()
      for fit in fit_group:
        fit.show()
        if (itvc_tab is not None):
          entry = itvc_tab.entries[wk.label()]
          sel = international_tables_stols <= fit.stol + 1.e-6
          gaussian_fit = scitbx.math.gaussian.fit(
            international_tables_stols.select(sel),
            entry.table_y.select(sel),
            entry.table_sigmas.select(sel),
            fit)
          max_errors.append(
            flex.max(gaussian_fit.significant_relative_errors()))
          labeled_fits.append(labeled_fit(wk.label(), gaussian_fit))
          print wk.label(), "n_terms=%d max_error: %.4f" % (
            fit.n_terms(), max_errors[-1])
      n_processed += 1
  print
  assert n_processed == len(fits.all)
  print "Summary:"
  perm = flex.sort_permutation(max_errors, 0001)
  max_errors = max_errors.select(perm)
  labeled_fits = flex.select(labeled_fits, perm)
  for me,lf in zip(max_errors, labeled_fits):
    print lf.label, "n_terms=%d max_error: %.4f" % (
      lf.gaussian_fit.n_terms(), me)
    if (me > 0.01):
      fit = lf.gaussian_fit
      re = fit.significant_relative_errors()
      for s,y,a,r in zip(fit.table_x(),fit.table_y(),fit.fitted_values(),re):
        comment = ""
        if (r > 0.01): comment = " large error"
        print "%4.2f %7.4f %7.4f %7.4f %7.4f%s" % (s,y,a,a-y,r,comment)
      print
  print

def main():
  parser = OptionParser(
    usage="usage: python %prog [options] file_name ...")
  parser.add_option("-t", "--itvc",
    action="store", dest="itvc", metavar="FILE",
    help="file name for international tables data")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  run(
    gaussian_fit_pickle_file_names=args,
    itvc_file_name=options.itvc)

if (__name__ == "__main__"):
  main()
