from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time
from scitbx.python_utils import easy_pickle
from libtbx.optparse_wrapper import OptionParser
import sys, os

def run(file_name, args, cutoff, max_n_terms,
        six_term=00000, full_fits=None, params=None,
        plots_dir="itvc_fits_plots", verbose=0):
  timer = user_plus_sys_time()
  tab = itvc_section61_io.read_table6111(file_name)
  if (params is None):
    params = cctbx.eltbx.gaussian_fit.fit_parameters(
      max_n_terms=max_n_terms)
  chunk_n = 1
  chunk_i = 0
  if (len(args) > 0 and len(args[0].split(",")) == 2):
    chunk_n, chunk_i = [int(i) for i in args[0].split(",")]
    args = args[1:]
  if (not six_term):
    if (not os.path.isdir(plots_dir)):
      print "No plots because target directory does not exist (mkdir %s)." % \
        plots_dir
      plots_dir = None
    if (chunk_n > 1):
      assert plots_dir is not None
  stols_more = cctbx.eltbx.gaussian_fit.international_tables_stols
  sel = stols_more <= cutoff + 1.e-6
  stols = stols_more.select(sel)
  results = {}
  results["fit_parameters"] = params
  i_chunk = 0
  for element in tab.elements:
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    if (len(args) > 0 and element not in args): continue
    wk = xray_scattering.wk1995(element, 1)
    entry = tab.entries[element]
    null_fit = scitbx.math.gaussian.fit(
      stols,
      entry.table_y[:stols.size()],
      entry.table_sigmas[:stols.size()],
      xray_scattering.gaussian(0, 00000))
    null_fit_more = scitbx.math.gaussian.fit(
      stols_more,
      entry.table_y[:stols_more.size()],
      entry.table_sigmas[:stols_more.size()],
      xray_scattering.gaussian(0, 00000))
    if (full_fits is not None):
      assert len(full_fits.all[wk.label()]) == 1
      results[wk.label()] = cctbx.eltbx.gaussian_fit.decremental_fits(
        label=wk.label(),
        null_fit=null_fit,
        full_fit=full_fits.all[wk.label()][0],
        params=params,
        plots_dir=plots_dir,
        verbose=verbose)
    elif (not six_term):
      results[wk.label()] = cctbx.eltbx.gaussian_fit.incremental_fits(
        label=wk.label(),
        null_fit=null_fit,
        params=params,
        plots_dir=plots_dir,
        verbose=verbose)
    else:
      print "label:", wk.label()
      sys.stdout.flush()
      best_min = scitbx.math.gaussian_fit.fit_with_golay_starts(
        label=wk.label(),
        null_fit=null_fit,
        null_fit_more=null_fit_more,
        n_terms=6,
        target_powers=params.target_powers,
        minimize_using_sigmas=params.minimize_using_sigmas,
        enforce_positive_b_mod_n=params.enforce_positive_b_mod_n,
        b_min=params.b_min,
        n_repeats_minimization=params.n_repeats_minimization)
      g = best_min.final_gaussian_fit
      results[wk.label()] = [xray_scattering.fitted_gaussian(
        stol=g.table_x()[-1], gaussian_sum=g)]
    sys.stdout.flush()
    easy_pickle.dump("fits_%02d.pickle" % chunk_i, results)
  print "CPU time: %.2f seconds" % timer.elapsed()

def main():
  parser = OptionParser(
    usage="usage: python %prog file_name [n_chunks,i_chunk] [scatterer...]")
  parser.add_option("-v", "--verbose",
    action="store_true", dest="verbose", default=0,
    help="show comparison table for each element")
  parser.add_option("-c", "--cutoff",
    type="float", dest="cutoff", default=6, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  parser.add_option("-n", "--max_n_terms",
    type="int", dest="max_n_terms", default=5, metavar="INT",
    help="maximum number of Gaussian terms")
  parser.add_option("-s", "--six_term",
    action="store_true", dest="six_term", default=0,
    help="fit six-term Gaussians using Golay based starts")
  parser.add_option("-r", "--full_fits",
    action="store", dest="full_fits",
    help="pickled six-term Gaussian fits")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  assert not (options.six_term and options.full_fits)
  if (options.full_fits is not None):
    full_fits = easy_pickle.load(options.full_fits)
  else:
    full_fits = None
  run(
    file_name=args[0],
    args=args[1:],
    cutoff=options.cutoff,
    max_n_terms=options.max_n_terms,
    six_term=options.six_term,
    full_fits=full_fits,
    verbose=options.verbose)

if (__name__ == "__main__"):
  main()
