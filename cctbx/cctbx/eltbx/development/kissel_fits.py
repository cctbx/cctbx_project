from cctbx.eltbx.development import kissel_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from scitbx.python_utils.math_utils import ifloor, iround
from scitbx.python_utils.misc import user_plus_sys_time
from scitbx.python_utils import easy_pickle
from libtbx.optparse_wrapper import OptionParser
import sys, os

def itvc_sampling_selection(kissel_tab):
  xi = cctbx.eltbx.gaussian_fit.international_tables_stols
  xk = kissel_tab.x
  selection = flex.bool(xk.size(), 00000)
  i_kissel = 0
  for i_itvc in xrange(xi.size()):
    while (xk[i_kissel] < xi[i_itvc]):
      i_kissel += 1
    if (xk[i_kissel] == xi[i_itvc]):
      selection[i_kissel] = 0001
    elif (i_kissel > 0 and xk[i_kissel-1] < xi[i_itvc] < xk[i_kissel]):
      if (xi[i_itvc] - xk[i_kissel-1] < xk[i_kissel] - xi[i_itvc]):
        selection[i_kissel-1] = 0001
      else:
        selection[i_kissel] = 0001
  return selection

def run(args, cutoff, six_term=00000, params=None,
        plots_dir="kissel_fits_plots", verbose=0):
  timer = user_plus_sys_time()
  if (params is None):
    params = cctbx.eltbx.gaussian_fit.fit_parameters()
  chunk_n = 1
  chunk_i = 0
  if (len(args) > 0 and len(args[0].split(",")) == 2):
    chunk_n, chunk_i = [int(i) for i in args[0].split(",")]
    args = args[1:]
  if (not os.path.isdir(plots_dir)):
    print "No plots because target directory does not exist (mkdir %s)." % \
      plots_dir
    plots_dir = None
  if (chunk_n > 1):
    assert plots_dir is not None
  results = {}
  results["fit_parameters"] = params
  i_chunk = 0
  for file_name in args:
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    tab = kissel_io.read_table(file_name)
    more_selection = itvc_sampling_selection(tab)
    fit_selection = more_selection & (tab.x <= cutoff + 1.e-6)
    null_fit = scitbx.math.gaussian.fit(
      tab.x.select(fit_selection),
      tab.y.select(fit_selection),
      tab.sigmas.select(fit_selection),
      xray_scattering.gaussian(0, 00000))
    null_fit_more = scitbx.math.gaussian.fit(
      tab.x.select(more_selection),
      tab.y.select(more_selection),
      tab.sigmas.select(more_selection),
      xray_scattering.gaussian(0, 00000))
    if (not six_term):
      results[tab.element] = cctbx.eltbx.gaussian_fit.incremental_fits(
        label=tab.element,
        null_fit=null_fit,
        params=params,
        plots_dir=plots_dir,
        verbose=verbose)
    else:
      print "label:", tab.element
      sys.stdout.flush()
      fit = scitbx.math.gaussian_fit.fit_with_golay_starts(
        label=tab.element,
        null_fit=null_fit,
        null_fit_more=null_fit_more,
        n_terms=6,
        target_powers=params.target_powers,
        minimize_using_sigmas=params.minimize_using_sigmas,
        enforce_positive_b_mod_n=params.enforce_positive_b_mod_n,
        b_min=params.b_min,
        n_repeats_minimization=params.n_repeats_minimization)
      g = fit.min.final_gaussian_fit
      results[tab.element] = [xray_scattering.fitted_gaussian(
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
  parser.add_option("-s", "--six_term",
    action="store_true", dest="six_term", default=0,
    help="fit six-term Gaussians using Golay based starts")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  run(
    args=args,
    cutoff=options.cutoff,
    six_term=options.six_term,
    verbose=options.verbose)

if (__name__ == "__main__"):
  main()
