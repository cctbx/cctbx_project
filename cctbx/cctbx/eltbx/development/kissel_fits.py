from cctbx.eltbx.development import kissel_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time
from scitbx.python_utils import easy_pickle
from libtbx.optparse_wrapper import OptionParser
import sys, os

def run(args, cutoff, params=None,
        plots_dir="kissel_fits_plots", verbose=0):
  timer = user_plus_sys_time()
  if (params is None):
    params = cctbx.eltbx.gaussian_fit.fit_parameters(
      table_rounding_error=1.e-5)
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
    sel = tab.x <= cutoff
    tab_x = tab.x.select(sel)
    tab_y = tab.y.select(sel)
    null_fit = scitbx.math.gaussian.fit(
      tab_x,
      tab_y,
      flex.double(tab_x.size(), 1.e-5),
      xray_scattering.gaussian(0, 00000))
    results[tab.element] = cctbx.eltbx.gaussian_fit.incremental_fits(
      label=tab.element,
      null_fit=null_fit,
      params=params,
      plots_dir=plots_dir,
      verbose=verbose)
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
    type="float", dest="cutoff", default=6.05, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  run(args=args, cutoff=options.cutoff, verbose=options.verbose)

if (__name__ == "__main__"):
  main()
