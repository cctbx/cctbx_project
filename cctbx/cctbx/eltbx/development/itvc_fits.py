from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx import xray_scattering
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time
from scitbx.python_utils import easy_pickle
from libtbx.optparse_wrapper import OptionParser
import sys, os

def run(file_name, args, params=None,
        plots_dir="itvc_fits_plots", verbose=0):
  timer = user_plus_sys_time()
  tab = itvc_section61_io.read_table6111(file_name)
  if (params is None): params = cctbx.eltbx.gaussian_fit.fit_parameters()
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
  for element in tab.elements:
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    if (len(args) > 0 and element not in args): continue
    entry = tab.entries[element]
    wk = xray_scattering.wk1995(element, 1)
    null_fit = scitbx.math.gaussian.fit(
      cctbx.eltbx.gaussian_fit.international_tables_stols,
      entry.table_y,
      cctbx.eltbx.gaussian_fit.international_tables_sigmas,
      xray_scattering.gaussian(0, 00000))
    results[wk.label()] = cctbx.eltbx.gaussian_fit.incremental_fits(
      label=element,
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
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  run(file_name=args[0], args=args[1:], verbose=options.verbose)

if (__name__ == "__main__"):
  main()
