from cctbx.eltbx.development import itvc_section61_io
from cctbx.eltbx.development import rez_rez_grant
from cctbx.eltbx.development.create_n_gaussian_raw_cpp import identifier
from cctbx.eltbx import xray_scattering
from scitbx.array_family import flex
import scitbx.math.gaussian_fit
import cctbx.eltbx.gaussian_fit
from libtbx.option_parser import OptionParser
from libtbx.utils import user_plus_sys_time
from libtbx import easy_pickle
import sys, os

def run(file_name, args, cutoff, params,
        zig_zag=False, six_term=False, full_fits=None,
        plots_dir="itvc_fits_plots", verbose=0):
  tab = itvc_section61_io.read_table6111(file_name)
  chunk_n = 1
  chunk_i = 0
  if (len(args) > 0 and len(args[0].split(",")) == 2):
    chunk_n, chunk_i = [int(i) for i in args[0].split(",")]
    args = args[1:]
  if (not six_term and not zig_zag):
    if (not os.path.isdir(plots_dir)):
      print "No plots because target directory does not exist (mkdir %s)." % \
        plots_dir
      plots_dir = None
    if (chunk_n > 1):
      assert plots_dir is not None
  stols_more = cctbx.eltbx.gaussian_fit.international_tables_stols
  sel = stols_more <= cutoff + 1.e-6
  stols = stols_more.select(sel)
  i_chunk = 0
  for element in tab.elements + ["O2-", "SDS"]:
    if (len(args) > 0 and element not in args): continue
    flag = i_chunk % chunk_n == chunk_i
    i_chunk += 1
    if (not flag):
      continue
    results = {}
    results["fit_parameters"] = params
    if (element == "SDS"):
      wrk_lbl = element
      from cctbx.eltbx.xray_scattering.hydrogen_plots import \
        itc_tab_6112_padded
      sds_stols, sds_data = [flex.double(vals)
        for vals in zip(*itc_tab_6112_padded)]
      sds_sigmas = flex.double(sds_data.size(), 0.00005)
      assert sorted(sds_stols) == list(sds_stols)
      sel = sds_stols <= cutoff + 1.e-6
      null_fit = scitbx.math.gaussian.fit(
        sds_stols.select(sel),
        sds_data.select(sel),
        sds_sigmas.select(sel),
        xray_scattering.gaussian(0, False))
      null_fit_more = scitbx.math.gaussian.fit(
        sds_stols,
        sds_data,
        sds_sigmas,
        xray_scattering.gaussian(0, False))
    else:
      wrk_lbl = xray_scattering.wk1995(element, True)
      if (element != "O2-"):
        entry = tab.entries[element]
        null_fit = scitbx.math.gaussian.fit(
          stols,
          entry.table_y[:stols.size()],
          entry.table_sigmas[:stols.size()],
          xray_scattering.gaussian(0, False))
        null_fit_more = scitbx.math.gaussian.fit(
          stols_more,
          entry.table_y[:stols_more.size()],
          entry.table_sigmas[:stols_more.size()],
          xray_scattering.gaussian(0, False))
      else:
        rrg_stols_more = rez_rez_grant.table_2_stol
        sel = rrg_stols_more <= cutoff + 1.e-6
        rrg_stols = rrg_stols_more.select(sel)
        null_fit = scitbx.math.gaussian.fit(
          rrg_stols,
          rez_rez_grant.table_2_o2minus[:rrg_stols.size()],
          rez_rez_grant.table_2_sigmas[:rrg_stols.size()],
          xray_scattering.gaussian(0, False))
        null_fit_more = scitbx.math.gaussian.fit(
          rrg_stols_more,
          rez_rez_grant.table_2_o2minus[:rrg_stols_more.size()],
          rez_rez_grant.table_2_sigmas[:rrg_stols_more.size()],
          xray_scattering.gaussian(0, False))
    if (zig_zag):
      results[wrk_lbl] = cctbx.eltbx.gaussian_fit.zig_zag_fits(
        label=wrk_lbl,
        null_fit=null_fit,
        null_fit_more=null_fit_more,
        params=params)
    elif (full_fits is not None):
      assert len(full_fits.all[wrk_lbl]) == 1
      results[wrk_lbl] = cctbx.eltbx.gaussian_fit.decremental_fits(
        label=wrk_lbl,
        null_fit=null_fit,
        full_fit=full_fits.all[wrk_lbl][0],
        params=params,
        plots_dir=plots_dir,
        verbose=verbose)
    elif (not six_term):
      results[wrk_lbl] = cctbx.eltbx.gaussian_fit.incremental_fits(
        label=wrk_lbl,
        null_fit=null_fit,
        params=params,
        plots_dir=plots_dir,
        verbose=verbose)
    else:
      best_min = scitbx.math.gaussian_fit.fit_with_golay_starts(
        label=wrk_lbl,
        null_fit=null_fit,
        null_fit_more=null_fit_more,
        params=params)
      g = best_min.final_gaussian_fit
      results[wrk_lbl] = [xray_scattering.fitted_gaussian(
        stol=g.table_x()[-1], gaussian_sum=g)]
    sys.stdout.flush()
    pickle_file_name = "%s_fits.pickle" % identifier(wrk_lbl)
    easy_pickle.dump(pickle_file_name, results)

def run_and_time(*args, **kw):
  timer = user_plus_sys_time()
  try:
    run(*args, **kw)
  finally:
    print "CPU time: %.2f seconds" % timer.elapsed()

def main():
  parser = OptionParser(
    usage="usage: python %prog file_name [n_chunks,i_chunk] [scatterer...]")
  parser.add_option("-v", "--verbose",
    action="store_true", default=0,
    help="show comparison table for each element")
  parser.add_option("-c", "--cutoff",
    type="float", default=6, metavar="FLOAT",
    help="maximum sin(theta)/lambda")
  parser.add_option("-q", "--quick",
    action="store_true", default=0,
    help="quick mode for debugging")
  parser.add_option("-n", "--max_n_terms",
    type="int", default=5, metavar="INT",
    help="maximum number of Gaussian terms")
  parser.add_option("-s", "--six_term",
    action="store_true", default=0,
    help="fit six-term Gaussians using Golay based starts")
  parser.add_option("-z", "--zig_zag",
    action="store_true", default=0,
    help="zig-zag fits starting from six-term Gaussians")
  parser.add_option("-r", "--full_fits",
    action="store",
    help="pickled six-term Gaussian fits")
  (options, args) = parser.parse_args()
  if (len(args) < 1):
    parser.print_help()
    return
  assert [options.six_term, options.zig_zag, options.full_fits].count(True) < 2
  if (options.full_fits is not None):
    full_fits = easy_pickle.load(options.full_fits)
  else:
    full_fits = None
  params = cctbx.eltbx.gaussian_fit.fit_parameters(
    max_n_terms=options.max_n_terms)
  if (options.quick):
    params = params.quick()
  run_and_time(
    file_name=args[0],
    args=args[1:],
    cutoff=options.cutoff,
    params=params,
    zig_zag=options.zig_zag,
    six_term=options.six_term,
    full_fits=full_fits,
    verbose=options.verbose)

if (__name__ == "__main__"):
  main()
