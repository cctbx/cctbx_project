from libtbx import itertbx
from libtbx.utils import Sorry
from scitbx.array_family import flex

from iotbx.option_parser import option_parser
from iotbx import reflection_file_reader
from iotbx.shelx import crystal_symmetry_from_ins

from cctbx import sgtbx
from cctbx import crystal
from cctbx import maptbx
from cctbx import translation_search
from cctbx import xray
from cctbx import miller

from smtbx.ab_initio import charge_flipping

def run(args, command_name='smtbx.charge_flipping'):
  command_line = (option_parser(
    usage="%s [options] reflection_file" % command_name,
    description="Examples:\n"
                "  %s data.mtz\n"
                "  %s hklf4=data.hkl (implicit --symmetry=data.ins)\n"
                "  %s --symmetry=file.ins data.hkl=hklf4\n"
                "  %s --unit_cell='1 2 3 90. 105.1 90.' "
                "--space_group=P21/n hklf3=data.hkl\n"
                % ((command_name,)*4)
    )
    .enable_symmetry_comprehensive()
    .option(None, "--delta",
      action="store",
      type="float",
      default=None,
      help="Parameter delta",
      metavar="FLOAT")
    .option(None, "--max_iterations",
      action="store",
      type="int",
      default=100,
      help="Number of iterations to perform",
      metavar="INT")
  ).process(args)
  if not command_line.args:
    command_line.parser.show_help()
    return 1
  for file_name in command_line.args:
    reflection_file = reflection_file_reader.any_reflection_file(file_name)
    miller_arrays = None
    if reflection_file.file_type() is not None:
      try:
        miller_arrays = reflection_file.as_miller_arrays(
          crystal_symmetry=command_line.symmetry,
          force_symmetry=True,
          merge_equivalents=False)
      except Sorry, KeyboardInterrupt: raise
      except: pass
    if miller_arrays is None:
      print >> sys.stderr, "Error: unknown file format:", file_name
      print >> sys.stderr
      return 1
    for i,mi in enumerate(miller_arrays):
      if mi.is_xray_intensity_array():
        print >> sys.stderr, "Warning: converting intensities to amplitudes"
        mi = mi.f_sq_as_f()
      elif not mi.is_xray_amplitude_array():
        continue
      mi = mi.merge_equivalents().array()
      solve(mi, command_line, i)

def solve(f_obs, command_line, i=0):
  from libtbx import easy_pickle
  co = command_line.options

  # Loop the charge flipping algorithm
  flipped = charge_flipping.weak_reflection_improved_iterator(f_obs=f_obs,
                                                              delta=co.delta)
  for i,state in enumerate(itertbx.islice(flipped,co.max_iterations)):
    pass
  if i == co.max_iterations - 1:
    print >> sys.stderr, "Maximum number of iterations reached: terminating"
  print "Converged! R_1 = %.4f" % state.r1_factor()

  cleaned = charge_flipping.low_density_elimination_iterator(
    f_obs,
    rho_map=flipped.rho_map)
  for i,state in enumerate(itertbx.islice(cleaned,20)):
    pass
  print "Cleaning: R_1 = %.4f" % state.r1_factor()

  easy_pickle.dump(os.path.expanduser("~/Desktop/charge_flipping.pickle"),
                 (f_obs.crystal_symmetry().unit_cell(),
                  state.g.fft_map().real_map()))

  # Translate the electron density so that its space-group is that of f_obs
  cc_strongest_peak = state.correlation_map_peak_cluster_analysis().next()
  tau = cc_strongest_peak.site
  state.apply_shift(tau)
  symmetric_sf = state.g.customized_copy(crystal_symmetry=f_obs) \
                        .merge_equivalents().array()
  symmetric_map = symmetric_sf.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry)

  easy_pickle.dump(os.path.expanduser("~/Desktop/charge_flipping.pickle"),
                 (f_obs.crystal_symmetry().unit_cell(),
                  symmetric_map.real_map()))

  peaks = symmetric_sf.fft_map(
    symmetry_flags=sgtbx.search_symmetry_flags(
      use_space_group_symmetry=True)
    ).peak_search(
      parameters=maptbx.peak_search_parameters(
        interpolate=True,
        min_distance_sym_equiv=1.5,
        max_clusters=15),
      verify_symmetry=False)

  xs = xray.structure(crystal_symmetry=f_obs.crystal_symmetry())
  for i,s in enumerate(peaks.sites()):
    xs.add_scatterer(xray.scatterer(site=s, label='Q%i' % i,
                                    scattering_type="C"))
  open(os.path.expanduser('~/Desktop/charge_flipping.pdb'),
       'w').write(xs.as_pdb_file())

if __name__ == '__main__':
  import sys
  import os.path
  sys.exit(run(args=sys.argv[1:]))
