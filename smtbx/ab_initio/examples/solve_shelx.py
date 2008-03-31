from smtbx.ab_initio import charge_flipping
from iotbx import reflection_file_reader
from cctbx import maptbx
from  libtbx.itertbx import izip
try:
  from crys3d import wx_map_viewer
except ImportError:
  wx_map_viewer = None


def run(hkl_path, verbose):
  # Get the reflections from the specified path
  reflections = reflection_file_reader.any_reflection_file(
    'hklf+ins/res=' + hkl_path)
  data = reflections.as_miller_arrays()[0]
  if data.is_xray_intensity_array():
    f_obs = data.f_sq_as_f()

  # merge them (essential!!)
  merging = f_obs.merge_equivalents()
  f_obs = merging.array()
  f_obs.show_summary()

  # charge flipping iterations
  flipping = charge_flipping.basic_iterator(f_obs, delta=None)
  solving = charge_flipping.solving_iterator(
    flipping,
    yield_during_delta_guessing=True,
    yield_solving_interval=1)
  charge_flipping.loop(solving, verbose=verbose)

  # play with the solutions
  if solving.f_calc_solutions:
    # actually only the supposedly best one
    f_calc, shift, cc_peak_height = solving.f_calc_solutions[0]
    fft_map = f_calc.fft_map(
        symmetry_flags=maptbx.use_space_group_symmetry)
    # 3D display of the electron density iso-contours
    if wx_map_viewer is not None:
      wx_map_viewer.display(fft_map)
    # search and print Fourier peaks
    peaks = fft_map.peak_search(
      parameters=maptbx.peak_search_parameters(
        min_distance_sym_equiv=1.0,
        max_clusters=30,),
      verify_symmetry=False
      ).all()
    for q,h in izip(peaks.sites(), peaks.heights()):
      print "(%.3f, %.3f, %.3f) -> %.3f" % (q+(h,))

if __name__ == '__main__':
  from libtbx.option_parser import option_parser
  import sys
  command_line = (option_parser()
                  .option(None, '--verbose',
                          action='store_true')
                  .option(None, '--verbosity',
                          action='store')
                  ).process(args=sys.argv[1:])
  verbose = command_line.options.verbosity or command_line.options.verbose
  run(command_line.args[0], verbose)
