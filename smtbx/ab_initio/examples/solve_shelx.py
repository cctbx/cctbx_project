import re, os
from smtbx.ab_initio import charge_flipping
from iotbx import reflection_file_reader
from iotbx.shelx.from_ins import from_ins
from cctbx import maptbx
from cctbx import euclidean_model_matching as emma

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
    if 0 and wx_map_viewer is not None:
      wx_map_viewer.display(fft_map)
    # search and print Fourier peaks
    nb_peaks_to_search = int(f_calc.unit_cell().volume()/18.6
                          / len(f_calc.space_group()))
    peaks = fft_map.peak_search(
      parameters=maptbx.peak_search_parameters(
        min_distance_sym_equiv=1.0,
        max_clusters=nb_peaks_to_search,),
      verify_symmetry=False
      ).all()
    for q,h in izip(peaks.sites(), peaks.heights()):
      print "(%.3f, %.3f, %.3f) -> %.3f" % (q+(h,))
    if verbose:
      shelx = re.sub(r'(\.(hkl|gz))+$', '', hkl_path)
      ins, res = shelx+'.ins', shelx+'.res'
      if os.path.isfile(ins): shelx = ins
      elif os.path.isfile(res): shelx = res
      else: raise RuntimeError('How did we manage to get here?')
      xs = from_ins(shelx)
      solution_peak_structure = emma.model(
        xs.special_position_settings(),
        positions=[ emma.position('Q%i' % i, x)
                    for i,x in enumerate(peaks.sites()) ])
      refined_matches = emma.model_matches(
        xs.as_emma_model(),
        solution_peak_structure,
        break_if_match_with_no_singles=True
        ).refined_matches
      assert refined_matches
      print refined_matches[0].show()
      

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
