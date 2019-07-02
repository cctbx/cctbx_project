from __future__ import absolute_import, division, print_function

import sys
import random

from cctbx import miller
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from cctbx import euclidean_model_matching as emma

from libtbx import group_args

import scitbx.matrix as mat
from six.moves import cStringIO as StringIO

from smtbx.ab_initio import charge_flipping
from six.moves import range

def randomly_exercise(flipping_type,
                      space_group_info, elements,
                      anomalous_flag,
                      d_min, grid_resolution_factor=1./2,
                      verbose=False,
                      amplitude_type="F",
                      ):

  # Generate a random structure in real space, that we will try to recover
  target_structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=True,
    random_u_iso=True,
    random_u_iso_scale=0.04,
    use_u_aniso=False,
  )
  exercise_one_structure(target_structure,
                         flipping_type,
                         anomalous_flag,
                         d_min,
                         grid_resolution_factor=grid_resolution_factor,
                         verbose=verbose,
                         amplitude_type=amplitude_type,
                         )

def exercise_one_structure(target_structure,
                           flipping_type,
                           anomalous_flag,
                           d_min, grid_resolution_factor=1./2,
                           verbose=False,
                           amplitude_type="F",
                           ):
  assert amplitude_type in ('F', 'E', 'quasi-E')

  # Generate its structure factors
  f_target = miller.build_set(
    crystal_symmetry=target_structure,
    anomalous_flag=target_structure.scatterers().count_anomalous() != 0,
    d_min=d_min
    ).structure_factors_from_scatterers(
      xray_structure=target_structure,
      algorithm="direct").f_calc()

  f_target_in_p1 = f_target.expand_to_p1()\
                           .as_non_anomalous_array()\
                           .merge_equivalents().array()
  f_obs = f_target.as_amplitude_array()

  # Unleash charge flipping on the amplitudes
  flipping = flipping_type(delta=None)
  extra = group_args()
  if amplitude_type == 'E':
    extra.normalisations_for = lambda f: f.amplitude_normalisations(
      target_structure.unit_cell_content(omit=('H','D')))
  elif amplitude_type == 'quasi-E':
    extra.normalisations_for = charge_flipping.amplitude_quasi_normalisations
  solving = charge_flipping.solving_iterator(
    flipping,
    f_obs,
    yield_during_delta_guessing=True,
    yield_solving_interval=1,
    **extra.__dict__
  )
  s = StringIO()
  charge_flipping.loop(solving, verbose="highly", out=s)
  if verbose:
    print(s.getvalue())

  # check whether a phase transition has occured
  assert solving.had_phase_transition

  flipping = solving.flipping_iterator
  f_result_in_p1 = solving.flipping_iterator.f_calc

  # Euclidean matching of the peaks from the obtained map
  # against those of the correct structure (in P1)
  target_structure = target_structure.select(
    target_structure.scattering_types() == "H", negate=True)
  target_structure_in_p1 = target_structure.expand_to_p1()
  search_parameters = maptbx.peak_search_parameters(
    interpolate=True,
    min_distance_sym_equiv=1.,
    max_clusters=int(target_structure_in_p1.scatterers().size()*1.2))
  peak_search_outcome = flipping.rho_map.peak_search(search_parameters)
  peak_structure = emma.model(
    target_structure_in_p1.crystal_symmetry().special_position_settings(),
    positions=[ emma.position('Q%i' % i, x)
                for i,x in enumerate(peak_search_outcome.all().sites()) ])
  refined_matches = emma.model_matches(
    target_structure_in_p1.as_emma_model(),
    peak_structure,
    tolerance=0.5,
    break_if_match_with_no_singles=False
    ).refined_matches
  m = refined_matches[0]
  assert m.rms < 0.2, m.rms # no farther than that
  assert m.rt.r in (mat.identity(3), mat.inversion(3))

  reference_shift = -refined_matches[0].rt.t

  # Find the translation to bring back the structure to the same space-group
  # setting as the starting f_obs from correlation map analysis
  is_allowed= lambda x: f_target.space_group_info().is_allowed_origin_shift(
    x, tolerance=0.1)
  first_correct_correlation_peak = None
  for i, (f_calc, shift, cc_peak_height) in enumerate(
                                                    solving.f_calc_solutions):
    if (   is_allowed(shift - reference_shift)
        or is_allowed(shift + reference_shift)):
      first_correct_correlation_peak = i
      break
    else:
      if verbose == "more":
        print("++ Incorrect peak: shift=(%.3f, %.3f, %.3f), height=%.2f"\
              % (tuple(shift)+(cc_peak_height,)))
        print("   Reference shift=(%.3f, %.3f, %.3f)" % tuple(reference_shift))
  assert first_correct_correlation_peak is not None
  if verbose and first_correct_correlation_peak != 0:
      print("** First correct correlation peak: #%i (%.3f) **"\
            % (first_correct_correlation_peak, cc_peak_height))

  # check Euclidean matching in the original space-group
  search_parameters = maptbx.peak_search_parameters(
    interpolate=True,
    min_distance_sym_equiv=1.,
    max_clusters=int(1.5*target_structure.scatterers().size()))
  solution_fft_map = f_calc.fft_map(
    symmetry_flags=maptbx.use_space_group_symmetry)
  solution_peaks = solution_fft_map.peak_search(search_parameters,
                                                verify_symmetry=False)
  solution_peak_structure = emma.model(
    target_structure.crystal_symmetry().special_position_settings(),
    positions=[ emma.position('Q%i' % i, x)
                for i,x in enumerate(solution_peaks.all().sites()) ])
  refined_matches = emma.model_matches(
    target_structure.as_emma_model(),
    solution_peak_structure,
    break_if_match_with_no_singles=False
    ).refined_matches
  assert refined_matches
  m = refined_matches[0]
  assert not m.singles1, m.show() # all sites match a peak
  assert m.rms < 0.15, m.rms  # no farther than that
  assert m.rt.r in (mat.identity(3), mat.inversion(3))

  # success!
  if verbose:
    print("@@ Success @@")

def exercise_sucrose(flipping_type,
                     anomalous_flag,
                     d_min,
                     verbose=False,
                     amplitude_type='quasi-E'):
  from smtbx import development
  target_structure = development.sucrose()

  print("Sucrose")
  exercise_one_structure(target_structure,
                         flipping_type,
                         anomalous_flag,
                         d_min,
                         grid_resolution_factor=1/2,
                         verbose=verbose,
                         amplitude_type=amplitude_type,
                         )

def exercise(flags, space_group_info):
  if not flags.repeats: flags.repeats = 1
  if not flags.algo: flags.algo = "weak_reflection_improved"
  if not flags.on: flags.on = "E"
  if flags.fix_seed:
    random.seed(1)
    flex.set_random_seed(1)

  n = len(space_group_info.group())
  print(space_group_info.type().hall_symbol(), end=' ')
  if not flags.high_symmetry and n > 24:
    print('  [ skipped ]')
    if flags.Verbose: print()
    return
  else:
    print()
  n_C = 12//n or 1
  n_O = 6//n
  n_N = 3//n
  if flags.Verbose:
    print("unit cell content: C%i O%i N%i" % (n_C*n, n_O*n, n_N*n))
    print("asu content: C%i O%i N%i" % (n_C, n_O, n_N))
    print("on %s's with %s" % (flags.on, flags.algo))
  flipping_type = eval("charge_flipping.%s_iterator" % flags.algo)
  for i in range(int(flags.repeats)):
    randomly_exercise(
      flipping_type=flipping_type,
      space_group_info=space_group_info,
      elements=["C"]*n_C + ["O"]*n_O + ["N"]*n_N,
      anomalous_flag=False,
      d_min=0.8,
      verbose=flags.Verbose,
      amplitude_type=flags.on
    )
  if flags.Verbose: print()

def exercise_charge_flipping():
  #exercise_sucrose(flipping_type=charge_flipping.weak_reflection_improved_iterator,
                   #anomalous_flag=False,
                   #d_min=0.7)
  import sys
  debug_utils.parse_options_loop_space_groups(
    sys.argv[1:],
    exercise,
    keywords=("repeats", 'on', 'algo', 'fix_seed', 'high_symmetry'),
    symbols_to_stderr=False,
  )

def run():
  exercise_charge_flipping()

if __name__ == '__main__':
  run()
