from __future__ import division

import os.path
import sys
import random
import math

from cctbx import sgtbx
from cctbx import miller
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from cctbx import euclidean_model_matching as emma

from iotbx import mtz

from libtbx.test_utils import approx_equal
from libtbx import itertbx
from libtbx.utils import flat_list
from libtbx import group_args
from libtbx.utils import format_cpu_times

import scitbx.matrix as mat

from smtbx.ab_initio import charge_flipping


def randomly_exercise(flipping_type,
                      space_group_info, elements,
                      anomalous_flag,
                      d_min, grid_resolution_factor=1./2,
                      verbose=False
                      ):
  # Generate a random structure in real space, that we will try to recover
  target_structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=False,
    use_u_aniso=False,
  )
  for s in target_structure.scatterers():
    assert s.flags.use_u_iso() and s.u_iso == 0
    assert not s.flags.use_u_aniso()

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
  f_obs = abs(f_target)

  # Unleash charge flipping on the amplitudes
  result = group_args(had_phase_transition=False,
                      emma_match_in_p1=False,
                      emma_match=False,
                      first_correct_correlation_peak=None,
                      inverted_solution=False,
                      succeeded=True)
  flipping = flipping_type(f_obs, delta=None)
  solving = charge_flipping.solving_iterator(
    flipping,
    yield_during_delta_guessing=True,
    yield_solving_interval=1)
  charge_flipping.loop(solving, verbose=verbose)

  # check whether a phase transition has occured
  result.had_phase_transition = solving.had_phase_transition()
  if not result.had_phase_transition:
    result.succeeded = False
    if verbose:
      print "@@ No phase transition @@"
    return result

  flipping = solving.flipping_iterator
  f_result_in_p1 = solving.flipping_iterator.f_calc

  # Euclidean matching of the peaks from the obtained map
  # against those of the correct structure (in P1)
  target_structure_in_p1 = target_structure.expand_to_p1()
  search_parameters = maptbx.peak_search_parameters(
    interpolate=True,
    min_distance_sym_equiv=1.,
    max_clusters=int(target_structure_in_p1.scatterers().size()*1.05))
  peak_search_outcome = flipping.rho_map.peak_search(search_parameters)
  peak_structure = emma.model(
    target_structure_in_p1.crystal_symmetry().special_position_settings(),
    positions=[ emma.position('Q%i' % i, x)
                for i,x in enumerate(peak_search_outcome.all().sites()) ])
  refined_matches = emma.model_matches(
    target_structure_in_p1.as_emma_model(),
    peak_structure,
    break_if_match_with_no_singles=False
    ).refined_matches
  result.emma_match_in_p1 = (
        refined_matches[0].rms < 0.15 # no farther than that
    and refined_matches[0].rt.r in (mat.identity(3), mat.inversion(3))
  )
  if not result.emma_match_in_p1:
    if 0:
      from crys3d import wx_map_viewer
      wx_map_viewer.display(title="Target",
                            fft_map=flipping.rho_map)
    result.succeeded = False
    if verbose:
      print "** no Euclidean matching in P1 **"
    if verbose == "debug":
      print refined_matches[0].show()
  result.inverted_solution = refined_matches[0].rt.r == mat.inversion(3)
  reference_shift = -refined_matches[0].rt.t

  if verbose and not result.succeeded:
    print "@@ Failure @@"
    return result

  # Find the translation to bring back the structure to the same space-group
  # setting as the starting f_obs from correlation map analysis
  is_allowed= lambda x: f_target.space_group_info().is_allowed_origin_shift(
    x, tolerance=0.1)
  for i, (f_calc, shift, cc_peak_height) in enumerate(
                                                    solving.f_calc_solutions):
    if (   is_allowed(shift - reference_shift)
        or is_allowed(shift + reference_shift)):
      result.first_correct_correlation_peak = i
      break
    else:
      if verbose == "more":
        print "++ Incorrect peak: shift=(%.3f, %.3f, %.3f), height=%.2f"\
              % (tuple(shift)+(cc_peak_height,))
        print "   Reference shift=(%.3f, %.3f, %.3f)" % tuple(reference_shift)
  if verbose:
    if result.first_correct_correlation_peak is None:
      result.succeeded = False
      print "** No correct correlation peak **"
    elif result.first_correct_correlation_peak != 0:
      print "** First correct correlation peak: #%i (%.3f) **"\
            % (result.first_correct_correlation_peak, cc_peak_height)

  if verbose and not result.succeeded:
    print "@@ Failure @@"
    return result

  # check Euclidean matching in the original space-group
  search_parameters = maptbx.peak_search_parameters(
    interpolate=True,
    min_distance_sym_equiv=1.,
    max_clusters=int(target_structure.scatterers().size()*1.05))
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
  if refined_matches:
    result.emma_match = (
          not refined_matches[0].singles1 # all sites match a peak
      and refined_matches[0].rms < 0.1 # no farther than that
      and refined_matches[0].rt.r in (mat.identity(3), mat.inversion(3))
    )
  else:
    result.emma_match = False
  if not result.emma_match:
    result.succeeded = False
    if verbose:
      print "** no Euclidean matching in original spacegroup **"

  # return a bunch of Boolean flags telling where it failed
  # or whether it succeeded
  if verbose:
    if result.succeeded:
      if result.first_correct_correlation_peak > 0:
        print "@@ Success (but not first solution) @@"
      else:
        print "@@ Success @@"
    else:
      print "@@ Failure @@"
  return result

def exercise(flags, space_group_info, flipping_type):
  if not flags.repeats: flags.repeats = 1
  results = []
  n = len(space_group_info.group())
  if n > 48: return False
  n_C = 12//n or 1
  n_O = 6//n
  n_N = 3//n
  #n_C = 3
  #n_O = 0
  #n_N = 0
  if flags.Verbose:
    print "C%i O%i N%i" % (n_C*n, n_O*n, n_N*n)
  for i in xrange(int(flags.repeats)):
    result = randomly_exercise(
      flipping_type=flipping_type,
      space_group_info=space_group_info,
      elements=["C"]*n_C + ["O"]*n_O + ["N"]*n_N,
      anomalous_flag=False,
      d_min=0.8,
      verbose=flags.Verbose
    )
    results.append(result)
  if flags.Verbose: print
  return True, results

def exercise_charge_flipping():
  import sys
  results = []
  results.extend(
    debug_utils.parse_options_loop_space_groups(
      sys.argv[1:],
      exercise,
      keywords=("repeats",),
      flipping_type=charge_flipping.weak_reflection_improved_iterator,
    ))
  results = flat_list([ x for x in results if x is not None ])
  n_tests = len(results)

  n_phase_transitions = [
    r.had_phase_transition for r in results ].count(True)
  n_emma_matches_in_p1 = [
    r.emma_match_in_p1 for r in results ].count(True)
  n_emma_matches = [
    r.emma_match for r in results ].count(True)
  n_found_shift = [
    r.first_correct_correlation_peak is not None for r in results ].count(True)
  n_success = [ r.succeeded for r in results ].count(True)
  print "%i test cases:" % n_tests
  print "\t%i successes" % n_success
  print "\t%i phase transitions" % n_phase_transitions
  print ("\t%i Euclidean matches with correct structure "
         "in P1" % n_emma_matches_in_p1)
  print "\t%i found shifts" % n_found_shift
  print ("\t%i Euclidean matches with correct structure "
         "in original spacegroup" % n_emma_matches)
  assert n_success/n_tests > 3/4

def run():
  exercise_charge_flipping()

if __name__ == '__main__':
  run()
