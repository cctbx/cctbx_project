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

import scitbx.matrix as mat

from smtbx.ab_initio import charge_flipping


def random_structure_factors(space_group_info, elements,
                             anomalous_flag,
                             d_min, grid_resolution_factor
                             ):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=False,
    use_u_aniso=False,
  )
  for s in structure.scatterers():
    assert s.flags.use_u_iso() and s.u_iso == 0
    assert not s.flags.use_u_aniso()

  return structure


def find_delta(rho_map, tol):
  """ Find delta as hinted on fig. 1 of ref. [1] in module charge_flipping """
  rho = rho_map.real_map_unpadded().as_1d()
  max_rho = flex.max(rho)
  rho /= max_rho
  sorting = flex.sort_permutation(rho)
  sorted_rho = rho.select(sorting)
  n = len(sorted_rho)
  p,q = n//4, 3*n//4
  indexes = flex.double_range(p,q)
  values = sorted_rho[p:q]
  c = flex.linear_correlation(indexes, values)
  assert c.is_well_defined() and c.coefficient() > 0.99
  r = flex.linear_regression(indexes, values)
  a,b = r.y_intercept(), r.slope()
  deviation = flex.abs(a + b*flex.double_range(n) - sorted_rho)
  non_linear_sel = deviation > tol
  low = flex.first_index(non_linear_sel, False)
  high = flex.last_index(non_linear_sel, False)
  assert non_linear_sel[low:high].count(False)/(high-low+1) > 0.99
  assert sorted_rho[low] < 0 and sorted_rho[high] > 0
  return min(sorted_rho[high], -sorted_rho[low]), max_rho

def write_sorted_moduli_as_mathematica_plot(f, filename):
  """ To obtain fig. 1 in ref [2] in module charge_flipping """
  abs_f = flex.abs(f.data())
  sorted = abs_f.select(flex.sort_permutation(abs_f))
  sorted /= flex.max(sorted)
  mf = open(os.path.expanduser(filename), 'w')
  print >> mf, 'fp1 = {'
  for f in sorted:
    print >> mf, "%f, " % f
  print >> mf, "1 };"
  print >> mf, "ListPlot[fp1]"
  mf.close()


def randomly_exercise(flipping_type,
                      space_group_info, elements,
                      anomalous_flag,
                      d_min, grid_resolution_factor=1./2,
                      verbose=False
                      ):
  target_structure = random_structure_factors(
    space_group_info, elements,
    anomalous_flag,
    d_min, grid_resolution_factor)
  return exercise_once(flipping_type=flipping_type,
                       target_structure=target_structure,
                       f_obs=None ,
                       d_min=d_min,
                       verbose=verbose)

def exercise_once(flipping_type,
                  target_structure, f_obs, d_min, verbose=False):
  result = group_args(had_phase_transition=False,
                      emma_match=False,
                      sharp_correlation_map=False,
                      first_correct_correlation_peak=None,
                      inverted_solution=False)

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
  if f_obs is None:
    f_obs = abs(f_target)

  flipping = flipping_type(f_obs, delta=None)
  solving = charge_flipping.solving_iterator(flipping,
                                             yield_during_delta_guessing=True,
                                             yield_solving_interval=1)

  while solving.state is not solving.finished:
    # Guessing a value of delta leading to subsequent good convergence
    if verbose=="highly":
      print "Guessing delta..."
      print "%10s | %10s | %10s | %10s | %10s" % ('delta', 'R', 'c_tot',
                                                  'c_flip', 'c_tot/c_flip')
      print "-"*64
    for flipping in solving:
      assert solving.state in (solving.guessing_delta, solving.solving)
      if verbose=="highly":
        rho = flipping.rho_map
        c_tot = rho.c_tot()
        c_flip = rho.c_flip(flipping.delta)
        # to compare with superflip output
        c_tot *= flipping.fft_scale; c_flip *= flipping.fft_scale
        print "%10.4f | %10.3f | %10.1f | %10.1f | %10.2f"\
              % (flipping.delta, flipping.r1_factor(),
                 c_tot, c_flip, c_tot/c_flip)
        if solving.state is not solving.guessing_delta: break

    # main charge flipping loop to solve the structure
    if verbose=="highly":
      print
      print "Solving..."
      print "with delta=%.4f" % flipping.delta
      print
      print "%5s | %10s | %10s" % ('#', 'R1', 'c_tot/c_flip')
      print '-'*33
    for flipping in solving:
      assert solving.state in (solving.solving,
                               solving.guessing_delta,
                               solving.polishing)
      r1 = flipping.r1_factor()
      r = flipping.c_tot_over_c_flip()
      if verbose=="highly":
        print "%5i | %10.3f | %10.3f" % (solving.iteration_index, r1, r)
      if solving.state is not solving.solving:
        break

    # need to restart
    if solving.state is solving.guessing_delta:
      msg = "%s Restarting %s" % (('****',)*2)
      print
      print "*" * len(msg)
      print msg
      print "*" * len(msg)
      print
      continue

    # polishing
    for flipping in solving:
      assert solving.state in (solving.polishing, solving.finished)
      if solving.state is solving.finished: break

  # check whether a phase transition has occured
  result.had_phase_transition = (solving.state == solving.finished
                                 and not solving.max_attemps_exceeded)
  if not result.had_phase_transition:
    if verbose: print "** no phase transition **"
    return result

  flipping = solving.flipping_iterator
  f_result_in_p1 = solving.flipping_iterator.f_calc

  # Euclidean matching of the peaks from the obtained map
  # against those of the correct structure
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
    break_if_match_with_no_singles=True
    ).refined_matches
  result.emma_match = (
        not refined_matches[0].singles1 # all sites match a peak
    and refined_matches[0].rms < 0.1 # no farther than that
    and refined_matches[0].rt.r in (mat.identity(3), mat.inversion(3))
  )
  if not result.emma_match:
    if verbose:
      print "** no Euclidean matching **"
    return result
  result.inverted_solution = refined_matches[0].rt.r == mat.inversion(3)
  reference_shift = -refined_matches[0].rt.t

  # Find the translation to bring back the structure to the same space-group
  # setting as the starting f_obs from correlation map analysis
  correlation_map_peak_search, correlation_map = flipping.search_origin(
    return_correlation_map_too=True)
  is_allowed= lambda x: f_target.space_group_info().is_allowed_origin_shift(
    x, tolerance=0.02)
  i=0
  while 1:
    peak = correlation_map_peak_search.next()
    result.sharp_correlation_map = peak.height >= 0.9
    if not result.sharp_correlation_map:
      if verbose:
        print "** weak correlation map highest peak **"
      return result
    shift = mat.col(peak.site)
    if (   is_allowed(shift - reference_shift)
        or is_allowed(shift + reference_shift)):
      result.first_correct_correlation_peak = i
      break
    i += 1
  if verbose and result.first_correct_correlation_peak != 0:
    print "First correct correlation peak: #%i"\
          % result.first_correct_correlation_peak

  # Shift the structure back accordingly
  multiplier = 1000
  shift *= multiplier
  cb_op = sgtbx.change_of_basis_op(
    sgtbx.rt_mx(sgtbx.tr_vec(shift.as_int(), multiplier)))
  f_calc = flipping.f_calc
  f_calc = f_calc.change_basis(cb_op)
  f_calc = f_calc.customized_copy(space_group_info=f_obs.space_group_info())\
                 .merge_equivalents().array()
  from crys3d import wx_map_viewer
  isinstance(f_calc, miller.array)
  wx_map_viewer.display(f_calc.fft_map())

  dphi = flex.fmod_positive(f_calc.phases().data() - f_target.phases().data(),
                            2*math.pi)


  # return a bunch of Boolean flags telling where it failed
  # or whether it succeeded
  return result

def exercise(flags, space_group_info, flipping_type):
  results = []
  if flags.Verbose:
    print flipping_type.__name__
  n = len(space_group_info.group())
  n_C = 12//n or 1
  n_O = 6//n
  n_N = 3//n
  if flags.Verbose:
    print "C%i O%i N%i" % (n_C*n, n_O*n, n_N*n)
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

def run():
  import sys
  results = []
  results.extend(
    debug_utils.parse_options_loop_space_groups(
      sys.argv[1:],
      exercise,
      #flipping_type=charge_flipping.weak_reflection_improved_iterator,
      flipping_type=charge_flipping.basic_iterator,
    ))
  results = flat_list(results)
  n_tests = len(results)

  n_phase_transitions = [
    r.had_phase_transition for r in results ].count(True)
  n_emma_matches = [
    r.emma_match for r in results ].count(True)
  n_sharp_correlation_map = [
    r.sharp_correlation_map for r in results ].count(True)
  n_found_shift = [
    r.first_correct_correlation_peak is not None for r in results ].count(True)
  print "%i test cases:" % n_tests
  print "\t%i phase transitions" % n_phase_transitions
  print "\t%i Euclidean matches with correct structure" % n_emma_matches
  print "\t%i sharp correlation maps" % n_sharp_correlation_map
  print "\t%i found shifts" % n_found_shift
  assert n_found_shift/n_tests > 2/3

if __name__ == '__main__':
  run()
