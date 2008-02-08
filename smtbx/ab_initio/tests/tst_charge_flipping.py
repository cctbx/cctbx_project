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
from libtbx import easy_pickle
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

  f_indices = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=anomalous_flag,
    d_min=d_min)
  f_target = f_indices.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc()

  f_000 = flex.miller_index(1)
  f_000[0] = (0,0,0)
  f_000 = miller.set(structure, f_000).structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc().data()[0]
  return f_target, f_000, structure


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
                      d_min=0.8, grid_resolution_factor=1./2,
                      verbose=False
                      ):
  result = group_args(had_phase_transition = False,
                      emma_match = False)

  f_target, f_000, target_structure = random_structure_factors(
    space_group_info, elements,
    anomalous_flag,
    d_min, grid_resolution_factor)
  f_target_in_p1 = f_target.expand_to_p1()\
                           .as_non_anomalous_array()\
                           .merge_equivalents().array()
  f_obs = abs(f_target)

  # Guessing a value of delta leading to subsequent good convergence
  if verbose:
    print "Guessing delta..."
    print "%10s | %10s | %10s | %10s | %10s" % ('delta', 'R', 'c_tot', 'c_flip',
                                        'c_tot/c_flip')
    print "-"*64
  delta = None
  flipping = flipping_type(f_obs, delta)
  while 1:
    for i in xrange(10): flipping.next()
    if verbose:
      rho = flipping.rho_map
      c_tot = rho.c_tot()
      c_flip = rho.c_flip(flipping.delta)
      # to compare with superflip output
      c_tot *= flipping.fft_scale; c_flip *= flipping.fft_scale
      print "%10.4f | %10.3f | %10.1f | %10.1f | %10.2f"\
            % (flipping.delta, flipping.r1_factor(),
               c_tot, c_flip, c_tot/c_flip)
    flipping.adjust_delta()
    if flipping.delta != flipping.old_delta:
      flipping.restart()
    else:
      break

  # main charge flipping loop to solve the structure
  if verbose:
    print
    print "Solving..."
    print "with delta=%.4f" % flipping.delta
    print
    print "%5s | %10s | %10s" % ('#', 'R1', 'c_tot/c_flip')
    print '-'*33
  r1_s = flex.double()
  c_tot_over_c_flip_s = flex.double()
  for i,state in itertbx.islice(enumerate(flipping), 0, 400, 10):
    r1 = state.r1_factor()
    r = state.c_tot_over_c_flip()
    r1_s.append(r1)
    c_tot_over_c_flip_s.append(r)
    if verbose:
      print "%5i | %10.3f | %10.3f" % (i, r1, r)
  f_result_in_p1 = state.f_calc

  # check whether a phase transition has occured
  diff_ctot_over_cflip = c_tot_over_c_flip_s[1:] - c_tot_over_c_flip_s[:-1]
  diff_r1 = r1_s[1:] - r1_s[:-1]
  i,j = flex.min_index(diff_ctot_over_cflip), flex.min_index(diff_r1)
  if 0 <= abs(i-j) <= 3:
    for k, diff in ((i,diff_ctot_over_cflip), (j,diff_r1)):
      normalised = diff/flex.min(diff)
      k = k+2
      small = (flex.abs(normalised[k:]) < 0.1).count(True)
      small /= len(normalised) - k
      if small < 0.8: break
    else:
      result.had_phase_transition = True
  if not result.had_phase_transition:
    if verbose: print "** no phase transition **"
    return result

  # sharpen the map
  polishing = charge_flipping.low_density_elimination_iterator(
    f_obs, f_calc=flipping.f_calc, f_000=0, rho_c=flipping.rho_c)
  for i,state in enumerate(itertbx.islice(polishing, 5)): pass

  # Euclidean matching of the peaks from the obtained map
  # against those of the correct structure
  target_structure_in_p1 = target_structure.expand_to_p1()
  search_parameters = maptbx.peak_search_parameters(
    interpolate=True,
    min_distance_sym_equiv=1.,
    max_clusters=int(target_structure_in_p1.scatterers().size()*1.05))
  peak_search_outcome = polishing.rho_map.peak_search(search_parameters)
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
  if not result.emma_match and verbose:
    print "** no Euclidean matching **"
  reference_shift = refined_matches[0].rt.t
  inverted_solution = refined_matches[0].rt.r == mat.identity(3)

  # Find the translation to bring back the structure to the same space-group
  # setting as the starting f_obs from correlation map analysis
  if inverted_solution:
    reference_shift = -reference_shift
  correlation_map_peak_search = polishing.search_origin()
  highest_peak = correlation_map_peak_search.next()
  assert highest_peak.height >= 0.99
  shift = -mat.col(highest_peak.site)
  delta = shift - reference_shift
  assert f_target.space_group_info().is_allowed_origin_shift(
    delta.elems, tolerance=0.04)

  # Shift the structure back accordingly


  # return a bunch of Boolean flags telling where it failed
  # or whether it succeeded
  return result

def exercise(flags, space_group_info, flipping_type):
  results = []
  if flags.Verbose:
    print flipping_type.__name__
  for i in xrange(2):
    n_C = random.randint(4,6)
    n_O = random.randint(1,2)
    n_N = random.randint(1,2)
    if flags.Verbose:
      print "C%i O%i N%i" % (n_C, n_O, n_N)
    result = randomly_exercise(
      flipping_type=flipping_type,
      space_group_info=space_group_info,
      elements=["C"]*n_C + ["O"]*n_O + ["N"]*n_N,
      anomalous_flag=False,
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
      flipping_type=charge_flipping.basic_iterator))
  results = flat_list(results)
  n_tests = len(results)
  n_phase_transitions = len(filter(lambda r: r.had_phase_transition, results))
  n_emma_matches = len(filter(lambda r: r.emma_match, results))
  print "%i test cases:" % n_tests
  print "\t%i phase transitions" % n_phase_transitions
  print "\t%i Euclidean matches with correct structure" % n_emma_matches
  assert n_emma_matches/n_tests > 0.8

if __name__ == '__main__':
  import sys
  if '--profile' in sys.argv:
    import profile
    import pstats
    sys.argv.remove('--profile')
    profile.run('run()', 'charge_flipping.prof')
    p = pstats.Stats('charge_flipping.prof')
    p.strip_dirs().sort_stats('time').print_stats(10)
  else:
    run()
