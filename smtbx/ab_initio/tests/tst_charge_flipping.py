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
                      emma_match_in_p1=False,
                      emma_match=False,
                      first_correct_correlation_peak=None,
                      inverted_solution=False,
                      succeeded=True)

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
  solving = charge_flipping.solving_iterator(
    flipping,
    yield_during_delta_guessing=True,
    yield_solving_interval=1)
  charge_flipping.loop(solving, verbose=verbose)

  # check whether a phase transition has occured
  result.had_phase_transition = (solving.state == solving.finished
                                 and not solving.max_attempts_exceeded)
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
    break_if_match_with_no_singles=True
    ).refined_matches
  result.emma_match_in_p1 = (
        not refined_matches[0].singles1 # all sites match a peak
    and refined_matches[0].rms < 0.1 # no farther than that
    and refined_matches[0].rt.r in (mat.identity(3), mat.inversion(3))
  )
  if not result.emma_match_in_p1:
    result.succeeded = False
    if verbose:
      print "** no Euclidean matching in P1 **"
    if verbose == "debug":
      print refined_matches[0].show()
  result.inverted_solution = refined_matches[0].rt.r == mat.inversion(3)
  reference_shift = -refined_matches[0].rt.t

  if not result.succeeded:
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

  if not result.succeeded:
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
    break_if_match_with_no_singles=True
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
  if result.succeeded:
    if result.first_correct_correlation_peak > 0:
      print "@@ Success (but not first solution) @@"
    else:
      print "@@ Success @@"
  else:
    print "@@ Failure @@"
  return result

def exercise(flags, space_group_info, flipping_type):
  results = []
  if flags.Verbose:
    print flipping_type.__name__
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
  loops = 1
  while loops > 0:
    result = randomly_exercise(
      flipping_type=flipping_type,
      space_group_info=space_group_info,
      elements=["C"]*n_C + ["O"]*n_O + ["N"]*n_N,
      anomalous_flag=False,
      d_min=0.8,
      verbose=flags.Verbose
    )
    results.append(result)
    if result.succeeded: loops -= 1
    else: loops += 2
  if flags.Verbose: print
  return True, results

def exercise_observable_evolution():
  values = [0.9243433862036855, 0.93888446148759841, 0.89297272756934098,
            0.91444304708908963, 0.89056379497301141, 0.88120944583371341,
            0.88225452586909803, 0.88653738668340243, 0.87070792588572521,
            0.87722692361678156, 0.88205556705027788, 0.85829946112348665,
            0.86698377695745421, 0.86864740321488576, 0.85487075022849845,
            0.86865609367547991, 0.86964112171406505, 0.85298264958450531,
            0.85385678312448687, 0.84464652057874379, 0.84641043987966835,
            0.83012857732351208, 0.82799936217038572, 0.83383947929782087,
            0.82100764553967154, 0.7915023822419347, 0.81317081214522768,
            0.81413024795359601, 0.80884067976073648, 0.79052504162717585,
            0.78985833674796491, 0.77799932305727437, 0.76101206532104115,
            0.77112506918389512, 0.73156878058867747, 0.77114312453174139,
            0.74230962083507179, 0.74943686259593456, 0.71967954299605441,
            0.71664377718493688, 0.67537908119607493, 0.68777789151315794,
            0.64678822971049599, 0.62324787899995138, 0.60309076681633222,
            0.5938063578971694, 0.55752209214443882, 0.53820599976390726,
            0.51571466469843452, 0.50227892570081023, 0.50006551090621876,
            0.48719104307110472, 0.47592932785238595, 0.47741992795222854,
            0.46191378819812401, 0.46437237548136973, 0.460152406059613,
            0.45696425621502484, 0.462905034445043, 0.45636278581522099,
            0.44523652675860226, 0.45265613293682072, 0.4470376817756041,
            0.45850692904946844, 0.433276726842306, 0.45188696158141439,
            0.44231654773228324, 0.43936464662363689, 0.43875493752518174,
            0.44542996478870139, 0.44392154956012092, 0.43546542725905513,
            0.44320196590966798, 0.4349473564169557, 0.44600231315360889,
            0.43749411002738853, 0.44119521859852201, 0.43648969576173563,
            0.44346674738070124, 0.43547858166530307, 0.44147312503457359,
            0.43606290339204423]
  obs = charge_flipping.observable_evolution(phase_transition_tail_len=12)
  for x in values: obs.append(x)
  assert obs.had_phase_transition()

  values = [0.2554279822789518, 0.25175775041554682, 0.25282900644628553,
            0.24956755929030103, 0.25134804082115031, 0.24823326310752472,
            0.25361460738591829, 0.25142091078242934, 0.25030717051497342,
            0.25284616018151584, 0.24950546799578305, 0.24558359248827596,
            0.25289591331275368, 0.24886195298568356, 0.24809057445962621,
            0.25053888372266298, 0.25168433950991959, 0.2484213864128641,
            0.24808634182911996, 0.24742500272491222, 0.25057819324558295,
            0.24888637679400497, 0.24695894935433052, 0.24391879976546502,
            0.24906807893674035, 0.24697315953853607, 0.246667742607214,
            0.2449340548413029, 0.25174581564122372, 0.24767715308474855,
            0.24720416630900058, 0.24684040315022765, 0.25265854389702563,
            0.24864591853855753, 0.25080291250497444, 0.24685069605912452,
            0.24371980947637228, 0.24552468030017974, 0.25123700327698323,
            0.24700853980274171, 0.24654863199905797, 0.24411686665334578,
            0.24752140679832801, 0.24457394837128227, 0.25126992296401723,
            0.24573975526583799, 0.24798651164559868, 0.24670568511735705,
            0.24703648010484003, 0.24554818654316407, 0.24698505996580195,
            0.24729504352609494, 0.24514172689583177, 0.24567600042439633,
            0.2467318184230852, 0.24550175344675165, 0.24583589472765699,
            0.2459125367412881, 0.25227023343741589, 0.2467613760878056,
            0.24733729418516187, 0.24979539453857372, 0.24760420633372374,
            0.24272444272902446, 0.24931003345777736, 0.2453974299383018,
            0.2449879877142789, 0.24325920237793267, 0.24801912339496765,
            0.2451545474330882, 0.24797624604665652, 0.24585678438635494,
            0.2457269774883338, 0.24282154299582528, 0.24414702252539228,
            0.24724842044537854, 0.24879139941023312, 0.24491862536506112,
            0.2448607320349229, 0.24257412892564048, 0.24359058238994219,
            0.24631369435027142]
  obs = charge_flipping.observable_evolution(phase_transition_tail_len=12)
  for x in values: obs.append(x)
  assert not obs.had_phase_transition()

def exercise_charge_flipping():
  import sys
  results = []
  results.extend(
    debug_utils.parse_options_loop_space_groups(
      sys.argv[1:],
      exercise,
      flipping_type=charge_flipping.weak_reflection_improved_iterator,
      #flipping_type=charge_flipping.basic_iterator,
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
  assert n_success/n_tests > 2/3

def run():
  exercise_observable_evolution()
  exercise_charge_flipping()
  print format_cpu_times()

if __name__ == '__main__':
  run()
