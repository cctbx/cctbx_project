from cctbx import dmtbx
from cctbx import maptbx
from cctbx import miller
from scitbx import fftpack
from scitbx.python_utils import complex_math
from cctbx.utils import phase_error
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
import random
import math
import sys

def direct_space_squaring(start, selection_fixed):
  map_gridding = miller.index_span(
    miller.set.expand_to_p1(start).indices()).map_grid()
  if (selection_fixed is None):
    fixed = start
    var = start
  else:
    fixed = start.apply_selection(selection_fixed)
    var = start.apply_selection(~selection_fixed)
  rfft = fftpack.real_to_complex_3d([n*3//2 for n in map_gridding])
  conjugate_flag = 0001
  structure_factor_map = maptbx.structure_factors.to_map(
    fixed.space_group(),
    fixed.anomalous_flag(),
    fixed.indices(),
    fixed.data(),
    rfft.n_real(),
    flex.grid(rfft.n_complex()),
    conjugate_flag)
  real_map = rfft.backward(structure_factor_map.complex_map())
  squared_map = flex.pow2(real_map)
  squared_sf_map = rfft.forward(squared_map)
  allow_miller_indices_outside_map = 00000
  from_map = maptbx.structure_factors.from_map(
    var.anomalous_flag(),
    var.indices(),
    squared_sf_map,
    conjugate_flag,
    allow_miller_indices_outside_map)
  if (selection_fixed is None):
    return from_map.data()
  result = start.data().deep_copy()
  result.set_selected(~selection_fixed, from_map.data())
  assert result.select(selection_fixed).all_eq(fixed.data())
  return result

def reciprocal_space_squaring(start, selection_fixed, verbose):
  tprs = dmtbx.triplet_generator(miller_set=start)
  if (0 or verbose):
    for ih in start.indices()[:1].indices():
      for relation in tprs.relations_for(ih):
        print relation.format(start.indices(), ih),
        if (not relation.is_sigma_2(ih)):
          print "not sigma-2",
        print
  amplitudes = abs(start).data()
  if (selection_fixed is not None):
    sel_var = ~selection_fixed
    amplitudes.set_selected(sel_var, flex.double(sel_var.count(0001), 0))
  input_phases = flex.arg(start.data())
  result = tprs.apply_tangent_formula(
    amplitudes=amplitudes,
    phases=input_phases,
    selection_fixed=selection_fixed,
    use_fixed_only=selection_fixed is not None)
  if (selection_fixed is not None):
    assert result.select(selection_fixed).all_eq(
      input_phases.select(selection_fixed))
  return result

def exercise(space_group_info, n_scatterers=8, d_min=2, verbose=0,
             e_min=1.5):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=200,
    min_distance=3.,
    general_positions_only=0001,
    u_iso=0.0)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=00000).f_calc()
  f_obs = abs(f_calc)
  q_obs = miller.array(
    miller_set=f_obs,
    data=f_obs.data()
        / math.sqrt(f_obs.space_group().order_p() * n_scatterers)
        / f_obs.space_group().n_ltr())
  q_obs = q_obs.sort(by_value="abs")
  q_obs.setup_binner(auto_binning=True)
  n_obs = q_obs.normalize_structure_factors(quasi=0001)
  r = flex.linear_regression(q_obs.data(), n_obs.data())
  if (0 or verbose):
    r.show_summary()
  assert r.is_well_defined()
  assert abs(r.y_intercept()) < 0.2
  assert abs(r.slope() - 1) < 0.3
  q_large = q_obs.apply_selection(
    q_obs.quasi_normalized_as_normalized().data() > e_min)
  if (0 or verbose):
    print "Number of e-values > %.6g: %d" % (e_min, q_large.size())
  other_structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=200,
    min_distance=3.,
    general_positions_only=0001,
    u_iso=0.0)
  assert other_structure.unit_cell().is_similar_to(structure.unit_cell())
  q_calc = q_large.structure_factors_from_scatterers(
    other_structure, direct=0001).f_calc()
  start = q_large.phase_transfer(q_calc.data())
  for selection_fixed in (
        None,
        flex.double([random.random() for i in xrange(start.size())]) < 0.4):
    from_map_data = direct_space_squaring(start, selection_fixed)
    direct_space_result = start.phase_transfer(phase_source=from_map_data)
    new_phases = reciprocal_space_squaring(start, selection_fixed, verbose)
    reciprocal_space_result = start.phase_transfer(
      phase_source=flex.polar(1,new_phases))
    mwpe = direct_space_result.mean_weighted_phase_error(
      reciprocal_space_result)
    if (0 or verbose):
      print "mwpe: %.2f" % mwpe, start.space_group_info()
    for i,h in direct_space_result.indices().items():
      amp_d,phi_d = complex_math.abs_arg(
        direct_space_result.data()[i], deg=1)
      amp_r,phi_r = complex_math.abs_arg(
        reciprocal_space_result.data()[i],deg=1)
      phase_err = phase_error(phi_d, phi_r, deg=1)
      assert phase_err < 1.0 or abs(from_map_data[i]) < 1.e-6

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
