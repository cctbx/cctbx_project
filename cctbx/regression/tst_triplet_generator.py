from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from cctbx import dmtbx
from cctbx import maptbx
from cctbx import miller
from scitbx import fftpack
from libtbx import complex_math
import scitbx.math
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx.development import debug_utils
import random
import math
import sys
from six.moves import range

def direct_space_squaring(start, selection_fixed):
  map_gridding = miller.index_span(
    miller.set.expand_to_p1(start).indices()).map_grid()
  if (selection_fixed is None):
    fixed = start
    var = start
  else:
    fixed = start.select(selection_fixed)
    var = start.select(~selection_fixed)
  rfft = fftpack.real_to_complex_3d([n*3//2 for n in map_gridding])
  conjugate_flag = True
  structure_factor_map = maptbx.structure_factors.to_map(
    space_group=fixed.space_group(),
    anomalous_flag=fixed.anomalous_flag(),
    miller_indices=fixed.indices(),
    structure_factors=fixed.data(),
    n_real=rfft.n_real(),
    map_grid=flex.grid(rfft.n_complex()),
    conjugate_flag=conjugate_flag)
  real_map = rfft.backward(structure_factor_map.complex_map())
  squared_map = flex.pow2(real_map)
  squared_sf_map = rfft.forward(squared_map)
  allow_miller_indices_outside_map = False
  from_map = maptbx.structure_factors.from_map(
    anomalous_flag=var.anomalous_flag(),
    miller_indices=var.indices(),
    complex_map=squared_sf_map,
    conjugate_flag=conjugate_flag,
    allow_miller_indices_outside_map=allow_miller_indices_outside_map)
  if (selection_fixed is None):
    return from_map.data()
  result = start.data().deep_copy()
  result.set_selected(~selection_fixed, from_map.data())
  assert result.select(selection_fixed).all_eq(fixed.data())
  return result

def reciprocal_space_squaring(start, selection_fixed, verbose):
  tprs = dmtbx.triplet_generator(miller_set=start)
  if (0 or verbose):
    for ih in range(start.indices()[:1].size()):
      for relation in tprs.relations_for(ih):
        print(relation.format(start.indices(), ih), end=' ')
        if (not relation.is_sigma_2(ih)):
          print("not sigma-2", end=' ')
        print()
  amplitudes = abs(start).data()
  if (selection_fixed is not None):
    amplitudes.set_selected(~selection_fixed, 0)
  input_phases = flex.arg(start.data())
  result = tprs.apply_tangent_formula(
    amplitudes=amplitudes,
    phases_rad=input_phases,
    selection_fixed=selection_fixed,
    use_fixed_only=selection_fixed is not None)
  if (selection_fixed is not None):
    assert result.select(selection_fixed).all_eq(
      input_phases.select(selection_fixed))
  return result

def exercise_truncate(q_large):
  tprs_full = dmtbx.triplet_generator(
    miller_set=q_large,
    discard_weights=True)
  tprs = dmtbx.triplet_generator(
    miller_set=q_large,
    amplitudes=q_large.data(),
    max_relations_per_reflection=0,
    discard_weights=True)
  assert tprs.n_relations().all_eq(tprs_full.n_relations())
  for n in (1,10,100,1000):
    tprs = dmtbx.triplet_generator(
      miller_set=q_large,
      amplitudes=q_large.data(),
      max_relations_per_reflection=n,
      discard_weights=True)
    assert (tprs.n_relations() >= n).all_eq(tprs.n_relations() == n)
  n = 3
  tprs = dmtbx.triplet_generator(
    miller_set=q_large,
    amplitudes=q_large.data(),
    max_relations_per_reflection=n,
    discard_weights=True)
  n_rel_full = tprs_full.n_relations()
  n_rel = tprs.n_relations()
  amp = q_large.data()
  for ih in range(q_large.indices().size()):
    if (n_rel[ih] == n_rel_full[ih]): continue
    aa_full = flex.double()
    for relation in tprs_full.relations_for(ih):
      aa_full.append(amp[relation.ik()] * amp[relation.ihmk()])
    aa = flex.double()
    for relation in tprs.relations_for(ih):
      aa.append(amp[relation.ik()] * amp[relation.ihmk()])
    aa_full = aa_full.select(flex.sort_permutation(data=aa_full, reverse=True))
    assert approx_equal(aa_full[:n], aa)

def exercise(space_group_info, n_scatterers=8, d_min=2, verbose=0,
             e_min=1.5):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=200,
    min_distance=3.,
    general_positions_only=True,
    u_iso=0.0)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  f_obs = abs(f_calc)
  q_obs = miller.array(
    miller_set=f_obs,
    data=f_obs.data()
        / math.sqrt(f_obs.space_group().order_p() * n_scatterers)
        / f_obs.space_group().n_ltr())
  q_obs = q_obs.sort(by_value="abs")
  q_obs.setup_binner(auto_binning=True)
  n_obs = q_obs.quasi_normalize_structure_factors()
  r = flex.linear_regression(q_obs.data(), n_obs.data())
  if (0 or verbose):
    r.show_summary()
  assert r.is_well_defined()
  assert abs(r.y_intercept()) < 0.1
  assert abs(r.slope() - 1) < 0.2
  q_large = q_obs.select(
    q_obs.quasi_normalized_as_normalized().data() > e_min)
  if (0 or verbose):
    print("Number of e-values > %.6g: %d" % (e_min, q_large.size()))
  other_structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=200,
    min_distance=3.,
    general_positions_only=True,
    u_iso=0.0)
  assert other_structure.unit_cell().is_similar_to(structure.unit_cell())
  q_calc = q_large.structure_factors_from_scatterers(
    other_structure, algorithm="direct").f_calc()
  start = q_large.phase_transfer(q_calc.data())
  for selection_fixed in (
        None,
        flex.double([random.random() for i in range(start.size())]) < 0.4):
    from_map_data = direct_space_squaring(start, selection_fixed)
    direct_space_result = start.phase_transfer(phase_source=from_map_data)
    new_phases = reciprocal_space_squaring(start, selection_fixed, verbose)
    reciprocal_space_result = start.phase_transfer(
      phase_source=flex.polar(1,new_phases))
    mwpe = direct_space_result.mean_weighted_phase_error(
      reciprocal_space_result)
    if (0 or verbose):
      print("mwpe: %.2f" % mwpe, start.space_group_info())
    for i,h in enumerate(direct_space_result.indices()):
      amp_d,phi_d = complex_math.abs_arg(
        direct_space_result.data()[i], deg=True)
      amp_r,phi_r = complex_math.abs_arg(
        reciprocal_space_result.data()[i],deg=True)
      phase_err = scitbx.math.phase_error(phi_d, phi_r, deg=True)
      assert phase_err < 1.0 or abs(from_map_data[i]) < 1.e-6
  exercise_truncate(q_large)

def run_call_back(flags, space_group_info):
  exercise(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
