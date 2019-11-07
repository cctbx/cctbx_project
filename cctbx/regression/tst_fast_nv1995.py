from __future__ import absolute_import, division, print_function
from cctbx import translation_search
from cctbx import crystal
from cctbx import miller
from cctbx import xray
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx import matrix
import random
import sys
from six.moves import range
from six.moves import zip

def run_fast_terms(structure_fixed, structure_p1,
                   f_obs, f_calc_fixed, f_calc_p1,
                   symmetry_flags, gridding, grid_tags,
                   n_sample_grid_points=10,
                   test_origin=False,
                   verbose=0):
  if (f_calc_fixed is None):
    f_part = flex.complex_double()
  else:
    f_part = f_calc_fixed.data()
  m = flex.double()
  for i in range(f_obs.indices().size()):
    m.append(random.random())
  assert f_obs.anomalous_flag() == f_calc_p1.anomalous_flag()
  fast_terms = translation_search.fast_terms(
    gridding=gridding,
    anomalous_flag=f_obs.anomalous_flag(),
    miller_indices_p1_f_calc=f_calc_p1.indices(),
    p1_f_calc=f_calc_p1.data())
  for squared_flag in (False, True):
    map = fast_terms.summation(
      space_group=f_obs.space_group(),
      miller_indices_f_obs=f_obs.indices(),
      m=m,
      f_part=f_part,
      squared_flag=squared_flag).fft().accu_real_copy()
    assert map.all() == gridding
    map_stats = maptbx.statistics(map)
    if (0 or verbose):
      map_stats.show_summary()
    grid_tags.build(f_obs.space_group_info().type(), symmetry_flags)
    assert grid_tags.n_grid_misses() == 0
    assert grid_tags.verify(map)
    for i_sample in range(n_sample_grid_points):
      run_away_counter = 0
      while 1:
        run_away_counter += 1
        assert run_away_counter < 1000
        if (i_sample == 0 and test_origin):
          grid_point = [0,0,0]
        else:
          grid_point = [random.randrange(g) for g in gridding]
        grid_site = [float(x)/g for x,g in zip(grid_point,gridding)]
        structure_shifted = structure_fixed.deep_copy_scatterers()
        assert structure_shifted.special_position_indices().size() == 0
        structure_shifted.add_scatterers(
          scatterers=structure_p1.apply_shift(grid_site).scatterers())
        if (structure_shifted.special_position_indices().size() == 0):
          break
        if (test_origin):
          assert i_sample != 0
      i_grid = flex.norm(f_obs.structure_factors_from_scatterers(
        xray_structure=structure_shifted, algorithm="direct").f_calc().data())
      if (squared_flag): p = 4
      else: p = 2
      map_value = map[grid_point] * f_obs.space_group().n_ltr()**p
      if (not squared_flag):
        sum_m_i_grid = flex.sum(m * i_grid)
      else:
        sum_m_i_grid = flex.sum(m * flex.pow2(i_grid))
      assert "%.6g" % sum_m_i_grid == "%.6g" % map_value, (
        sum_m_i_grid, map_value)

def run_fast_nv1995(f_obs, f_calc_fixed, f_calc_p1,
                    symmetry_flags, gridding, grid_tags, verbose):
  if (f_calc_fixed is None):
    f_part = flex.complex_double()
  else:
    f_part = f_calc_fixed.data()
  assert f_obs.anomalous_flag() == f_calc_p1.anomalous_flag()
  fast_nv1995 = translation_search.fast_nv1995(
    gridding=gridding,
    space_group=f_obs.space_group(),
    anomalous_flag=f_obs.anomalous_flag(),
    miller_indices_f_obs=f_obs.indices(),
    f_obs=f_obs.data(),
    f_part=f_part,
    miller_indices_p1_f_calc=f_calc_p1.indices(),
    p1_f_calc=f_calc_p1.data())
  assert fast_nv1995.target_map().all() == gridding
  map_stats = maptbx.statistics(fast_nv1995.target_map())
  if (0 or verbose):
    map_stats.show_summary()
  grid_tags.build(f_obs.space_group_info().type(), symmetry_flags)
  assert grid_tags.n_grid_misses() == 0
  assert grid_tags.verify(fast_nv1995.target_map())
  peak_list = maptbx.peak_list(
    data=fast_nv1995.target_map(),
    tags=grid_tags.tag_array(),
    peak_search_level=1,
    max_peaks=10,
    interpolate=True)
  if (0 or verbose):
    print("gridding:", gridding)
    for i,site in enumerate(peak_list.sites()):
      print("(%.4f,%.4f,%.4f)" % site, "%.6g" % peak_list.heights()[i])
  assert approx_equal(map_stats.max(), flex.max(peak_list.grid_heights()))
  return peak_list

def test_atom(space_group_info, use_primitive_setting,
              n_elements=3, d_min=3.,
              grid_resolution_factor=0.48, max_prime=5, verbose=0):
  if (use_primitive_setting):
    space_group_info = space_group_info.primitive_setting()
  structure = random_structure.xray_structure(
    space_group_info,
    n_scatterers=n_elements,
    volume_per_atom=150,
    min_distance=1.,
    general_positions_only=True)
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=(random.random() < 0.5),
    d_min=d_min)
  symmetry_flags = translation_search.symmetry_flags(
    is_isotropic_search_model=True,
    have_f_part=(n_elements>=2))
  gridding = miller_set_f_obs.crystal_gridding(
    symmetry_flags=symmetry_flags,
    resolution_factor=grid_resolution_factor,
    max_prime=max_prime).n_real()
  structure.build_scatterers(
    elements=["Se"]*n_elements,
    grid=gridding)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_obs = abs(miller_set_f_obs.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc())
  if (0 or verbose):
    f_obs.show_summary()
  if (0 or verbose):
    f_obs.show_array()
  miller_set_p1 = miller.set.expand_to_p1(f_obs)
  special_position_settings_p1 = crystal.special_position_settings(
    crystal_symmetry=miller_set_p1)
  structure_fixed = xray.structure(special_position_settings=structure)
  for scatterer in structure.scatterers():
    structure_p1 = xray.structure(
      special_position_settings=special_position_settings_p1)
    scatterer_at_origin = scatterer.customized_copy(site=(0,0,0))
    structure_p1.add_scatterer(scatterer_at_origin)
    if (0 or verbose):
      structure_p1.show_summary().show_scatterers()
    f_calc_p1 = miller_set_p1.structure_factors_from_scatterers(
      xray_structure=structure_p1,
      algorithm="direct").f_calc()
    if (0 or verbose):
      f_calc_p1.show_array()
    f_calc_fixed = None
    if (structure_fixed.scatterers().size() > 0):
      f_calc_fixed = f_obs.structure_factors_from_scatterers(
        xray_structure=structure_fixed,
        algorithm="direct").f_calc()
    symmetry_flags = translation_search.symmetry_flags(
      is_isotropic_search_model=True,
      have_f_part=(f_calc_fixed is not None))
    if (structure_fixed.scatterers().size() <= 1):
      gridding = miller_set_f_obs.crystal_gridding(
        symmetry_flags=symmetry_flags,
        resolution_factor=grid_resolution_factor,
        max_prime=max_prime).n_real()
      grid_tags = maptbx.grid_tags(gridding)
    run_fast_terms(
      structure_fixed, structure_p1,
      f_obs, f_calc_fixed, f_calc_p1,
      symmetry_flags, gridding, grid_tags,
      verbose=verbose)
    peak_list = run_fast_nv1995(
      f_obs, f_calc_fixed, f_calc_p1,
      symmetry_flags, gridding, grid_tags, verbose)
    structure_fixed.add_scatterer(scatterer)
    if (0 or verbose):
      structure_fixed.show_summary().show_scatterers()
    if (structure_fixed.scatterers().size() < n_elements):
      assert peak_list.heights()[0] < 1
    else:
      assert peak_list.heights()[0] > 0.99
  assert peak_list.heights()[0] > 0.99

def test_molecule(space_group_info, use_primitive_setting, flag_f_part,
                  d_min=3., grid_resolution_factor=0.48, max_prime=5,
                  verbose=0):
  if (use_primitive_setting):
    space_group_info = space_group_info.primitive_setting()
  elements = ("N", "C", "C", "O", "N", "C", "C", "O")
  structure = random_structure.xray_structure(
    space_group_info,
    elements=elements,
    volume_per_atom=50,
    min_distance=1.,
    general_positions_only=True,
    random_u_iso=True,
    random_occupancy=True)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=(random.random() < 0.5),
    d_min=d_min)
  f_obs = abs(miller_set_f_obs.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc())
  if (0 or verbose):
    f_obs.show_summary()
  if (0 or verbose):
    f_obs.show_array()
  miller_set_p1 = miller.set.expand_to_p1(f_obs)
  special_position_settings_p1 = crystal.special_position_settings(
    crystal_symmetry=miller_set_p1)
  structure_p1 = xray.structure(
    special_position_settings=special_position_settings_p1)
  structure_fixed = xray.structure(special_position_settings=structure)
  for scatterer in structure.scatterers():
    if (flag_f_part and   structure_fixed.scatterers().size()
                        < structure.scatterers().size()//2):
      structure_fixed.add_scatterer(scatterer)
    else:
      structure_p1.add_scatterer(scatterer)
  if (0 or verbose):
    if (flag_f_part):
      structure_fixed.show_summary().show_scatterers()
    structure_p1.show_summary().show_scatterers()
  f_calc_fixed = None
  if (flag_f_part):
    f_calc_fixed = f_obs.structure_factors_from_scatterers(
      xray_structure=structure_fixed,
      algorithm="direct").f_calc()
  f_calc_p1 = miller_set_p1.structure_factors_from_scatterers(
    xray_structure=structure_p1,
    algorithm="direct").f_calc()
  symmetry_flags = translation_search.symmetry_flags(
    is_isotropic_search_model=False,
    have_f_part=flag_f_part)
  gridding = miller_set_f_obs.crystal_gridding(
    symmetry_flags=symmetry_flags,
    resolution_factor=grid_resolution_factor,
    max_prime=max_prime).n_real()
  grid_tags = maptbx.grid_tags(gridding)
  run_fast_terms(
    structure_fixed, structure_p1,
    f_obs, f_calc_fixed, f_calc_p1,
    symmetry_flags, gridding, grid_tags,
    test_origin=True,
    verbose=verbose)
  peak_list = run_fast_nv1995(
    f_obs, f_calc_fixed, f_calc_p1,
    symmetry_flags, gridding, grid_tags, verbose)
  assert peak_list.heights()[0] > 0.99

def test_shift(space_group_info,
               d_min=0.8, grid_resolution_factor=0.48, max_prime=5,
               verbose=0):
  n = 12 // len(space_group_info.group()) or 1
  target_structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=['C']*n,
    use_u_iso=False,
    use_u_aniso=False,
  )
  f_target = miller.build_set(
    crystal_symmetry=target_structure,
    anomalous_flag=False,
    d_min=d_min
  ).structure_factors_from_scatterers(
    xray_structure=target_structure,
    algorithm="direct").f_calc()
  f_obs = abs(f_target)

  indices_in_p1 = miller.set.expand_to_p1(f_target)
  target_structure_in_p1 = target_structure.expand_to_p1()

  reference_translation = matrix.col((0.1, 0.2, 0.7))
  structure_in_p1 = target_structure_in_p1.apply_shift(reference_translation)
  f_structure_in_p1 = indices_in_p1.structure_factors_from_scatterers(
    xray_structure=structure_in_p1,
    algorithm="direct").f_calc()

  symmetry_flags = translation_search.symmetry_flags(
    is_isotropic_search_model=False,
    have_f_part=False)
  gridding = f_target.crystal_gridding(
    symmetry_flags=symmetry_flags,
    resolution_factor=grid_resolution_factor,
    max_prime=max_prime).n_real()
  grid_tags = maptbx.grid_tags(gridding)
  for f_calc_in_p1 in (f_structure_in_p1,):
    peak_list = run_fast_nv1995(
      f_obs=f_obs, f_calc_fixed=None, f_calc_p1=f_calc_in_p1,
      symmetry_flags=symmetry_flags, gridding=gridding,
      grid_tags=grid_tags, verbose=verbose)
    assert peak_list.heights()[0] > 0.9
    shift = matrix.col(peak_list.sites()[0])
    assert f_target.space_group_info().is_allowed_origin_shift(
      shift + reference_translation, tolerance=0.04)

def run_call_back(flags, space_group_info):
  if (space_group_info.group().order_p() > 24 and not flags.HighSymmetry):
    print("High symmetry space group skipped.")
    return
  if (not (flags.Atom or flags.Molecule or flags.Shift)):
    flags.Atom = True
    flags.Molecule = True
  use_primitive_setting_flags = [False]
  if (space_group_info.group().conventional_centring_type_symbol() != "P"):
    use_primitive_setting_flags.append(True)
  if (flags.Atom):
    for use_primitive_setting in use_primitive_setting_flags:
      test_atom(space_group_info, use_primitive_setting,
                verbose=flags.Verbose)
  if (flags.Molecule):
    for flag_f_part in (False, True)[:]: #SWITCH
      for use_primitive_setting in use_primitive_setting_flags:
        test_molecule(space_group_info, use_primitive_setting, flag_f_part,
                      verbose=flags.Verbose)
  if flags.Shift:
    for i in range(1):
      test_shift(space_group_info, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "HighSymmetry",
    "Atom",
    "Molecule",
    "Shift"))

if (__name__ == "__main__"):
  run()
