from cctbx import translation_search
from cctbx import crystal
from cctbx import miller
from cctbx import xray
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal
import sys

def run_fast_nv1995(f_obs, f_calc_fixed, f_calc_p1,
                    symmetry_flags, gridding, grid_tags, verbose):
  if (f_calc_fixed == None):
    f_part = flex.complex_double()
  else:
    f_part = f_calc_fixed.data()
  fast_nv1995 = translation_search.fast_nv1995(
    gridding=gridding,
    space_group=f_obs.space_group(),
    anomalous_flag=00000,
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
    max_peaks=10)
  if (0 or verbose):
    print "gridding:", gridding
    for entry in peak_list.entries():
      print "(%d,%d,%d)" % entry.index, "%.6g" % (entry.value,)
  assert approx_equal(map_stats.max(), peak_list.entries()[0].value)
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
    min_distance=2.,
    general_positions_only=0001)
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=00000,
    d_min=d_min)
  symmetry_flags = translation_search.symmetry_flags(
    is_isotropic_search_model=0001,
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
    direct=0001).f_calc())
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
    scatterer_at_origin = scatterer.copy()
    scatterer_at_origin.site = (0,0,0)
    structure_p1.add_scatterer(scatterer_at_origin)
    if (0 or verbose):
      structure_p1.show_summary().show_scatterers()
    f_calc_p1 = miller_set_p1.structure_factors_from_scatterers(
      xray_structure=structure_p1,
      direct=0001).f_calc()
    if (0 or verbose):
      f_calc_p1.show_array()
    f_calc_fixed = None
    if (structure_fixed.scatterers().size() > 0):
      f_calc_fixed = f_obs.structure_factors_from_scatterers(
        xray_structure=structure_fixed,
        direct=0001).f_calc()
    symmetry_flags = translation_search.symmetry_flags(
      is_isotropic_search_model=0001,
      have_f_part=(f_calc_fixed != None))
    if (structure_fixed.scatterers().size() <= 1):
      gridding = miller_set_f_obs.crystal_gridding(
        symmetry_flags=symmetry_flags,
        resolution_factor=grid_resolution_factor,
        max_prime=max_prime).n_real()
      grid_tags = maptbx.grid_tags(gridding)
    structure_fixed.add_scatterer(scatterer)
    if (0 or verbose):
      structure_fixed.show_summary().show_scatterers()
    peak_list = run_fast_nv1995(
      f_obs, f_calc_fixed, f_calc_p1,
      symmetry_flags, gridding, grid_tags, verbose)
    if (structure_fixed.scatterers().size() < n_elements):
      assert peak_list.entries()[0].value < 1
    else:
      assert peak_list.entries()[0].value > 0.99
  assert peak_list.entries()[0].value > 0.99

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
    min_distance=2.,
    general_positions_only=0001,
    random_u_iso=0001,
    random_occupancy=0001)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=00000,
    d_min=d_min)
  f_obs = abs(miller_set_f_obs.structure_factors_from_scatterers(
    xray_structure=structure,
    direct=0001).f_calc())
  if (0 or verbose):
    f_obs.show_summary()
  if (0 or verbose):
    f_obs.show_array()
  miller_set_p1 = miller.set.expand_to_p1(f_obs)
  special_position_settings_p1 = crystal.special_position_settings(
    crystal_symmetry=miller_set_p1)
  structure_p1 = xray.structure(
    special_position_settings=special_position_settings_p1)
  if (flag_f_part):
    structure_fixed = xray.structure(
      special_position_settings=structure)
  for scatterer in structure.scatterers():
    if (flag_f_part and   structure_fixed.scatterers().size()
                        < structure.scatterers().size()/2):
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
      direct=0001).f_calc()
  f_calc_p1 = miller_set_p1.structure_factors_from_scatterers(
    xray_structure=structure_p1,
    direct=0001).f_calc()
  symmetry_flags = translation_search.symmetry_flags(
    is_isotropic_search_model=00000,
    have_f_part=flag_f_part)
  gridding = miller_set_f_obs.crystal_gridding(
    symmetry_flags=symmetry_flags,
    resolution_factor=grid_resolution_factor,
    max_prime=max_prime).n_real()
  grid_tags = maptbx.grid_tags(gridding)
  peak_list = run_fast_nv1995(
    f_obs, f_calc_fixed, f_calc_p1,
    symmetry_flags, gridding, grid_tags, verbose)
  assert peak_list.entries()[0].value > 0.99

def run_call_back(flags, space_group_info):
  if (space_group_info.group().order_p() > 8 and not flags.HighSymmetry):
    print "High symmetry space group skipped."
    return
  if (not (flags.Atom or flags.Molecule)):
    flags.Atom = 0001
    flags.Molecule = 0001
  use_primitive_setting_flags = [00000]
  if (space_group_info.group().conventional_centring_type_symbol() != "P"):
    use_primitive_setting_flags.append(0001)
  if (flags.Atom):
    for use_primitive_setting in use_primitive_setting_flags:
      test_atom(space_group_info, use_primitive_setting, verbose=flags.Verbose)
  if (flags.Molecule):
    for flag_f_part in (00000, 0001)[:]: #SWITCH
      for use_primitive_setting in use_primitive_setting_flags:
        test_molecule(space_group_info, use_primitive_setting, flag_f_part,
                      verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "HighSymmetry",
    "Atom",
    "Molecule"))

if (__name__ == "__main__"):
  run()
