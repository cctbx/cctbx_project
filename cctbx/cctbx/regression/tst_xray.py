from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
from cctbx import crystal
from cctbx import adptbx
from cctbx.array_family import flex
from scitbx.test_utils import approx_equal
import sys

def exercise_structure():
  cs = crystal.symmetry((5.01, 5.01, 5.47, 90, 90, 120), "P 62 2 2")
  sp = crystal.special_position_settings(cs)
  scatterers = flex.xray_scatterer((
    xray.scatterer("Si1", (1./2, 1./2, 1./3)),
    xray.scatterer("O1", (0.19700, -0.19700, 0.83333))))
  xs = xray.structure(sp, scatterers)
  assert xs.scatterers().size() == 2
  assert tuple(xs.special_position_indices()) == (0, 1)
  xs.all_apply_symmetry()
  assert tuple(xs.special_position_indices()) == (0, 1)
  ys = xs.deep_copy_scatterers()
  ys.add_scatterers(ys.scatterers())
  assert ys.scatterers().size() == 4
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3)
  ys.add_scatterer(ys.scatterers()[0])
  assert ys.scatterers().size() == 5
  assert tuple(ys.special_position_indices()) == (0, 1, 2, 3, 4)
  sx = xs.primitive_setting()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 62 2 2"
  sx = xs.change_hand()
  assert sx.unit_cell().is_similar_to(xs.unit_cell())
  assert str(sx.space_group_info()) == "P 64 2 2"
  assert approx_equal(sx.scatterers()[0].site, (-1./2, -1./2, -1./3))
  assert approx_equal(sx.scatterers()[1].site, (-0.19700, 0.19700, -0.833333))
  p1 = xs.expand_to_p1()
  assert p1.scatterers().size() == 9
  sh = p1.apply_shift((0.2,0.3,-1/6.))
  assert approx_equal(sh.scatterers()[0].site, (0.7,0.8,1/6.))
  assert approx_equal(sh.scatterers()[3].site, (0.3970,0.1030,2/3.))
  sl = sh[:1]
  assert sl.scatterers().size() == 1
  assert sl.scatterers()[0].label == sh.scatterers()[0].label
  sl = sh[1:4]
  assert sl.scatterers().size() == 3
  for i in xrange(3):
    assert sl.scatterers()[i].label == sh.scatterers()[i+1].label
  xs.scatterers().set_occupancies(flex.double((0.5,0.2)))
  s = xs.sort(by_value="occupancy")
  assert approx_equal(s.scatterers().extract_occupancies(), (0.2,0.5))
  assert s.scatterers()[0].label == "O1"
  assert s.scatterers()[1].label == "Si1"
  aw = xs.atomic_weights()
  assert approx_equal(aw, (28.086, 15.999))
  center_of_mass = xs.center_of_mass(atomic_weights=aw)
  assert approx_equal(center_of_mass.elems, (1.335228, 1.071897, 2.815899))
  center_of_mass = xs.center_of_mass()
  assert approx_equal(center_of_mass.elems, (1.335228, 1.071897, 2.815899))
  ys = xs.apply_shift(xs.unit_cell().fractionalize((-center_of_mass).elems))
  assert approx_equal(ys.center_of_mass().elems, (0,0,0))
  ys = xray.structure(xs)
  assert ys.atomic_weights().size() == 0
  assert ys.center_of_mass().elems == (0,0,0)
  ys = xray.structure(sp, scatterers)
  ys.set_occupancy(1, 0.5)
  assert approx_equal(ys.scatterers()[1].weight(),0.25)
  ys.shift_occupancy(1, -0.1)
  assert approx_equal(ys.scatterers()[1].weight(),0.2)

def exercise_u_extra():
  d_min = 9
  grid_resolution_factor = 1/3.
  for quality_factor in (1,2,4,8,10,100,200,1000):
    u_extra = xray.calc_u_extra(d_min, grid_resolution_factor, quality_factor)
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        d_min=d_min,
        grid_resolution_factor=grid_resolution_factor,
        u_extra=u_extra))
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        d_min=d_min,
        grid_resolution_factor=grid_resolution_factor,
        b_extra=adptbx.u_as_b(u_extra)))
    assert approx_equal(
      quality_factor,
      xray.structure_factors.quality_factor_from_any(
        quality_factor=quality_factor))

def exercise_from_scatterers_direct(space_group_info,
                                    n_elements=3,
                                    volume_per_atom=1000,
                                    d_min=2,
                                    fdp_flag=0,
                                    anisotropic_flag=0,
                                    verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=1,
    random_f_prime_d_min=d_min-1,
    random_f_prime_scale=0.6,
    random_f_double_prime=fdp_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=0001,
    random_u_iso_scale=.3,
    random_u_cart_scale=.3,
    random_occupancy=0001)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_obs_exact = abs(structure.structure_factors(
    d_min=d_min, anomalous_flag=fdp_flag, direct=0001,
    cos_sin_table=00000).f_calc())
  f_obs_table = f_obs_exact.structure_factors_from_scatterers(
    xray_structure=structure,
    direct=0001,
    cos_sin_table=0001).f_calc()
  ls = xray.targets_least_squares_residual(
    f_obs_exact.data(), f_obs_table.data(), 00000, 1)
  if (0 or verbose):
    print "r-factor:", ls.target()
  assert ls.target() < 1.e-4

def run_call_back(flags, space_group_info):
  for fdp_flag in [0,1]:
    for anisotropic_flag in [0,1]:
      exercise_from_scatterers_direct(
        space_group_info=space_group_info,
        fdp_flag=fdp_flag,
        anisotropic_flag=anisotropic_flag,
        verbose=flags.Verbose)

def run():
  exercise_structure()
  exercise_u_extra()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
