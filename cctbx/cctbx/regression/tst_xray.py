from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx import xray
from cctbx import crystal
from cctbx import adptbx
import cctbx.eltbx.xray_scattering
from cctbx import eltbx
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
  p1 = xs.asymmetric_unit_in_p1()
  assert p1.scatterers().size() == 2
  for i in xrange(2):
    assert p1.scatterers()[i].weight() == xs.scatterers()[i].weight()
  assert str(p1.space_group_info()) == "P 1"
  p1 = xs.expand_to_p1()
  assert p1.scatterers().size() == 9
  for i in xrange(2):
    assert p1.scatterers()[i].weight() != xs.scatterers()[i].weight()
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
  assert xs.n_parameters(xray.structure_factors.gradient_flags(default=0001)) \
         == 14
  g = flex.vec3_double(((0.1,0.2,0.3),(0.2,0.3,0.4)))
  xs.apply_special_position_ops_d_target_d_site(g)
  assert approx_equal(g[0], (0,0,0))
  assert approx_equal(g[1], (-0.05,0.05,0))
  xs.replace_scatterers(xs.scatterers()[:1])
  assert xs.scatterers().size() == 1
  assert tuple(xs.special_position_indices()) == (0,)
  sd = ys.scattering_dict()
  assert sd.lookup("Si").gaussian.n_ab() == 5
  sd = ys.scattering_dict(d_min=3)
  assert sd.lookup("Si").gaussian.n_ab() == 4
  sd = ys.scattering_dict(custom_dict={"Si":eltbx.xray_scattering.gaussian(1)})
  assert sd.lookup("Si").gaussian.n_ab() == 0

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
                                    element_type,
                                    n_elements=3,
                                    volume_per_atom=1000,
                                    d_min=3,
                                    anomalous_flag=0,
                                    anisotropic_flag=0,
                                    verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=[element_type]*n_elements,
    volume_per_atom=volume_per_atom,
    min_distance=5,
    general_positions_only=1,
    random_f_prime_d_min=d_min-1,
    random_f_prime_scale=0.6,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=anisotropic_flag,
    random_u_iso=0001,
    random_u_iso_scale=.3,
    random_u_cart_scale=.3,
    random_occupancy=0001)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_obs_exact = structure.structure_factors(
    d_min=d_min, algorithm="direct",
    cos_sin_table=00000).f_calc()
  assert f_obs_exact.anomalous_flag() == anomalous_flag
  f_obs_simple = xray.ext.structure_factors_simple(
    f_obs_exact.unit_cell(),
    f_obs_exact.space_group(),
    f_obs_exact.indices(),
    structure.scatterers(),
    structure.scattering_dict()).f_calc()
  if (0 or verbose):
    for i,h in f_obs_exact.indices().items():
      print h
      print f_obs_simple[i]
      print f_obs_exact.data()[i]
      if (abs(f_obs_simple[i]-f_obs_exact.data()[i]) >= 1.e-10):
        print "MISMATCH"
      print
  mismatch = flex.max(flex.abs(f_obs_exact.data() - f_obs_simple))
  assert mismatch < 1.e-10, mismatch
  f_obs_table = f_obs_exact.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct",
    cos_sin_table=0001).f_calc()
  ls = xray.targets_least_squares_residual(
    abs(f_obs_exact).data(), f_obs_table.data(), 00000, 1)
  if (0 or verbose):
    print "r-factor:", ls.target()
  assert ls.target() < 1.e-4

def run_call_back(flags, space_group_info):
  for element_type in ("Se", "const"):
    for anomalous_flag in [0,1]:
      for anisotropic_flag in [0,1]:
        for with_shift in [0,1]:
          if (with_shift):
            sgi = debug_utils.random_origin_shift(space_group_info)
          else:
            sgi = space_group_info
          exercise_from_scatterers_direct(
            space_group_info=sgi,
            element_type=element_type,
            anomalous_flag=anomalous_flag,
            anisotropic_flag=anisotropic_flag,
            verbose=flags.Verbose)

def run():
  exercise_structure()
  exercise_u_extra()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print "OK"

if (__name__ == "__main__"):
  run()
