from cctbx import translation_search
from cctbx import uctbx
from cctbx import sgtbx
from cctbx.array_family import flex

def exercise_symmetry_flags():
  for i_flags in xrange(4):
    is_isotropic_search_model = (i_flags % 2 != 0)
    have_f_part = ((i_flags/2) % 2 != 0)
    f = translation_search.symmetry_flags(
      is_isotropic_search_model, have_f_part)
    assert f.is_isotropic_search_model() == is_isotropic_search_model
    assert f.have_f_part() == have_f_part
    assert f.use_space_group_symmetry() == is_isotropic_search_model
    assert f.use_normalizer_k2l() \
        == (is_isotropic_search_model and (not have_f_part))
    assert f.use_structure_seminvariants() == (not have_f_part)

def exercise_map_gridding():
  space_group_type = sgtbx.space_group_info("P 21 21 21").type()
  miller_indices_f_obs = flex.miller_index(((3,4,5),(4,5,6)))
  g = translation_search.map_gridding(
    unit_cell=uctbx.unit_cell((10,13,17)),
    space_group_type=space_group_type,
    symmetry_flags=translation_search.symmetry_flags(True, False),
    resolution_factor=1./3,
    miller_indices_f_obs=miller_indices_f_obs,
    max_prime=5)
  assert g.target() == (20,30,36)
  assert g.quarter() == (40,60,72)
  assert g.eighth() == (60,90,108)
  return space_group_type.group(), miller_indices_f_obs, g

def exercise_fast_nv1995():
  space_group, miller_indices_f_obs, gridding = exercise_map_gridding()
  f = translation_search.fast_nv1995(
    gridding=gridding,
    space_group=space_group,
    anomalous_flag=False,
    miller_indices_f_obs=miller_indices_f_obs,
    f_obs=flex.double((1,2)),
    f_part=flex.complex_double(),
    miller_indices_p1_f_calc=flex.miller_index(((1,2,3),)),
    p1_f_calc=flex.complex_double((12,)))
  assert f.target_map().all() == gridding.target()

def run():
  exercise_symmetry_flags()
  exercise_map_gridding()
  exercise_fast_nv1995()
  print "OK"

if (__name__ == "__main__"):
  run()
