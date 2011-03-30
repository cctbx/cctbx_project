from cctbx import translation_search
from cctbx import sgtbx
from cctbx.array_family import flex

def exercise_symmetry_flags():
  for i_flags in xrange(4):
    is_isotropic_search_model = (i_flags % 2 != 0)
    have_f_part = ((i_flags//2) % 2 != 0)
    f = translation_search.symmetry_flags(
      is_isotropic_search_model=is_isotropic_search_model,
      have_f_part=have_f_part)
    assert f.is_isotropic_search_model() == is_isotropic_search_model
    assert f.have_f_part() == have_f_part
    assert f.use_space_group_symmetry() == is_isotropic_search_model
    assert f.use_normalizer_k2l() \
        == (is_isotropic_search_model and (not have_f_part))
    assert f.use_seminvariants() == (not have_f_part)

def exercise_fast_nv1995():
  space_group = sgtbx.space_group_info("P 21 21 21").group()
  miller_indices_f_obs = flex.miller_index(((3,4,5),(4,5,6)))
  f = translation_search.fast_nv1995(
    gridding=(20,20,36),
    space_group=space_group,
    anomalous_flag=False,
    miller_indices_f_obs=miller_indices_f_obs,
    f_obs=flex.double((1,2)),
    f_part=flex.complex_double(),
    miller_indices_p1_f_calc=flex.miller_index(((1,2,3),)),
    p1_f_calc=flex.complex_double((12,)))
  assert f.target_map().all() == (20,20,36)

def exercise_fast_terms():
  space_group = sgtbx.space_group_info("P 21 21 21").group()
  miller_indices_f_obs = flex.miller_index(((3,4,5),(4,5,6)))
  f = translation_search.fast_terms(
    gridding=(20,20,36),
    anomalous_flag=False,
    miller_indices_p1_f_calc=flex.miller_index(((1,2,3),)),
    p1_f_calc=flex.complex_double((12,)))
  assert f.summation(
    space_group=space_group,
    miller_indices_f_obs=miller_indices_f_obs,
    m=flex.double((1,2)),
    f_part=None,
    squared_flag=False).fft().accu_real_copy().all() == (20,20,36)

def run():
  exercise_symmetry_flags()
  exercise_fast_nv1995()
  exercise_fast_terms()
  print "OK"

if (__name__ == "__main__"):
  run()
