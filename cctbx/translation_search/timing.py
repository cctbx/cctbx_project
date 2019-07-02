from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx import miller
from cctbx import translation_search
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.utils import show_times
import omptbx

def run():
  target_structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("I432"),
    elements=['C']*6+['O'],
    use_u_iso=False,
    use_u_aniso=True,
  )
  shift = tuple(flex.random_double(3))
  print("shift to be found: (%.3f, %.3f, %.3f)" % shift)
  target_structure_in_p1 = target_structure.expand_to_p1().apply_shift(shift)
  miller_indices = miller.build_set(
    crystal_symmetry=target_structure,
    anomalous_flag=True,
    d_min=0.8)
  f_obs = miller_indices.structure_factors_from_scatterers(
    xray_structure=target_structure,
    algorithm="direct").f_calc().amplitudes()
  miller_indices_in_p1 = miller.build_set(
    crystal_symmetry=target_structure_in_p1,
    anomalous_flag=True,
    d_min=0.8)
  f_calc = miller_indices_in_p1.structure_factors_from_scatterers(
      xray_structure=target_structure_in_p1,
      algorithm="direct").f_calc()
  crystal_gridding = f_calc.crystal_gridding(
    symmetry_flags=translation_search.symmetry_flags(
      is_isotropic_search_model=False,
      have_f_part=False),
    resolution_factor=1/2
  )
  omptbx.env.num_threads = 1
  t_fast_tf = show_times()
  fast_tf_map = translation_search.fast_nv1995(
    gridding=crystal_gridding.n_real(),
    space_group=f_obs.space_group(),
    anomalous_flag=f_obs.anomalous_flag(),
    miller_indices_f_obs=f_obs.indices(),
    f_obs=f_obs.data(),
    f_part=flex.complex_double(), ## no sub-structure is already fixed
    miller_indices_p1_f_calc=f_calc.indices(),
    p1_f_calc=f_calc.data()).target_map()
  print()
  print("Fast translation function")
  t_fast_tf()
  t_cross_corr = show_times()
  for op in target_structure.space_group():
    f, op_times_f = f_calc.original_and_transformed(op)
    cross_corr_map = miller.fft_map(crystal_gridding,
                                    f * op_times_f.conjugate().data())
  print()
  print("Traditional cross-correlation")
  t_cross_corr()





if __name__ == '__main__':
  run()
