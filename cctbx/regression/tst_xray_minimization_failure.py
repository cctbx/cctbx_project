from cctbx import xray
from cctbx.array_family import flex
from cctbx.regression import tst_xray_minimization
from cctbx import sgtbx


def test1():
  import libtbx.easy_pickle
  import os.path
  structure_ideal = libtbx.easy_pickle.load(
    os.path.expanduser('structure_ideal.pickle'))
  structure_start = libtbx.easy_pickle.load(
    os.path.expanduser('structure_start.pickle'))
  if 1:
    structure_ideal.show_scatterers()
    structure_start.show_scatterers()
  rnd_f_calc = structure_ideal.structure_factors(
    anomalous_flag=False,
    d_min=2.5,
    algorithm="direct",
    cos_sin_table=True).f_calc()
  y_obs = abs(rnd_f_calc)
  structure_shake = structure_start.deep_copy_scatterers()
  structure_shake_bis = structure_ideal.deep_copy_scatterers()
  tst_xray_minimization.shift_u_aniso(structure_shake_bis, 0.001)
  minimizer = xray.minimization.lbfgs(
    target_functor=xray.target_functors.intensity_correlation(y_obs),
    xray_structure=structure_shake,
    use_special_position_constraints=False,
    occupancy_penalty=None,
    structure_factor_algorithm="direct")
  assert minimizer.final_target_value < minimizer.first_target_value
  if 1:
    structure_shake.show_scatterers()
  f_final = y_obs.structure_factors_from_scatterers(
    xray_structure=structure_shake,
    algorithm="direct",
    cos_sin_table=True).f_calc()
  f_final = abs(f_final)
  c = flex.linear_correlation(y_obs.data(), f_final.data())
  assert c.is_well_defined()
  c_coefficient = c.coefficient()
  assert c_coefficient > 0.999

def test2():
  import random
  import libtbx.easy_pickle
  random.setstate(libtbx.easy_pickle.load('random_generator_state.pickle'))
  tst_xray_minimization.exercise(
    target_functor=xray.target_functors.intensity_correlation,
    data_type='F',
    space_group_info= sgtbx.space_group_info('R -3 c :H'),
    anomalous_flag=False,
    gradient_flags=xray.structure_factors.gradient_flags(
                                        site      = 0,
                                        u_iso     = 0,
                                        u_aniso   = 1,
                                        occupancy = 0),
    occupancy_penalty=None,
    d_min=2.5 )

def run():
  test2()
  test1()
  print 'OK'


if __name__ == '__main__':
  run()
