from cctbx import xray
from cctbx.array_family import flex

def run():
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
  minimizer = xray.minimization.lbfgs(
    target_functor=xray.target_functors.intensity_correlation(y_obs),
    xray_structure=structure_shake,
    #use_special_position_constraints=False,
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


if __name__ == '__main__':
  run()
