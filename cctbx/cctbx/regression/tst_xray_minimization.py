from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys

def exercise(target_functor, space_group_info, anomalous_flag,
             gradient_flags,
             n_elements=3, d_min=3., shake_sigma=0.1,
             verbose=0):
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200,
    random_u_iso=True,
    random_occupancy=True)
  f_obs = abs(structure_ideal.structure_factors(
    anomalous_flag=anomalous_flag,
    d_min=d_min,
    algorithm="direct",
    cos_sin_table=True).f_calc())
  if (0 or verbose):
    print "structure_ideal:"
    structure_ideal.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_ideal.special_position_indices().size()
    print
  structure_shake = structure_ideal
  if (gradient_flags.site):
    structure_shake = structure_shake.random_modify_parmeters(
      "site", shake_sigma)
  if (gradient_flags.u_iso):
    structure_shake = structure_shake.random_modify_parmeters(
      "u_iso", shake_sigma)
  if (gradient_flags.occupancy):
    structure_shake = structure_shake.random_modify_parmeters(
      "occupancy", shake_sigma)
  assert tuple(structure_ideal.special_position_indices()) \
      == tuple(structure_shake.special_position_indices())
  if (0 or verbose):
    print "structure_shake:"
    structure_shake.show_summary().show_scatterers()
    print
  minimizer = xray.minimization.lbfgs(
    target_functor=target_functor(f_obs),
    gradient_flags=gradient_flags,
    xray_structure=structure_shake,
    structure_factor_algorithm="direct")
  if (0 or verbose):
    print "first:", minimizer.first_target_value
    print "final:", minimizer.final_target_value
    print
  assert minimizer.final_target_value < minimizer.first_target_value
  if (0 or verbose):
    print "minimized structure_shake:"
    structure_shake.show_summary().show_scatterers()
    print
  f_final = abs(f_obs.structure_factors_from_scatterers(
    xray_structure=structure_shake,
    algorithm="direct",
    cos_sin_table=True).f_calc())
  c = flex.linear_correlation(f_obs.data(), f_final.data())
  assert c.is_well_defined()
  if (0 or verbose):
    print "correlation:", c.coefficient()
    print
  assert c.coefficient() > 0.99

def run_call_back(flags, space_group_info):
  for target_functor in xray.target_functors.registry().values():
    for i_options in (1,2,4): #SWITCH
      gradient_flags = xray.structure_factors.gradient_flags(
        site=(i_options % 2),
        u_iso=(i_options/2 % 2),
        occupancy=(i_options/4 % 2))
      for anomalous_flag in (False, True)[:]: #SWITCH
        exercise(target_functor, space_group_info, anomalous_flag,
                 gradient_flags, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
