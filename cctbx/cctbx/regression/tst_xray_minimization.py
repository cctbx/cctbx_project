from cctbx.xray import minimization
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys, os

def exercise(target_functor, space_group_info, anomalous_flag,
             minimization_options,
             n_elements=3, d_min=3., shake_sigma=0.1,
             verbose=0):
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200)
  f_obs_array = abs(structure_ideal.structure_factors_direct(
    anomalous_flag=anomalous_flag, d_min=d_min).f_calc_array())
  if (0 or verbose):
    structure_ideal.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_ideal.special_position_indices().size()
  structure_shake = structure_ideal
  if (minimization_options.site):
    structure_shake = structure_shake.random_modify_parmeters(
      "site", shake_sigma)
  if (minimization_options.u_iso):
    structure_shake = structure_shake.random_modify_parmeters(
      "u_iso", shake_sigma)
  if (minimization_options.occupancy):
    structure_shake = structure_shake.random_modify_parmeters(
      "occupancy", shake_sigma)
  if (0 or verbose):
    structure_shake.show_summary().show_scatterers()
  assert tuple(structure_ideal.special_position_indices()) \
      == tuple(structure_shake.special_position_indices())
  if (0 or verbose):
    structure_shake.show_summary().show_scatterers()
  minimizer = minimization.lbfgs(
    target_functor(f_obs_array),
    minimization_options,
    structure_shake)
  if (0 or verbose):
    print "first:", minimizer.first_target_value
    print "final:", minimizer.final_target_value
  assert minimizer.final_target_value < minimizer.first_target_value
  f_final_array = abs(xray.structure_factors_direct(
    structure_shake, f_obs_array).f_calc_array())
  r = flex.linear_regression(f_obs_array.data(), f_final_array.data())
  assert r.is_well_defined()
  if (0 or verbose):
    print "correlation:", r.correlation()
  assert r.correlation() > 0.999

def run_call_back(flags, space_group_info):
  for target_functor in xray.target_functors.values():
    for i_options in (1,2,4): #SWITCH
      minimization_options = minimization.options(
        site=(i_options % 2),
        u_iso=(i_options/2 % 2),
        occupancy=(i_options/4 % 2))
      for anomalous_flag in (False, True)[:]: #SWITCH
        exercise(target_functor, space_group_info, anomalous_flag,
                 minimization_options, verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
