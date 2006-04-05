from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys

def exercise(target_functor, space_group_info, anomalous_flag,
             gradient_flags, u_penalty, occupancy_penalty,
             n_elements=3, d_min=3., shake_sigma=0.1,
             verbose=0):
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=("Se",)*n_elements,
    volume_per_atom=200,
    random_u_iso=True,
    random_u_iso_min=0.05,
    random_occupancy=True)
  for scatterer in structure_ideal.scatterers():
      scatterer.flags.set_grad_site(gradient_flags.site)
      scatterer.flags.set_grad_u_iso(gradient_flags.u_iso)
      scatterer.flags.set_grad_u_aniso(gradient_flags.u_aniso)
      scatterer.flags.set_grad_occupancy(gradient_flags.occupancy)
      scatterer.flags.set_grad_fp(gradient_flags.fp)
      scatterer.flags.set_grad_fdp(gradient_flags.fdp)

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
    structure_shake = structure_shake.random_modify_parameters(
      "site", shake_sigma)
  if (gradient_flags.u_iso):
    structure_shake = structure_shake.random_modify_parameters(
      "u_iso", shake_sigma)
  if (gradient_flags.occupancy):
    structure_shake = structure_shake.random_modify_parameters(
      "occupancy", shake_sigma)
    if (occupancy_penalty is not None):
      structure_shake.scatterers()[-1].occupancy = 0
  assert tuple(structure_ideal.special_position_indices()) \
      == tuple(structure_shake.special_position_indices())
  if (0 or verbose):
    print "structure_shake:"
    structure_shake.show_summary().show_scatterers()
    print
  minimizer = xray.minimization.lbfgs(
    target_functor=target_functor(f_obs),
    xray_structure=structure_shake,
    u_penalty=u_penalty,
    occupancy_penalty=occupancy_penalty,
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
    label = gradient_flags.string_of_true()
    if (u_penalty is not None):
      label += ","+u_penalty.__class__.__name__
    if (anomalous_flag):
      label += ",anomalous"
    print "correlation: %10.8f" % c.coefficient(), label
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
        u_penalty_types = [None]
        sqrt_u_isos = [False]
        if (gradient_flags.u_iso):
          u_penalty_types.extend([
            xray.minimization.u_penalty_exp(),
            xray.minimization.u_penalty_singular_at_zero()])
          sqrt_u_isos.append(True)
        occupancy_penalty_types = [None]
        if (gradient_flags.occupancy):
          occupancy_penalty_types.append(
            xray.minimization.occupancy_penalty_exp())
        for sqrt_u_iso in sqrt_u_isos:
          for u_penalty in u_penalty_types:
            gradient_flags.sqrt_u_iso = sqrt_u_iso
            for occupancy_penalty in occupancy_penalty_types:
              exercise(target_functor, space_group_info, anomalous_flag,
                       gradient_flags,
                       u_penalty=u_penalty,
                       occupancy_penalty=occupancy_penalty,
                       verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)

if (__name__ == "__main__"):
  run()
