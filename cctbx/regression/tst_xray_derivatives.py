from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import random
import math
import sys
from cctbx import adptbx
from six.moves import range
from six.moves import zip

def finite_differences_site(cartesian_flag, target_ftor, structure,
                            delta=0.00001):
  unit_cell = structure.unit_cell()
  abc = unit_cell.parameters()[:3]
  derivatives = flex.vec3_double()
  for i_scatterer in range(structure.scatterers().size()):
    d_target_d_site = [0,0,0]
    for ix in range(3):
      target_values = []
      for d_sign in (-1, 1):
        modified_structure = structure.deep_copy_scatterers()
        ms = modified_structure.scatterers()[i_scatterer]
        site = list(ms.site)
        if (not cartesian_flag):
          site[ix] += d_sign * delta / abc[ix]
        else:
          site_cart = list(unit_cell.orthogonalize(site))
          site_cart[ix] += d_sign * delta
          site = unit_cell.fractionalize(site_cart)
        ms.site = site
        f_calc = target_ftor.f_obs().structure_factors_from_scatterers(
            xray_structure=modified_structure,
            algorithm="direct").f_calc()
        target_result = target_ftor(f_calc, compute_derivatives=False)
        target_values.append(target_result.target())
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      if (not cartesian_flag): derivative *= abc[ix]
      d_target_d_site[ix] = derivative
    derivatives.append(d_target_d_site)
  return derivatives

def finite_differences_u_star(target_ftor, structure,
                              delta=0.000001):
  derivatives = flex.sym_mat3_double()
  for i_scatterer in range(structure.scatterers().size()):
    d_target_d_u_star = [0,0,0,0,0,0]
    if(structure.scatterers()[i_scatterer].flags.use_u_aniso()):
       for iu in range(6):
         target_values = []
         for d_sign in (-1, 1):
           modified_structure = structure.deep_copy_scatterers()
           ms = modified_structure.scatterers()[i_scatterer]
           u_star = list(ms.u_star)
           u_star[iu] += d_sign * delta
           ms.u_star = u_star
           f_calc = target_ftor.f_obs().structure_factors_from_scatterers(
             xray_structure=modified_structure,
             algorithm="direct").f_calc()
           target_result = target_ftor(f_calc, compute_derivatives=False)
           target_values.append(target_result.target())
         derivative = (target_values[1] - target_values[0]) / (2 * delta)
         d_target_d_u_star[iu] = derivative
       derivatives.append(d_target_d_u_star)
    else:
       derivatives.append(d_target_d_u_star)
  return derivatives

def finite_differences_scalar(parameter_name, target_ftor, structure,
                              delta=0.00001):
  derivatives = flex.double()
  for i_scatterer in range(structure.scatterers().size()):
    target_values = []
    for d_sign in (-1, 1):
      modified_structure = structure.deep_copy_scatterers()
      ms = modified_structure.scatterers()[i_scatterer]
      if   (parameter_name == "u_iso" and ms.flags.use_u_iso()):
        ms.u_iso += d_sign * delta
      elif (parameter_name == "tan_u_iso" and ms.flags.use_u_iso()):
        assert ms.u_iso >= 0
        x = math.tan(ms.u_iso*math.pi/adptbx.b_as_u(ms.flags.param)-math.pi/2) + d_sign * delta
        ms.u_iso = adptbx.b_as_u(ms.flags.param)/math.pi*(math.atan(x)+math.pi/2)
      elif (parameter_name == "occupancy"):
        ms.occupancy += d_sign * delta
      elif (parameter_name == "fp"):
        ms.fp = ms.fp + d_sign * delta
      elif (parameter_name == "fdp"):
        ms.fdp = ms.fdp + d_sign * delta
      else:
        raise RuntimeError
      f_calc = target_ftor.f_obs().structure_factors_from_scatterers(
        xray_structure=modified_structure,
        algorithm="direct").f_calc()
      target_result = target_ftor(f_calc, compute_derivatives=False)
      target_values.append(target_result.target())
    derivative = (target_values[1] - target_values[0]) / (2 * delta)
    derivatives.append(derivative)
  return derivatives

def flex_tuple_as_flex_double(flex_tuple):
  result = flex.double()
  for t in flex_tuple:
    result.extend(flex.double(tuple(t)))
  return result

def linear_regression_test(d_analytical, d_numerical, test_hard=True,
                           slope_tolerance=1.e-3,
                           correlation_min=0.999,
                           verbose=0):
  if (type(d_analytical) != type(flex.double())):
    d_analytical = flex_tuple_as_flex_double(d_analytical)
  if (type(d_numerical) != type(flex.double())):
    d_numerical = flex_tuple_as_flex_double(d_numerical)
  if (0 or verbose):
    print("analytical:", tuple(d_analytical))
    print("numerical: ", tuple(d_numerical))
  if (    flex.max(flex.abs(d_analytical)) == 0
      and flex.max(flex.abs(d_numerical)) == 0):
    return
  regr = flex.linear_regression(d_analytical, d_numerical)
  corr = flex.linear_correlation(d_analytical, d_numerical).coefficient()
  assert regr.is_well_defined()
  if (abs(regr.slope() - 1) > slope_tolerance or corr < correlation_min):
    print("Error: finite difference mismatch:")
    print("slope:", regr.slope())
    print("correlation:", corr)
    if (0 or verbose):
      for a, n in zip(d_analytical, d_numerical):
        print(a, n)
    assert not test_hard

def exercise(target_functor, data_type, parameter_name, space_group_info,
             anomalous_flag,
             cartesian_flag,
             n_elements=9,
             d_min=2.5,
             shake_sigma=0.25,
             test_hard=True, verbose=0):
  assert data_type == 'F' or data_type == 'F^2'
  if (data_type == 'F^2'
       and not target_functor == xray.unified_least_squares_residual): return
  if (parameter_name != "site" and cartesian_flag == True): return
  if (parameter_name == "fdp" and not anomalous_flag): return
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=(("O","N","C")*(n_elements))[:n_elements],
    volume_per_atom=100,
    random_f_prime_d_min=d_min,
    random_f_double_prime=anomalous_flag,
    use_u_aniso = True,
    use_u_iso = False,
    random_u_cart_scale = 0.3,
    random_u_iso = False,
    random_occupancy=True)
  if(parameter_name in ["u_star", "u_iso"]): shake_sigma = shake_sigma / 2.
  random_structure.random_modify_adp_and_adp_flags(
                             scatterers         = structure_ideal.scatterers(),
                             random_u_iso_scale = 0.3,
                             random_u_iso_min   = 0.0,
                             parameter_name     = parameter_name)
  rnd_f_calc = structure_ideal.structure_factors(
    anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc()
  if data_type == 'F':
    y_obs = abs(rnd_f_calc)
  elif data_type == 'F^2':
    y_obs = rnd_f_calc.norm()
    y_obs.set_observation_type_xray_intensity()
  structure_shake = structure_ideal.random_modify_parameters(
        parameter_name, shake_sigma, vary_z_only=False)
  assert tuple(structure_ideal.special_position_indices()) \
      == tuple(structure_shake.special_position_indices())
  target_ftor = target_functor(y_obs)
  for structure in (structure_ideal, structure_shake)[:]: #SWITCH
    f_calc = y_obs.structure_factors_from_scatterers(
      xray_structure=structure,
      algorithm="direct").f_calc()
    target_result = target_ftor(f_calc, compute_derivatives=True)
    if (structure == structure_ideal):
      assert abs(target_result.target()) < 1.e-5

  gradient_flags=xray.structure_factors.gradient_flags(
     site=(parameter_name=="site" or random.choice((False,True))),
     u_iso=(parameter_name=="u_iso" or random.choice((False,True))),
     u_aniso=(parameter_name=="u_star" or random.choice((False,True))),
     occupancy=(parameter_name=="occupancy" or random.choice((False,True))),
     fp=(parameter_name=="fp" or random.choice((False,True))),
     fdp=(parameter_name=="fdp" or (anomalous_flag
                                    and random.choice((False,True)))))
  xray.set_scatterer_grad_flags(scatterers = structure.scatterers(),
                                site       = gradient_flags.site,
                                u_iso      = gradient_flags.u_iso,
                                u_aniso    = gradient_flags.u_aniso,
                                occupancy  = gradient_flags.occupancy,
                                fp         = gradient_flags.fp,
                                fdp        = gradient_flags.fdp)
  grad_flags_counts = xray.scatterer_grad_flags_counts(structure.scatterers())
  sf = xray.structure_factors.gradients_direct(
    xray_structure=structure,
    u_iso_refinable_params=None,
    miller_set=y_obs,
    d_target_d_f_calc=target_result.derivatives(),
    n_parameters=0)
  if (parameter_name == "site"):
    d_analytical = sf.d_target_d_site_frac()
    if (cartesian_flag): # average d_analytical or d_numerical, but not both
      structure_ideal.apply_special_position_ops_d_target_d_site(d_analytical)
    if (cartesian_flag):
      d_analytical = sf.d_target_d_site_cart()
    d_numerical = finite_differences_site(
      cartesian_flag, target_ftor, structure)
    if (not cartesian_flag): # aver. d_analytical or d_numerical, but not both
      structure_ideal.apply_special_position_ops_d_target_d_site(d_numerical)
  elif (parameter_name == "u_star" and grad_flags_counts.use_u_aniso > 0):
    d_analytical = sf.d_target_d_u_star()
    d_numerical = finite_differences_u_star(target_ftor, structure)
  else:
   if (parameter_name == "occupancy"):
     d_analytical = sf.d_target_d_occupancy()
   elif (parameter_name == "u_iso" and grad_flags_counts.use_u_iso > 0):
     d_analytical = sf.d_target_d_u_iso()
   elif (parameter_name == "fp"):
     d_analytical = sf.d_target_d_fp()
   elif (parameter_name == "fdp"):
     d_analytical = sf.d_target_d_fdp()
   else:
     raise RuntimeError
   d_numerical = finite_differences_scalar(
     parameter_name, target_ftor, structure)
  linear_regression_test(d_analytical, d_numerical, test_hard,
                           verbose=verbose)
  if (parameter_name == "u_iso" and grad_flags_counts.use_u_iso > 0):
    u_iso_refinable_params = flex.double()
    for scatterer in structure.scatterers():
        scatterer.flags.set_tan_u_iso(True)
        scatterer.flags.set_grad_u_iso(gradient_flags.u_iso)
        param = random.randint(90,120)
        scatterer.flags.param= param
        value = math.tan(scatterer.u_iso*math.pi/adptbx.b_as_u(param)-math.pi/2)
        u_iso_refinable_params.append(value)
    sf = xray.structure_factors.gradients_direct(
      xray_structure=structure,
      u_iso_refinable_params = u_iso_refinable_params,
      miller_set=y_obs,
      d_target_d_f_calc=target_result.derivatives(),
      n_parameters=0)
    d_analytical = sf.d_target_d_u_iso()
    d_numerical = finite_differences_scalar("tan_u_iso",target_ftor, structure)
    linear_regression_test(d_analytical, d_numerical, test_hard,
                           verbose=verbose)

def run_call_back(flags, space_group_info):
  coordinate_systems = []
  if (flags.Frac): coordinate_systems.append(False)
  if (flags.Cart): coordinate_systems.append(True)
  if (len(coordinate_systems) == 0):
    coordinate_systems = [False, True]
  for parameter_name in ("site", "u_iso", "u_star", "occupancy",
                         "fp", "fdp")[:]: #SWITCH
    for anomalous_flag in (False, True)[:]: #SWITCH
       for cartesian_flag in coordinate_systems:
         for target_functor in xray.target_functors.registry().values():
           if(parameter_name == "u_iso"):  use_u_iso = True
           if(parameter_name == "u_star"): use_u_aniso = True
           for data_type in ('F', 'F^2'):
            exercise(target_functor,
                     data_type,
                     parameter_name,
                     space_group_info,
                     anomalous_flag,
                     cartesian_flag,
                     verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "Frac",
    "Cart"))

if (__name__ == "__main__"):
  run()
