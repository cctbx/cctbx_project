from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys
import random

def finite_differences_site(cartesian_flag, target_ftor, structure,
                            delta=0.00001):
  unit_cell = structure.unit_cell()
  abc = unit_cell.parameters()[:3]
  derivatives = flex.vec3_double()
  for i_scatterer in structure.scatterers().indices():
    d_target_d_site = [0,0,0]
    for ix in xrange(3):
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
        target_result = target_ftor(f_calc, compute_derivatives=00000)
        target_values.append(target_result.target())
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      if (not cartesian_flag): derivative *= abc[ix]
      d_target_d_site[ix] = derivative
    derivatives.append(d_target_d_site)
  return derivatives

def finite_differences_u_star(target_ftor, structure,
                              delta=0.000001):
  derivatives = flex.sym_mat3_double()
  for i_scatterer in structure.scatterers().indices():
    d_target_d_u_star = [0,0,0,0,0,0]
    for iu in xrange(6):
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
        target_result = target_ftor(f_calc, compute_derivatives=00000)
        target_values.append(target_result.target())
      derivative = (target_values[1] - target_values[0]) / (2 * delta)
      d_target_d_u_star[iu] = derivative
    derivatives.append(d_target_d_u_star)
  return derivatives

def finite_differences_scalar(parameter_name, target_ftor, structure,
                              delta=0.00001):
  derivatives = flex.double()
  for i_scatterer in structure.scatterers().indices():
    target_values = []
    for d_sign in (-1, 1):
      modified_structure = structure.deep_copy_scatterers()
      ms = modified_structure.scatterers()[i_scatterer]
      if   (parameter_name == "u_iso"):
        ms.u_iso += d_sign * delta
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
      target_result = target_ftor(f_calc, compute_derivatives=00000)
      target_values.append(target_result.target())
    derivative = (target_values[1] - target_values[0]) / (2 * delta)
    derivatives.append(derivative)
  return derivatives

def flex_tuple_as_flex_double(flex_tuple):
  result = flex.double()
  for t in flex_tuple:
    result.extend(flex.double(tuple(t)))
  return result

def linear_regression_test(d_analytical, d_numerical, test_hard=0001,
                           slope_tolerance=1.e-3,
                           correlation_min=0.999,
                           verbose=0):
  if (type(d_analytical) != type(flex.double())):
    d_analytical = flex_tuple_as_flex_double(d_analytical)
  if (type(d_numerical) != type(flex.double())):
    d_numerical = flex_tuple_as_flex_double(d_numerical)
  if (0 or verbose):
    print "analytical:", tuple(d_analytical)
    print "numerical: ", tuple(d_numerical)
  if (    flex.max(flex.abs(d_analytical)) == 0
      and flex.max(flex.abs(d_numerical)) == 0):
    return
  regr = flex.linear_regression(d_analytical, d_numerical)
  corr = flex.linear_correlation(d_analytical, d_numerical).coefficient()
  assert regr.is_well_defined()
  if (abs(regr.slope() - 1) > slope_tolerance or corr < correlation_min):
    print "Error: finite difference mismatch:"
    print "slope:", regr.slope()
    print "correlation:", corr
    if (0 or verbose):
      for a, n in zip(d_analytical, d_numerical):
        print a, n
    assert not test_hard

def exercise(target_functor, parameter_name, space_group_info,
             anomalous_flag, cartesian_flag,
             n_elements=4, d_min=3., shake_sigma=0.25,
             test_hard=0001, verbose=0):
  if (parameter_name != "site" and cartesian_flag == 0001): return
  if (parameter_name == "fdp" and not anomalous_flag): return
  structure_ideal = random_structure.xray_structure(
    space_group_info,
    elements=(("N","C","O")*(n_elements/3+1))[:n_elements],
    volume_per_atom=100,
    random_f_prime_d_min=d_min,
    random_f_double_prime=anomalous_flag,
    anisotropic_flag=(parameter_name == "u_star"),
    random_u_iso=0001,
    random_occupancy=0001)
  f_obs = abs(structure_ideal.structure_factors(
    anomalous_flag=anomalous_flag, d_min=d_min, algorithm="direct").f_calc())
  if (0 or verbose):
    structure_ideal.show_summary().show_scatterers()
    print "n_special_positions:", \
          structure_ideal.special_position_indices().size()
  if (0 or verbose):
    f_obs.show_summary()
  structure_shake = structure_ideal.random_modify_parmeters(
    parameter_name, shake_sigma, vary_z_only=00000)
  assert tuple(structure_ideal.special_position_indices()) \
      == tuple(structure_shake.special_position_indices())
  if (0 or verbose):
    structure_shake.show_summary().show_scatterers()
  target_ftor = target_functor(f_obs)
  for structure in (structure_ideal, structure_shake)[:]: #SWITCH
    f_calc = f_obs.structure_factors_from_scatterers(
      xray_structure=structure,
      algorithm="direct").f_calc()
    target_result = target_ftor(f_calc, compute_derivatives=0001)
    if (structure == structure_ideal):
      assert abs(target_result.target()) < 1.e-5
  if (0 or verbose):
    try: print "scale factor = %.6g" % (target_result.scale_factor(),)
    except: pass
    try: print "correlation = %.6g" % (target_result.correlation(),)
    except: pass
    print "target = %.6g" % (target_result.target(),)
  sf = xray.structure_factors.gradients_direct(
    xray_structure=structure,
    miller_set=f_obs,
    d_target_d_f_calc=target_result.derivatives(),
    gradient_flags=xray.structure_factors.gradient_flags(
      site=(parameter_name=="site" or random.choice((0,1))),
      u_iso=(parameter_name=="u_iso" or random.choice((0,1))),
      u_aniso=(parameter_name=="u_star" or random.choice((0,1))),
      occupancy=(parameter_name=="occupancy" or random.choice((0,1))),
      fp=(parameter_name=="fp" or random.choice((0,1))),
      fdp=(parameter_name=="fdp"
           or (anomalous_flag and random.choice((0,1))))),
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
  elif (parameter_name == "u_star"):
    d_analytical = sf.d_target_d_u_star()
    d_numerical = finite_differences_u_star(target_ftor, structure)
  else:
    if (parameter_name == "occupancy"):
      d_analytical = sf.d_target_d_occupancy()
    elif (parameter_name == "u_iso"):
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

def run_call_back(flags, space_group_info):
  coordinate_systems = []
  if (flags.Frac): coordinate_systems.append(00000)
  if (flags.Cart): coordinate_systems.append(0001)
  if (len(coordinate_systems) == 0):
    coordinate_systems = [00000, 0001]
  for target_functor in xray.target_functors.registry().values():
    for parameter_name in ("site", "u_iso", "u_star", "occupancy",
                           "fp", "fdp")[:]: #SWITCH
      for anomalous_flag in (00000, 0001)[:]: #SWITCH
        for cartesian_flag in coordinate_systems:
          exercise(target_functor, parameter_name, space_group_info,
                   anomalous_flag, cartesian_flag,
                   verbose=flags.Verbose)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back, (
    "Frac",
    "Cart"))

if (__name__ == "__main__"):
  run()
