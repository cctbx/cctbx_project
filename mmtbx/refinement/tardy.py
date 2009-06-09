from __future__ import division
from mmtbx.dynamics import \
  kinetic_energy_as_temperature, \
  temperature_as_kinetic_energy
from mmtbx.dynamics.constants import akma_time_as_pico_seconds
from cctbx import xray
import cctbx.geometry_restraints
from cctbx.array_family import flex
import scitbx.rigid_body.essence.tardy
from scitbx.graph import tardy_tree
from scitbx import matrix
from libtbx.utils import Sorry
from libtbx.str_utils import format_value
from libtbx import group_args

master_phil_str = """\
  xray_weight_factor = 10
    .type = float
  start_temperature_kelvin = 2500
    .type = float
  final_temperature_kelvin = 300
    .type = float
  velocity_scaling = True
    .type = bool
  temperature_cap_factor = 1.5
    .type = float
  excessive_temperature_factor = 5
    .type = float
  number_of_cooling_steps = 500
    .type = int
  number_of_time_steps = 1
    .type = int
  time_step_pico_seconds = 0.001
    .type = float
  temperature_degrees_of_freedom = *cartesian constrained
    .type = choice
    .optional = False
  minimization_max_iterations = 0
    .type = int
  nonbonded_attenuation_factor = 0.75
    .type = float
  omit_bonds_with_slack_greater_than = 0
    .type = float
  constrain_dihedrals_with_sigma_less_than = 10
    .type = float
  near_singular_hinges_angular_tolerance_deg = 5
    .type = float
"""

class potential_object(object):

  def __init__(O,
        xray_weight_factor,
        nonbonded_attenuation_factor,
        fmodels,
        model,
        target_weights,
        reduced_geo_manager):
    O.xray_weight_factor = xray_weight_factor
    O.fmodels = fmodels
    O.model = model
    O.weights = target_weights.xyz_weights_result
    O.reduced_geo_manager = reduced_geo_manager
    if (nonbonded_attenuation_factor is None):
      O.custom_nonbonded_function = None
    else:
      assert nonbonded_attenuation_factor > 0
      assert nonbonded_attenuation_factor <= 1
      nonbonded_function = model.restraints_manager.geometry.nonbonded_function
      assert isinstance(
        nonbonded_function,
        cctbx.geometry_restraints.prolsq_repulsion_function)
      O.custom_nonbonded_function = nonbonded_function.customized_copy(
        k_rep=nonbonded_attenuation_factor)
    O.fmodels.create_target_functors()
    O.fmodels.prepare_target_functors_for_minimization()
    O.last_sites_moved = None
    O.f = None
    O.g = None
    O.last_grms = None

  def e_pot(O, sites_moved):
    if (O.last_sites_moved is not sites_moved):
      O.last_sites_moved = sites_moved
      xs = O.fmodels.fmodel_xray().xray_structure
      assert len(sites_moved) == xs.scatterers().size()
      sites_cart = flex.vec3_double(sites_moved)
      xs.set_sites_cart(sites_cart=sites_cart)
      O.fmodels.update_xray_structure(update_f_calc=True)
      xs.scatterers().flags_set_grads(state=False)
      xs.scatterers().flags_set_grad_site(
        iselection=flex.size_t_range(xs.scatterers().size()))
      tg = O.fmodels.target_and_gradients(
        weights=O.weights,
        compute_gradients=True)
      O.f = tg.target() * O.xray_weight_factor
      O.g = tg.gradients() * O.xray_weight_factor
      assert O.g.size() == len(sites_moved) * 3
      reduced_geo_energies = O.reduced_geo_manager.energies_sites(
        sites_cart=sites_cart,
        compute_gradients=True)
      other_energies = O.model.restraints_manager.energies_sites(
        sites_cart=sites_cart,
        geometry_flags=cctbx.geometry_restraints.flags.flags(nonbonded=True),
        custom_nonbonded_function=O.custom_nonbonded_function,
        compute_gradients=True)
      nfw = other_energies.normalization_factor * O.weights.w
      O.f += other_energies.target * O.weights.w
      O.f += reduced_geo_energies.target * nfw
      gg = other_energies.gradients * O.weights.w
      gg += reduced_geo_energies.gradients * nfw
      assert nfw != 0
      scale = 1 / nfw
      O.last_grms = group_args(
        geo=scale*flex.mean_sq(gg.as_double())**0.5,
        xray=scale*flex.mean_sq(O.g)**0.5,
        real_or_xray="xray")
      xray.minimization.add_gradients(
        scatterers=xs.scatterers(),
        xray_gradients=O.g,
        site_gradients=gg)
      O.f *= scale
      O.g *= scale
      O.last_grms.total = flex.mean_sq(O.g)**0.5
      O.g = flex.vec3_double(O.g)
    return O.f

  def d_e_pot_d_sites(O, sites_moved):
    O.e_pot(sites_moved=sites_moved)
    return matrix.col_list(O.g)

def run(fmodels, model, target_weights, params, log):
  assert fmodels.fmodel_neutron() is None # not implemented
  assert model.ias_selection is None # tardy+ias is not a useful combination
  xs = fmodels.fmodel_xray().xray_structure
  sites_cart_start = xs.sites_cart()
  sites = matrix.col_list(sites_cart_start)
  labels = [sc.label for sc in xs.scatterers()]
  tt = model.restraints_manager.geometry.construct_tardy_tree(
    sites=sites,
    omit_bonds_with_slack_greater_than
      =params.omit_bonds_with_slack_greater_than,
    constrain_dihedrals_with_sigma_less_than
      =params.constrain_dihedrals_with_sigma_less_than,
    near_singular_hinges_angular_tolerance_deg
      =params.near_singular_hinges_angular_tolerance_deg)
  print >> log, "tardy_tree summary:"
  tt.show_summary(vertex_labels=labels, out=log, prefix="  ")
  print >> log
  log.flush()
  potential_obj = potential_object(
    xray_weight_factor=params.xray_weight_factor,
    nonbonded_attenuation_factor=params.nonbonded_attenuation_factor,
    fmodels=fmodels,
    model=model,
    target_weights=target_weights,
    reduced_geo_manager=model.restraints_manager.geometry
      .reduce_for_tardy(
        tardy_tree=tt,
        omit_bonds_with_slack_greater_than
          =params.omit_bonds_with_slack_greater_than))
  tardy_model = scitbx.rigid_body.essence.tardy.model(
    labels=labels,
    sites=sites,
    masses=xs.atomic_weights(),
    tardy_tree=tt,
    potential_obj=potential_obj)
  action(tardy_model=tardy_model, params=params, callback=None, log=log)

def action(tardy_model, params, callback, log):
  sites_cart_start = flex.vec3_double(tardy_model.sites_moved())
  qd_e_kin_scales = tardy_model.assign_random_velocities()
  cartesian_dof = sites_cart_start.size() * 3
  if   (params.temperature_degrees_of_freedom == "cartesian"):
    temperature_dof = cartesian_dof
  elif (params.temperature_degrees_of_freedom == "constrained"):
    temperature_dof = tardy_model.degrees_of_freedom
  else:
    raise RuntimeError(
      "Unknown temperature_degrees_of_freedom: %s"
        % params.temperature_degrees_of_freedom)
  def e_as_t(e):
    return kinetic_energy_as_temperature(dof=temperature_dof, e=e)
  def t_as_e(t):
    return temperature_as_kinetic_energy(dof=temperature_dof, t=t)
  time_step_akma = params.time_step_pico_seconds / akma_time_as_pico_seconds
  print >> log, "tardy dynamics:"
  print >> log, "  number of bodies:", len(tardy_model.bodies)
  fmt = "%%%dd" % len(str(cartesian_dof))
  print >> log, "  number of degrees of freedom:", \
    fmt % tardy_model.degrees_of_freedom
  print >> log, "           number of atoms * 3:", fmt % cartesian_dof
  print >> log, "                         ratio: %.3f = 1/%.2f" % (
    tardy_model.degrees_of_freedom / max(1, cartesian_dof),
    cartesian_dof / max(1, tardy_model.degrees_of_freedom))
  print >> log, "  temperature degrees of freedom: %s (%d)" % (
    params.temperature_degrees_of_freedom, temperature_dof)
  print >> log, "  kinetic energy sensitivity to generalized velocities:"
  qd_e_kin_scales.min_max_mean().show(out=log, prefix="    ")
  print >> log, "  time step: %7.5f pico seconds" % (
    params.time_step_pico_seconds)
  print >> log, "  velocity scaling:", params.velocity_scaling
  log.flush()
  if (callback is not None):
    if (callback(tardy_model=tardy_model) == False): return
  n_time_steps = 0
  for i_cool_step in xrange(params.number_of_cooling_steps+1):
    if (params.number_of_cooling_steps == 0):
      if (   params.start_temperature_kelvin
          != params.final_temperature_kelvin):
        break
      t_target = params.start_temperature_kelvin
    else:
      t_target = params.start_temperature_kelvin \
               - i_cool_step * (  params.start_temperature_kelvin
                                - params.final_temperature_kelvin) \
                             / params.number_of_cooling_steps
    e_kin_target = t_as_e(t=t_target)
    def reset_e_kin(msg):
      tardy_model.reset_e_kin(e_kin_target=e_kin_target)
      print >> log, "  %s temperature: %8.2f K" % (
        msg, e_as_t(e=tardy_model.e_kin()))
      log.flush()
      return True
    if (   n_time_steps == 0
        or params.number_of_time_steps > 1):
      show_column_headings = reset_e_kin("new target")
    for i_time_step in xrange(params.number_of_time_steps):
      assert params.temperature_cap_factor > 1.0
      assert params.excessive_temperature_factor > params.temperature_cap_factor
      if (tardy_model.e_kin() > e_kin_target * params.temperature_cap_factor):
        print >> log, "  system temperature is too high:"
        if (tardy_model.e_kin()
              > e_kin_target * params.excessive_temperature_factor):
          print >> log, "    excessive_temperature_factor: %.6g" % \
            params.excessive_temperature_factor
          print >> log, "    excessive temperature limit: %.2f K" % (
            t_target * params.excessive_temperature_factor)
          print >> log, "    time_step_pico_seconds: %.6g" % (
            params.time_step_pico_seconds)
          log.flush()
          raise Sorry(
            "Excessive system temperature in torsion angle dynamics:\n"
            "  Please try again with a smaller time_step_pico_seconds.")
        print >> log, "     temperature_cap_factor: %.6g" % \
          params.temperature_cap_factor
        print >> log, "     temperature cap: %.2f K" % (
          t_target * params.temperature_cap_factor)
        show_column_headings = reset_e_kin("resetting")
      e_kin_before, e_tot_before = tardy_model.e_kin(), tardy_model.e_tot()
      tardy_model.dynamics_step(delta_t=time_step_akma)
      e_kin_after, e_tot_after = tardy_model.e_kin(), tardy_model.e_tot()
      if (e_tot_before > e_tot_after):
        fluct_e_tot = -e_tot_after / e_tot_before
      elif (e_tot_after > e_tot_before):
        fluct_e_tot = e_tot_before / e_tot_after
      else:
        fluct_e_tot = None
      if (params.velocity_scaling):
        tardy_model.reset_e_kin(e_kin_target=e_kin_target)
      n_time_steps += 1
      grms = tardy_model.potential_obj.last_grms
      if (show_column_headings):
        show_column_headings = False
        log.write("""\
          coordinate                   fluctuations        gradient rms
    step      rmsd        temp        temp   e_total     geo %    7s   total
""" % grms.real_or_xray)
      print >> log, "    %4d  %8.4f A  %8.2f K  %8.2f K  %s" \
        "  %6.2f  %6.2f  %6.2f" % (
          n_time_steps,
          flex.vec3_double(tardy_model.sites_moved()).rms_difference(
            sites_cart_start),
          e_as_t(e=tardy_model.e_kin()),
          e_as_t(e=e_kin_after-e_kin_before),
          format_value(format="%6.3f", value=fluct_e_tot),
          grms.geo,
          getattr(grms, grms.real_or_xray),
          grms.total)
      log.flush()
      if (callback is not None):
        if (callback() == False): return
  if (   params.minimization_max_iterations is None
      or params.minimization_max_iterations > 0):
    print >> log, "tardy gradient-driven minimization:"
    log.flush()
    def show_rms(minimizer=None):
      print >> log, "  coor. rmsd: %8.4f" % (
        flex.vec3_double(tardy_model.sites_moved()).rms_difference(
          sites_cart_start))
      log.flush()
      if (callback is not None):
        if (callback() == False): raise StopIteration
    try:
      refinery = tardy_model.minimization(
        max_iterations=params.minimization_max_iterations,
        callback_after_step=show_rms)
    except StopIteration:
      return
    print >> log, "After tardy minimization:"
    show_rms()
    print >> log, "  number of function evaluations:", \
      refinery.function_evaluations_total
    print >> log, "  number of lbfgs steps:", refinery.lbfgs_steps_total
    print >> log, "  number of lbfgs restarts:", refinery.lbfgs_restarts
    print >> log
