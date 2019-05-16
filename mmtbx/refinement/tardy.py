from __future__ import absolute_import, division, print_function
from mmtbx.dynamics import \
  kinetic_energy_as_temperature, \
  temperature_as_kinetic_energy
from mmtbx.dynamics.constants import akma_time_as_pico_seconds
from cctbx import xray
import cctbx.geometry_restraints
from cctbx.array_family import flex
import scitbx.rigid_body
from scitbx.graph import tardy_tree
from scitbx import matrix
from libtbx.utils import Sorry
from libtbx.str_utils import format_value, show_string
from libtbx import group_args
import os
from six.moves import range
op = os.path

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
  prolsq_repulsion_function_changes
    .help = "energy(delta) = "
            "c_rep*(max(0,(k_rep*vdw_distance)**irexp-delta**irexp))**rexp"
  {
    c_rep = None
      .type = float
      .help = "Usual value: 16"
    k_rep = 0.75
      .type = float
      .help = "Usual value: 1."
              "Smaller values reduce the distance threshold at which"
              "the repulsive force becomes active."
    irexp = None
      .type = float
      .help = "Usual value: 1"
    rexp = None
      .type = float
      .help = "Usual value: 4"
  }
  omit_bonds_with_slack_greater_than = 0
    .type = float
  constrain_dihedrals_with_sigma_less_than = 10
    .type = float
  near_singular_hinges_angular_tolerance_deg = 5
    .type = float
  emulate_cartesian = False
    .type = bool
  trajectory_directory = None
    .type = path
    .expert_level = 3
"""

class potential_object(object):

  def __init__(O,
        xray_weight_factor,
        prolsq_repulsion_function_changes,
        fmodels,
        model,
        target_weights,
        reduced_geo_manager):
    O.xray_weight_factor = xray_weight_factor
    O.fmodels = fmodels
    O.model = model
    O.weights = target_weights.xyz_weights_result
    O.reduced_geo_manager = reduced_geo_manager
    c = prolsq_repulsion_function_changes
    if (    c.c_rep is None
        and c.k_rep is None
        and c.irexp is None
        and c.rexp is None):
      O.custom_nonbonded_function = None
    else:
      nonbonded_function = model.restraints_manager.geometry.nonbonded_function
      assert isinstance(
        nonbonded_function,
        cctbx.geometry_restraints.prolsq_repulsion_function)
      O.custom_nonbonded_function = nonbonded_function.customized_copy(
        c_rep=c.c_rep,
        k_rep=c.k_rep,
        irexp=c.irexp,
        rexp=c.rexp)
    O.fmodels.create_target_functors()
    O.fmodels.prepare_target_functors_for_minimization()
    O.allowed_origin_shifts_need_to_be_suppressed = \
      O.fmodels.target_functions_are_invariant_under_allowed_origin_shifts()
    O.last_sites_moved = None
    O.f = None
    O.g = None
    O.last_grms = None

  def crystal_symmetry(O):
    return O.model.restraints_manager.geometry.crystal_symmetry

  def e_pot(O, sites_moved):
    if (   O.last_sites_moved is None
        or O.last_sites_moved.id() != sites_moved.id()):
      O.last_sites_moved = sites_moved
      xs = O.fmodels.fmodel_xray().xray_structure
      assert len(sites_moved) == xs.scatterers().size()
      sites_cart = sites_moved
      xs.set_sites_cart(sites_cart=sites_cart)
      O.fmodels.update_xray_structure(update_f_calc=True)
      xs.scatterers().flags_set_grads(state=False)
      xs.scatterers().flags_set_grad_site(
        iselection=flex.size_t_range(xs.scatterers().size()))
      expected_g_size = len(sites_moved) * 3
      if (O.xray_weight_factor is not None):
        tg = O.fmodels.target_and_gradients(
          weights=O.weights,
          compute_gradients=True)
        O.f = tg.target() * O.xray_weight_factor
        O.g = tg.gradients() * O.xray_weight_factor
        assert O.g.size() == expected_g_size
      else:
        O.f = 0.
        O.g = flex.double(expected_g_size, 0)
      if (O.reduced_geo_manager is None):
        reduced_geo_energies = None
      else:
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
      gg = other_energies.gradients * O.weights.w
      if (reduced_geo_energies is not None):
        O.f += reduced_geo_energies.target * nfw
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
    return O.g

def run(fmodels, model, target_weights, params, log,
        format_for_phenix_refine=False, monitor=None,
        call_back_after_step=True):
  assert fmodels.fmodel_neutron() is None # not implemented
  assert model.ias_manager is None # tardy+ias is not a useful combination
  xs = fmodels.fmodel_xray().xray_structure
  sites_cart_start = xs.sites_cart()
  sites = sites_cart_start
  labels = [sc.label for sc in xs.scatterers()]
  if (params.emulate_cartesian):
    tt = scitbx.graph.tardy_tree.construct(
      n_vertices=len(sites), edge_list=[])
    tt.build_tree()
  else:
    tt = model.restraints_manager.geometry.construct_tardy_tree(
      sites=sites,
      selection=model.refinement_flags.sites_torsion_angles,
      omit_bonds_with_slack_greater_than
        =params.omit_bonds_with_slack_greater_than,
      constrain_dihedrals_with_sigma_less_than
        =params.constrain_dihedrals_with_sigma_less_than,
      near_singular_hinges_angular_tolerance_deg
        =params.near_singular_hinges_angular_tolerance_deg)
  print("tardy_tree summary:", file=log)
  tt.show_summary(vertex_labels=labels, out=log, prefix="  ")
  print(file=log)
  log.flush()
  if (params.emulate_cartesian):
    reduced_geo_manager = None
  else:
    include_den_restraints = False
    if model.restraints_manager.geometry.den_manager is not None:
      if "torsion" in model.restraints_manager.geometry. \
         den_manager.params.annealing_type:
        include_den_restraints = True
    reduced_geo_manager = model.restraints_manager.geometry \
      .reduce_for_tardy(
        tardy_tree=tt,
        omit_bonds_with_slack_greater_than
          =params.omit_bonds_with_slack_greater_than,
        include_den_restraints=include_den_restraints)
  potential_obj = potential_object(
    xray_weight_factor=params.xray_weight_factor,
    prolsq_repulsion_function_changes=params.prolsq_repulsion_function_changes,
    fmodels=fmodels,
    model=model,
    target_weights=target_weights,
    reduced_geo_manager=reduced_geo_manager)
  tardy_model = scitbx.rigid_body.tardy_model(
    labels=labels,
    sites=sites,
    masses=xs.atomic_weights(),
    tardy_tree=tt,
    potential_obj=potential_obj,
    near_singular_hinges_angular_tolerance_deg=
      params.near_singular_hinges_angular_tolerance_deg)
  def refinement_callback(fmodel):
    if (monitor is not None) and (call_back_after_step):
      monitor.call_back(model, fmodel, "torsion_dynamics")
  action( # XXX neutron
    tardy_model=tardy_model,
    params=params,
    rmsd_calculator=tt.rmsd_calculator(),
    callback=None,
    log=log,
    fmodel=fmodels.fmodel_xray(),
    format_for_phenix_refine=format_for_phenix_refine,
    refinement_callback=refinement_callback)

def action(
      tardy_model,
      params,
      rmsd_calculator,
      callback,
      log,
      fmodel=None,
      format_for_phenix_refine=False,
      refinement_callback=None):
  sites_cart_start = tardy_model.sites_moved()
  qd_e_kin_scales = tardy_model.assign_random_velocities()
  traj_dir = params.trajectory_directory
  if (traj_dir is not None):
    print("Creating trajectory directory: %s" % show_string(traj_dir), file=log)
    from libtbx.path import move_old_create_new_directory
    move_old_create_new_directory(path=traj_dir)
    from libtbx import easy_pickle
    easy_pickle.dump(
      file_name=op.join(traj_dir, "labels"),
      obj=tardy_model.labels)
    traj_serial_fmt = "sites_%%0%dd" % len(
      "%d" % (params.number_of_cooling_steps * params.number_of_time_steps))
    easy_pickle.dump(
      file_name=op.join(traj_dir, traj_serial_fmt % 0),
      obj=sites_cart_start)
    if (fmodel is not None):
      easy_pickle.dump(
        file_name=op.join(traj_dir, "xray_structure"),
        obj=fmodel.xray_structure)
    print(file=log)
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
  print("tardy dynamics:", file=log)
  print("  number of bodies:", tardy_model.bodies_size(), file=log)
  fmt = "%%%dd" % len(str(cartesian_dof))
  print("  number of degrees of freedom:", \
    fmt % tardy_model.degrees_of_freedom, file=log)
  print("           number of atoms * 3:", fmt % cartesian_dof, file=log)
  print("                         ratio: %.3f = 1/%.2f" % (
    tardy_model.degrees_of_freedom / max(1, cartesian_dof),
    cartesian_dof / max(1, tardy_model.degrees_of_freedom)), file=log)
  print("  temperature degrees of freedom: %s (%d)" % (
    params.temperature_degrees_of_freedom, temperature_dof), file=log)
  print("  kinetic energy sensitivity to generalized velocities:", file=log)
  qd_e_kin_scales.min_max_mean().show(out=log, prefix="    ")
  print("  time step: %7.5f pico seconds" % (
    params.time_step_pico_seconds), file=log)
  print("  velocity scaling:", params.velocity_scaling, file=log)
  allowed_origin_shifts_need_to_be_suppressed = tardy_model.potential_obj \
    .allowed_origin_shifts_need_to_be_suppressed
  print("  suppressing allowed origin shifts:", \
    allowed_origin_shifts_need_to_be_suppressed, file=log)
  log.flush()
  if (callback is not None):
    if (callback(
          tardy_model=tardy_model,
          rmsd_calculator=rmsd_calculator) == False):
      return
  if (allowed_origin_shifts_need_to_be_suppressed):
    number_of_sites_in_each_tree = tardy_model.number_of_sites_in_each_tree()
    crystal_symmetry = tardy_model.potential_obj.crystal_symmetry()
    assert crystal_symmetry is not None
    allowed_origin_shift_velocity_corrections = flex.double()
    sum_of_allowed_origin_shift_velocity_corrections = [matrix.col((0,0,0))]
  def suppress_allowed_origin_shifts(collect_stats):
    if (not allowed_origin_shifts_need_to_be_suppressed):
      return
    mlv = tardy_model.mean_linear_velocity(
      number_of_sites_in_each_tree=number_of_sites_in_each_tree)
    if (mlv is not None):
      mlv = matrix.col(mlv)
      mlv_perp = matrix.col(
        crystal_symmetry.subtract_continuous_allowed_origin_shifts(
          translation_cart=mlv))
      correction = mlv - mlv_perp
      tardy_model.subtract_from_linear_velocities(
        number_of_sites_in_each_tree=number_of_sites_in_each_tree,
        value=correction)
      if (collect_stats):
        allowed_origin_shift_velocity_corrections.append(abs(correction))
        sum_of_allowed_origin_shift_velocity_corrections[0] += correction
  suppress_allowed_origin_shifts(collect_stats=False)
  n_time_steps = 0
  den_update_interval = None
  if (tardy_model.potential_obj.reduced_geo_manager is not None):
    den_update_interval = (50 *params.number_of_cooling_steps /
        (params.start_temperature_kelvin - params.final_temperature_kelvin))
    den_update_interval = int(round(den_update_interval))
  for i_cool_step in range(params.number_of_cooling_steps+1):
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
    if (tardy_model.potential_obj.reduced_geo_manager is not None):
      if (tardy_model.potential_obj.reduced_geo_manager. \
          den_manager is not None) and \
         i_cool_step > 0 and not (i_cool_step)%den_update_interval:
        print("   update DEN eq distances at step %d, temp=%.1f" % \
          ( (n_time_steps), t_target), file=log)
        tardy_model.potential_obj.reduced_geo_manager. \
          den_manager.update_eq_distances(
            sites_cart=tardy_model.sites_moved())
    e_kin_target = t_as_e(t=t_target)
    def reset_e_kin(msg):
      tardy_model.reset_e_kin(e_kin_target=e_kin_target)
      print("  %s temperature: %8.2f K" % (
        msg, e_as_t(e=tardy_model.e_kin())), file=log)
      log.flush()
      return True
    if (   n_time_steps == 0
        or params.number_of_time_steps > 1):
      show_column_headings = reset_e_kin("new target")
    for i_time_step in range(params.number_of_time_steps):
      assert params.temperature_cap_factor > 1.0
      assert params.excessive_temperature_factor > params.temperature_cap_factor
      if (tardy_model.e_kin() > e_kin_target * params.temperature_cap_factor):
        print("  system temperature is too high:", file=log)
        if (tardy_model.e_kin()
              > e_kin_target * params.excessive_temperature_factor):
          print("    excessive_temperature_factor: %.6g" % \
            params.excessive_temperature_factor, file=log)
          print("    excessive temperature limit: %.2f K" % (
            t_target * params.excessive_temperature_factor), file=log)
          print("    time_step_pico_seconds: %.6g" % (
            params.time_step_pico_seconds), file=log)
          log.flush()
          raise Sorry(
            "Excessive system temperature in torsion angle dynamics:\n"
            "  Please try again with a smaller time_step_pico_seconds.")
        print("     temperature_cap_factor: %.6g" % \
          params.temperature_cap_factor, file=log)
        print("     temperature cap: %.2f K" % (
          t_target * params.temperature_cap_factor), file=log)
        show_column_headings = reset_e_kin("resetting")
      e_kin_before, e_tot_before = tardy_model.e_kin(), tardy_model.e_tot()
      tardy_model.dynamics_step(delta_t=time_step_akma)
      suppress_allowed_origin_shifts(collect_stats=True)
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
      if (traj_dir is not None):
        easy_pickle.dump(
          file_name=op.join(traj_dir, traj_serial_fmt % n_time_steps),
          obj=tardy_model.sites_moved())
      grms = tardy_model.potential_obj.last_grms
      if(format_for_phenix_refine):
        if(n_time_steps==1 or not n_time_steps%25):
          fmtr = "   step=%s temperature=%s rmsd=%s r_work=%s r_free=%s"
          rmsd = rmsd_calculator(tardy_model.sites_moved(), sites_cart_start)
          print(fmtr%(
            format_value("%5d", n_time_steps),
            format_value("%7.1f", e_as_t(e=tardy_model.e_kin())),
            format_value("%6.4f", rmsd),
            format_value("%6.4f", fmodel.r_work()),
            format_value("%6.4f", fmodel.r_free())), file=log)
          if hasattr(refinement_callback, "__call__"):
            refinement_callback(fmodel) # for phenix.refine GUI
      else:
        if (show_column_headings):
          show_column_headings = False
          log.write("""\
            coordinate                   fluctuations        gradient rms
    step      rmsd        temp        temp   e_total     geo %    7s   total
""" % grms.real_or_xray)
        print("    %4d  %8.4f A  %8.2f K  %8.2f K  %s" \
          "  %6.2f  %6.2f  %6.2f" % (
            n_time_steps,
            rmsd_calculator(tardy_model.sites_moved(), sites_cart_start),
            e_as_t(e=tardy_model.e_kin()),
            e_as_t(e=e_kin_after-e_kin_before),
            format_value(format="%6.3f", value=fluct_e_tot),
            grms.geo,
            getattr(grms, grms.real_or_xray),
            grms.total), file=log)
      log.flush()
      if (callback is not None):
        if (callback() == False): return
  if (allowed_origin_shifts_need_to_be_suppressed):
    print("  allowed origin shift velocity corrections applied (magnitudes):", file=log)
    allowed_origin_shift_velocity_corrections.min_max_mean().show(
      out=log, prefix="    ")
    print("  sum of allowed origin shift velocity corrections (vectors):", file=log)
    print("   ", numstr(
      values=sum_of_allowed_origin_shift_velocity_corrections[0].elems,
      fmt="%.5f",
      brackets=("(",")")), file=log)
  if (   params.minimization_max_iterations is None
      or params.minimization_max_iterations > 0):
    print("tardy gradient-driven minimization:", file=log)
    log.flush()
    def show_rms(minimizer=None):
      print("  coor. rmsd: %8.4f" % (
        rmsd_calculator(tardy_model.sites_moved(), sites_cart_start)), file=log)
      log.flush()
      if (callback is not None):
        if (callback() == False): raise StopIteration
    try:
      refinery = tardy_model.minimization(
        max_iterations=params.minimization_max_iterations,
        callback_after_step=show_rms)
    except StopIteration:
      return
    print("After tardy minimization:", file=log)
    show_rms()
    print("  number of function evaluations:", \
      refinery.function_evaluations_total, file=log)
    print("  number of lbfgs steps:", refinery.lbfgs_steps_total, file=log)
    print("  number of lbfgs restarts:", refinery.lbfgs_restarts, file=log)
    print(file=log)
