from cctbx import xray
from mmtbx.refinement import print_statistics
from mmtbx.dynamics import cartesian_dynamics
import mmtbx.refinement.minimization


def manager(simulated_annealing_params,
            bulk_solvent_parameters,
            refinement_parameters,
            alpha_beta_parameters,
            mask_parameters,
            target_weights,
            all_params,
            macro_cycle,
            tan_b_iso_max,
            h_params,
            fmodels,
            model,
            out = None,
            monitor=None):
  if(out is None): out = sys.stdout
  print_statistics.make_header("simulated annealing refinement", out = out)
  model.set_refine_individual_sites()
  fmodel = fmodels.fmodel_xray() # XXX use only xray data
  fmodel.xray_structure = model.xray_structure # XXX use only xray data
  if (simulated_annealing_params.max_number_of_iterations >= 0):
    print_statistics.make_sub_header(
      "lbfgs minimization: before simulated annealing", out = out)
    is_neutron_scat_table = False
    if(all_params.main.scattering_table == "neutron"):
      is_neutron_scat_table = True
    import scitbx.lbfgs
    minimized = mmtbx.refinement.minimization.lbfgs(
      restraints_manager       = model.restraints_manager,
      refine_xyz               = True,
      fmodels                  = fmodels,
      is_neutron_scat_table    = is_neutron_scat_table,
      model                    = model,
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = simulated_annealing_params.max_number_of_iterations),
      target_weights           = target_weights,
      h_params                 = h_params,
      verbose                  = 0)
    minimized.monitor.show(message = "LBFGS minimization", log  = out)
  fmodel.update_xray_structure(xray_structure = model.xray_structure,
                               update_f_calc  = True,
                               update_f_mask  = True,
                               out            = out)
  print_statistics.make_header("simulated annealing", out = out)
  wx = target_weights.xyz_weights_result.wx * \
    target_weights.xyz_weights_result.wx_scale
  run_simulated_annealing(
    simulated_annealing_params = simulated_annealing_params,
    model                      = model,
    fmodel                     = fmodel,
    wx                         = wx, # XXX
    neutron_refinement         = fmodels.neutron_refinement,
    bulk_solvent_parameters    = bulk_solvent_parameters,
    alpha_beta_parameters      = alpha_beta_parameters,
    mask_parameters            = mask_parameters,
    wc                         = target_weights.xyz_weights_result.w,
    out                        = out,
    monitor                    = monitor)

def run_simulated_annealing(simulated_annealing_params,
                            model,
                            fmodel,
                            wx,
                            wc,
                            neutron_refinement,
                            bulk_solvent_parameters,
                            alpha_beta_parameters,
                            mask_parameters,
                            out,
                            monitor):
  xray_structure_last_updated = model.xray_structure.deep_copy_scatterers()
  sa_temp = simulated_annealing_params.start_temperature
  xray_gradient = None
  reset_velocities = True
  vxyz = None
  cd_manager = None
  while simulated_annealing_params.final_temperature <= sa_temp:
    print >> out
    if(sa_temp==simulated_annealing_params.start_temperature):
      cmremove=True
    else: cmremove=False
    cd_manager = cartesian_dynamics.cartesian_dynamics(
      structure                   = model.xray_structure,
      restraints_manager          = model.restraints_manager,
      temperature                 = sa_temp,
      vxyz                        = vxyz,
      n_steps                     = simulated_annealing_params.number_of_steps,
      time_step                   = simulated_annealing_params.time_step,
      initial_velocities_zero_fraction \
        = simulated_annealing_params.initial_velocities_zero_fraction,
      interleaved_minimization_params \
        = simulated_annealing_params.interleaved_minimization,
      n_print                     = simulated_annealing_params.n_print,
      fmodel                      = fmodel,
      stop_cm_motion              = cmremove,
      xray_target_weight          = wx,
      chem_target_weight          = wc,
      xray_structure_last_updated = xray_structure_last_updated,
      shift_update                = simulated_annealing_params.update_grads_shift,
      xray_gradient               = xray_gradient,
      reset_velocities            = reset_velocities,
      log=out,
      verbose=simulated_annealing_params.verbose)
    reset_velocities = False
    xray_structure_last_updated = \
                  cd_manager.xray_structure_last_updated.deep_copy_scatterers()
    xray_gradient = cd_manager.xray_gradient
    fmodel.update_xray_structure(xray_structure  = model.xray_structure,
                                        update_f_calc = True,
                                        update_f_mask = True,
                                        out           = out)
    fmodel.info().show_rfactors_targets_scales_overall(
      header = "2:SA temperature = "+str(sa_temp), out = out)

    geom_stat = model.show_geometry_statistics(
      ignore_hd = not neutron_refinement,
      message = "SA temperature = "+str(sa_temp))
    if monitor is not None :
      monitor.call_back(model, fmodel, "simulated_annealing")
    sa_temp -= simulated_annealing_params.cool_rate
