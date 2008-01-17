from cctbx.array_family import flex
import random
import time, math
from iotbx import pdb
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from scitbx import matrix
import scitbx.math
from cctbx import crystal
import iotbx.pdb
from cctbx import xray
from mmtbx import dynamics
from mmtbx.refinement import print_statistics
from cctbx import miller
import cctbx.xray.structure_factors
from mmtbx import bulk_solvent
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from mmtbx import dynamics
from mmtbx.dynamics import cartesian_dynamics
import mmtbx.refinement.minimization


def manager(simulated_annealing_params,
            bulk_solvent_parameters,
            refinement_parameters,
            alpha_beta_parameters,
            mask_parameters,
            target_weights,
            macro_cycle,
            tan_b_iso_max,
            h_params,
            fmodels,
            model,
            out = None):
  if(out is None): out = sys.stdout
  print_statistics.make_header("simulated annealing refinement", out = out)
  print_statistics.make_sub_header(
                   "lbfgs minimization: before simulated annealing", out = out)
  model.set_refine_individual_sites()
  minimized = mmtbx.refinement.minimization.lbfgs(
    restraints_manager       = model.restraints_manager,
    refine_xyz               = True,
    fmodels                  = fmodels,
    model                    = model,
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
              max_iterations = refinement_parameters.max_number_of_iterations),
    target_weights           = target_weights,
    h_params                 = h_params,
    verbose                  = 0)
  fmodel = fmodels.fmodel_xray() # XXX use only xray data
  minimized.monitor.show(message = "LBFGS minimization", log  = out)
  fmodel.update_xray_structure(xray_structure           = model.xray_structure,
                               update_f_calc            = True,
                               update_f_mask            = True,
                               update_f_ordered_solvent = False,
                               out                      = out)


  print_statistics.make_header("simulated annealing", out = out)
  wx = target_weights.xyz_weights_result.wx * \
    target_weights.xyz_weights_result.wx_scale
  run_simulated_annealing(
    simulated_annealing_params = simulated_annealing_params,
    model                      = model,
    fmodel                     = fmodel,
    wx                         = wx, # XXX
    bulk_solvent_parameters    = bulk_solvent_parameters,
    alpha_beta_parameters      = alpha_beta_parameters,
    mask_parameters            = mask_parameters,
    wc                         = target_weights.xyz_weights_result.w,
    out                        = out)

def run_simulated_annealing(simulated_annealing_params,
                            model,
                            fmodel,
                            wx,
                            wc,
                            bulk_solvent_parameters,
                            alpha_beta_parameters,
                            mask_parameters,
                            out):
  xray_structure_start        = model.xray_structure.deep_copy_scatterers()
  xray_structure_last_updated = model.xray_structure.deep_copy_scatterers()
  sa_temp = simulated_annealing_params.start_temperature
  xray_gradient = None
  reset_velocities = True
  while simulated_annealing_params.final_temperature <= sa_temp:
    print >> out
    cd_manager = cartesian_dynamics.cartesian_dynamics(
      structure                   = model.xray_structure,
      restraints_manager          = model.restraints_manager,
      temperature                 = sa_temp,
      n_steps                     = simulated_annealing_params.number_of_steps,
      time_step                   = simulated_annealing_params.time_step,
      initial_velocities_zero_fraction \
        = simulated_annealing_params.initial_velocities_zero_fraction,
      interleaved_minimization_params \
        = simulated_annealing_params.interleaved_minimization,
      n_print                     = simulated_annealing_params.n_print,
      fmodel                      = fmodel,
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
                                        update_f_calc            = True,
                                        update_f_mask            = True,
                                        update_f_ordered_solvent = False,
                                        out                      = out)
    fmodel.info().show_rfactors_targets_scales_overall(
      header = "2:SA temperature = "+str(sa_temp), out = out)

    geom_stat = model.show_geometry_statistics(
      message = "SA temperature = "+str(sa_temp))
    sa_temp -= simulated_annealing_params.cool_rate
