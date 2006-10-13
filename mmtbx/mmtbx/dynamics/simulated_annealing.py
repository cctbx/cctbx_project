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
            wxnc_scale,
            tan_b_iso_max,
            monitor,
            fmodel,
            model,
            out = None):
  if(out is None): out = sys.stdout
  print_statistics.make_header("simulated annealing refinement", out = out)
  print_statistics.make_sub_header(
                   "lbfgs minimization: before simulated annealing", out = out)
  minimized = mmtbx.refinement.minimization.lbfgs(
    restraints_manager       = model.restraints_manager,
    refine_xyz               = True,
    fmodel                   = fmodel,
    model                    = model,
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
              max_iterations = refinement_parameters.max_number_of_iterations),
    wx                          = target_weights.wx_xyz(),
    wc                          = target_weights.wc(),
    wxnc_scale                  = wxnc_scale,
    verbose                     = 0)

  minimized.collector.show(text = "lbfgs minimization", out = out)
  fmodel.update_xray_structure(xray_structure           = model.xray_structure,
                               update_f_calc            = True,
                               update_f_mask            = False,
                               update_f_ordered_solvent = False,
                               out                      = out)

  monitor.collect(step           = str(macro_cycle) + "_xyz:",
                  model          = model,
                  fmodel         = fmodel,
                  tan_b_iso_max  = tan_b_iso_max,
                  target_weights = target_weights,
                  wilson_b       = None)

  print_statistics.make_header("simulated annealing", out = out)
  run_simulated_annealing(
                       simulated_annealing_params = simulated_annealing_params,
                       model                      = model,
                       fmodel                     = fmodel,
                       wx                         = target_weights.wx_xyz(),
                       bulk_solvent_parameters    = bulk_solvent_parameters,
                       alpha_beta_parameters      = alpha_beta_parameters,
                       mask_parameters            = mask_parameters,
                       wc                         = target_weights.wc(),
                       out                        = out)

  fmodel.update_xray_structure(xray_structure           = model.xray_structure,
                               update_f_calc            = True,
                               update_f_mask            = False,
                               update_f_ordered_solvent = False,
                               out                      = out)

  monitor.collect(step           = str(macro_cycle) + "_sar:",
                  model          = model,
                  fmodel         = fmodel,
                  tan_b_iso_max  = tan_b_iso_max,
                  target_weights = target_weights,
                  wilson_b       = None)

def run_simulated_annealing(simulated_annealing_params,
                            model,
                            fmodel,
                            wx,
                            wc,
                            bulk_solvent_parameters,
                            alpha_beta_parameters,
                            mask_parameters,
                            out):
  assert fmodel.sf_algorithm is not None
  sf_algorithm = fmodel.sf_algorithm
  fmodel_copy = fmodel.deep_copy()
  fmodel_copy_1 = fmodel.deep_copy()
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
      n_print                     = simulated_annealing_params.n_print,
      verbose                     = simulated_annealing_params.verbose,
      fmodel                      = fmodel_copy,
      xray_target_weight          = wx,
      chem_target_weight          = wc,
      xray_structure_last_updated = xray_structure_last_updated,
      shift_update                = simulated_annealing_params.update_grads_shift,
      xray_gradient               = xray_gradient,
      reset_velocities            = reset_velocities)
    reset_velocities = False

    xray_structure_last_updated = \
                  cd_manager.xray_structure_last_updated.deep_copy_scatterers()
    xray_gradient = cd_manager.xray_gradient

    fmodel_copy_1.update_xray_structure(xray_structure  = model.xray_structure,
                                        update_f_calc            = True,
                                        update_f_mask            = False,
                                        update_f_ordered_solvent = False,
                                        out                      = out)

    fmodel_copy_1.show_essential(header = "2:SA temperatrure = "+str(sa_temp),
                                 out    = out)

    geom_stat = model.geometry_statistics(
                                      show = True,
                                      text = "SA temperatrure = "+str(sa_temp))
    sa_temp -= simulated_annealing_params.cool_rate
