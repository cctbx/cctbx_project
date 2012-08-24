from __future__ import division
from cctbx import xray
from mmtbx.refinement import print_statistics
from mmtbx.dynamics import cartesian_dynamics
import mmtbx.refinement.minimization
from scitbx.array_family import flex
import iotbx.phil
import sys

master_params_str = """\
start_temperature = 5000
  .type = float
final_temperature = 300
  .type = float
cool_rate = 100
  .type = float
number_of_steps = 50
  .type = int
time_step = 0.0005
  .type = float
  .expert_level=2
initial_velocities_zero_fraction = 0
  .type = float
  .expert_level=2
interleaved_minimization
  .expert_level=2
  .style = box
{
  number_of_iterations = 0
    .type = float
  time_step_factor = 10
    .type = float
  restraints = *bonds *angles
    .type = choice(multi=True)
}
n_print = 100
  .type = int
  .short_caption = Steps between log output
  .expert_level=2
update_grads_shift = 0.3
  .type = float
  .short_caption = Update gradient shifts
  .expert_level=2
refine_sites = True
  .caption = "lbfgs refinement of atomic coordinates before sa"
  .short_caption=Refine sites first
  .type = bool
refine_adp = False
  .caption = "lbfgs refinement of adp before sa"
  .short_caption=Refine ADPs first
  .type = bool
max_number_of_iterations = 25
  .type = int
mode = every_macro_cycle *second_and_before_last once first first_half
  .type = choice
verbose = -1
  .type = int
  .short_caption = Verbosity level
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def manager(params,
            target_weights,
            all_params,
            macro_cycle,
            h_params,
            fmodels,
            model,
            out = None):
  if(out is None): out = sys.stdout
  print_statistics.make_header("simulated annealing refinement", out = out)
  model.set_refine_individual_sites()
  fmodel = fmodels.fmodel_xray() # XXX use only xray data
  fmodel.xray_structure = model.xray_structure # XXX use only xray data
  if (params.max_number_of_iterations >= 0):
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
        max_iterations = params.max_number_of_iterations),
      target_weights           = target_weights,
      h_params                 = h_params,
      verbose                  = 0)
  fmodel.update_xray_structure(xray_structure = model.xray_structure,
                               update_f_calc  = True,
                               update_f_mask  = True)
  print_statistics.make_header("simulated annealing", out = out)
  wx = target_weights.xyz_weights_result.wx * \
    target_weights.xyz_weights_result.wx_scale
  run_simulated_annealing(
    params = params,
    restraints_manager = model.restraints_manager,
    fmodel             = fmodel,
    wx                 = wx,
    wc                 = target_weights.xyz_weights_result.w,
    out                = out)

def run_simulated_annealing(params,
                            fmodel,
                            restraints_manager,
                            wx,
                            wc,
                            out = None,
                            verbose=True):
  if(out is None): out = sys.stdout
  sites_cart_start = fmodel.xray_structure.sites_cart()
  sa_temp = params.start_temperature
  verbose = params.verbose
  reset_velocities = True
  vxyz = None
  cd_manager = None
  den_manager = getattr(restraints_manager.geometry.generic_restraints_manager,
    "den_manager", None)
  cartesian_den_restraints = False
  if(den_manager is not None):
    if("cartesian" in den_manager.params.annealing_type):
      restraints_manager.geometry.generic_restraints_manager.flags.den = True
      cartesian_den_restraints = True
      verbose = False
  if(verbose):
    print >> out, "  sa_temp r_work r_free distance_moved rmsd_bond rmsd_angle"
  while params.final_temperature <= sa_temp:
    if(sa_temp==params.start_temperature):
      cmremove=True
    else: cmremove=False
    #
    from mmtbx.dynamics import cartesian_dynamics as cd
    gradients_calculator = cd.gradients_calculator_reciprocal_space(
      restraints_manager        = restraints_manager,
      fmodel                    = fmodel,
      sites_cart                = fmodel.xray_structure.sites_cart(),
      wx                        = wx,
      wc                        = wc,
      update_gradient_threshold = params.update_grads_shift)
    #
    cd_manager = cartesian_dynamics.cartesian_dynamics(
      xray_structure              = fmodel.xray_structure,
      gradients_calculator        = gradients_calculator,
      temperature                 = sa_temp,
      vxyz                        = vxyz,
      n_steps                     = params.number_of_steps,
      time_step                   = params.time_step,
      initial_velocities_zero_fraction \
        = params.initial_velocities_zero_fraction,
      n_print                     = params.n_print,
      stop_cm_motion              = cmremove,
      reset_velocities            = reset_velocities,
      log=out,
      verbose=verbose)
    reset_velocities = False
    fmodel.update_xray_structure(
      xray_structure = fmodel.xray_structure,
      update_f_calc  = True,
      update_f_mask  = True)
    if(verbose):
      sites_cart = fmodel.xray_structure.sites_cart()
      es = restraints_manager.geometry.energies_sites(sites_cart = sites_cart)
      dist = flex.mean(flex.sqrt((sites_cart_start - sites_cart).dot()))
      fmt="  %7.1f %6.4f %6.4f         %6.2f    %6.3f     %6.2f"
      print >> out, fmt%(sa_temp, fmodel.r_work(), fmodel.r_free(), dist,
        es.bond_deviations()[2], es.angle_deviations()[2])
    sa_temp -= params.cool_rate
    if(cartesian_den_restraints):
      print >> out, "update DEN eq distances at temp=%.1f" % sa_temp
      den_manager.update_eq_distances(
        sites_cart=fmodel.xray_structure.sites_cart())
  if(den_manager is not None):
    restraints_manager.geometry.generic_restraints_manager.flags.den = False
