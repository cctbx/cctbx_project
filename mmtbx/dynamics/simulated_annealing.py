from __future__ import absolute_import, division, print_function
from cctbx import xray
from mmtbx.refinement import print_statistics
from mmtbx.dynamics import cartesian_dynamics
import mmtbx.refinement.minimization
from scitbx.array_family import flex
import iotbx.phil
import sys
from libtbx import adopt_init_args

main_params_str = """\
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
interleave_minimization = False
  .type = bool
verbose = -1
  .type = int
  .short_caption = Verbosity level
n_print = 100
  .type = int
  .short_caption = Steps between log output
  .expert_level=2
update_grads_shift = 0.3
  .type = float
  .short_caption = Update gradient shifts
  .expert_level=2
random_seed = None
  .type = int
"""

master_params_str = """\
%s
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
"""%main_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def manager(params,
            target_weights,
            all_params,
            macro_cycle,
            h_params,
            fmodels,
            model,
            out = None,
            states_collector=None,
            callback=None):
  if(out is None): out = sys.stdout
  if (states_collector is not None):
    assert hasattr(states_collector, "add")
  print_statistics.make_header("simulated annealing refinement", out = out)
  model.set_refine_individual_sites()
  fmodel = fmodels.fmodel_xray() # XXX use only xray data
  fmodel.xray_structure = model.get_xray_structure() # XXX use only xray data
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
  fmodel.update_xray_structure(xray_structure = model.get_xray_structure(),
                               update_f_calc  = True,
                               update_f_mask  = True)
  print_statistics.make_header("simulated annealing", out = out)
  wx = target_weights.xyz_weights_result.wx * \
    target_weights.xyz_weights_result.wx_scale
  run(
    params = params,
    restraints_manager = model.restraints_manager,
    fmodel             = fmodel,
    wx                 = wx,
    wc                 = target_weights.xyz_weights_result.w,
    states_collector   = states_collector,
    callback           = callback,
    log                = out)

class run(object):
  def __init__(self,
               params,
               restraints_manager,
               xray_structure = None,
               wx = None,
               wc = None,
               fmodel = None,
               target_map = None,
               real_space = False,
               log = None,
               states_collector = None,
               callback = None,
               verbose=True):
    adopt_init_args(self, locals())
    assert (callback is None) or hasattr(callback, "__call__")
    if(self.params is None): self.params = master_params().extract()
    if(log is None): self.log = sys.stdout
    if(self.fmodel is not None):
      self.xray_structure = self.fmodel.xray_structure
    self.sites_cart_start = self.xray_structure.sites_cart()
    self.curr_temp = params.start_temperature
    verbose = params.verbose
    reset_velocities = True
    vxyz = None
    den_manager = restraints_manager.geometry.den_manager
    cartesian_den_restraints = False
    if(den_manager is not None):
      if("cartesian" in den_manager.params.annealing_type):
        # restraints_manager.geometry.generic_restraints_manager.flags.den = True
        cartesian_den_restraints = True
        verbose = False
    while params.final_temperature <= self.curr_temp:
      #if(self.curr_temp == params.start_temperature):
      #  cmremove=True
      #else: cmremove=False
      cmremove=True
      cd_manager = cartesian_dynamics.run(
        xray_structure           = self.xray_structure,
        gradients_calculator     = self.gradients_calculator(),
        temperature              = self.curr_temp,
        interleaved_minimization = self.params.interleave_minimization,
        vxyz                 = vxyz,
        n_steps              = self.params.number_of_steps,
        time_step            = self.params.time_step,
        n_print              = self.params.n_print,
        random_seed          = self.params.random_seed,
        stop_cm_motion       = cmremove,
        reset_velocities     = reset_velocities,
        log                  = self.log,
        verbose              = verbose)
      reset_velocities = False
      vxyz = cd_manager.vxyz
      self.xray_structure = cd_manager.xray_structure
      if(self.fmodel is not None):
        self.fmodel.update_xray_structure(
          xray_structure = self.xray_structure,
          update_f_calc  = True,
          update_f_mask  = True)
      if(states_collector is not None):
        self.states_collector.add(sites_cart = cd_manager.xray_structure.sites_cart())
      if (callback is not None):
        callback(fmodel=self.fmodel)
      self.show(curr_temp = self.curr_temp)
      self.curr_temp -= params.cool_rate
      if(cartesian_den_restraints):
        print("update DEN eq distances at temp=%.1f" % self.curr_temp, file=self.log)
        den_manager.update_eq_distances(
          sites_cart=fmodel.xray_structure.sites_cart())
    # if(den_manager is not None):
    #   restraints_manager.geometry.generic_restraints_manager.flags.den = False

  def show(self, curr_temp):
    if(self.verbose):
      sites_cart = self.xray_structure.sites_cart()
      es=self.restraints_manager.geometry.energies_sites(sites_cart=sites_cart)
      a,b = es.bond_deviations()[2], es.angle_deviations()[2]
      dist = flex.mean(flex.sqrt((self.sites_cart_start - sites_cart).dot()))
      if(self.fmodel is not None):
        fmt="  temp=%7.1f r_work=%6.4f r_free=%6.4f dist_moved=%6.2f angles=%6.2f bonds=%6.3f"
        print(fmt%(curr_temp, self.fmodel.r_work(),
          self.fmodel.r_free(), dist, b, a), file=self.log)
      else:
        fmt="  temp=%7.1f dist_moved=%6.2f angles=%6.2f bonds=%6.3f"
        print(fmt%(curr_temp, dist, b, a), file=self.log)

  def gradients_calculator(self):
    if(not self.real_space):
      if(self.fmodel is not None):
        grad_calc = cartesian_dynamics.gradients_calculator_reciprocal_space(
          restraints_manager        = self.restraints_manager, # XXX WHY?
          fmodel                    = self.fmodel,
          sites_cart                = self.fmodel.xray_structure.sites_cart(),
          wx                        = self.wx,
          wc                        = self.wc,
          update_gradient_threshold = self.params.update_grads_shift)
      else:
        grad_calc = cartesian_dynamics.gradients_calculator_geometry_restraints(
          restraints_manager = self.restraints_manager)
    else:
      grad_calc = cartesian_dynamics.gradients_calculator_real_space_simple(
        restraints_manager        = self.restraints_manager.geometry, # XXX WHY?
        target_map                = self.target_map,
        unit_cell                 = self.xray_structure.unit_cell(),
        sites_cart                = self.xray_structure.sites_cart(),
        wx                        = self.wx,
        wc                        = self.wc,
        update_gradient_threshold = 0)
    return grad_calc
