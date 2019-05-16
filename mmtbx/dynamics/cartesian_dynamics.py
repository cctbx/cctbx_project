from __future__ import absolute_import, division, print_function
from mmtbx import dynamics
from mmtbx.dynamics.constants import \
  boltzmann_constant_akma, \
  akma_time_as_pico_seconds
from cctbx import xray
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
import random
import time
import math
import iotbx.phil
from cctbx import maptbx
from six.moves import range

def random_velocities(
      masses,
      target_temperature,
      zero_fraction=0,
      random_gauss=None,
      random_random=None,
      seed = None):
  result = flex.vec3_double()
  result.reserve(masses.size())
  if seed is not None:
    random.seed(seed)
  if (random_gauss is None): random_gauss = random.gauss
  if (random_random is None): random_random = random.random
  kt = boltzmann_constant_akma * target_temperature
  for mass in masses:
    assert mass > 0
    if (zero_fraction == 0 or random_random() >= zero_fraction):
      sigma = (kt / mass)**0.5
      result.append([random_gauss(0, sigma) for i in (1,2,3)])
    else:
      result.append([0,0,0])
  return result

class interleaved_lbfgs_minimization(object):

  def __init__(self,
        restraints_manager,
        sites_cart,
        max_iterations):
    self.restraints_manager = restraints_manager
    self.x = sites_cart.as_double()
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    sites_cart.clear()
    sites_cart.extend(flex.vec3_double(self.x))

  def compute_functional_and_gradients(self):
    sites_cart = flex.vec3_double(self.x)
    tmp = self.restraints_manager.energies_sites(sites_cart = sites_cart,
      compute_gradients=True)
    f = tmp.target
    g = tmp.gradients
    return f, g.as_double()

master_params = iotbx.phil.parse("""\
  temperature = 300
    .type = int
  number_of_steps = 200
    .type = int
  time_step = 0.0005
    .type = float
  initial_velocities_zero_fraction = 0
    .type = float
  n_print = 100
    .type = int
  verbose = -1
    .type = int
  random_seed = None
    .type = int
  n_collect = 10
    .type = int
  stop_cm_motion = True
    .type = bool
""")

class gradients_calculator_geometry_restraints(object):
  def __init__(self, restraints_manager = None):
    adopt_init_args(self, locals())
    self.gc = 0, 0

  def gradients(self, xray_structure, force_update_mask=False):
    factor = 1.0
    sites_cart = xray_structure.sites_cart()
    c = self.restraints_manager.energies_sites(sites_cart = sites_cart,
      compute_gradients=True)
    if(c.normalization_factor is not None): factor *= c.normalization_factor
    result = c.gradients
    #factor = 0.0001
    if(factor != 1.0): result *= 1.0 / factor
    #print flex.min(result.as_double()), flex.max(result.as_double()), \
    #  flex.mean(result.as_double()), result.norm()
    return result

class gradients_calculator_reciprocal_space(object):
  def __init__(self,
               restraints_manager        = None,
               fmodel                    = None,
               sites_cart                = None,
               wx                        = None,
               wc                        = None,
               update_gradient_threshold = 0):
    adopt_init_args(self, locals())
    assert [self.fmodel,             self.wx].count(None) in [0,2]
    assert [self.restraints_manager, self.wc].count(None) in [0,2]
    self.gx, self.gc = 0, 0
    if(self.fmodel is not None):
      self.x_target_functor = self.fmodel.target_functor()
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
      self.gx = flex.vec3_double(self.x_target_functor(compute_gradients=True).\
        gradients_wrt_atomic_parameters(site=True).packed())

  def gradients(self, xray_structure, force_update_mask=False):
    factor = 1.0
    sites_cart = xray_structure.sites_cart()
    if(self.fmodel is not None):
      max_shift = flex.max(flex.sqrt((self.sites_cart - sites_cart).dot()))
      if(max_shift > self.update_gradient_threshold):
        self.fmodel.update_xray_structure(
          xray_structure = xray_structure,
          update_f_calc  = True,
          update_f_mask  = False)
        self.gx = flex.vec3_double(self.x_target_functor(compute_gradients=True).\
          gradients_wrt_atomic_parameters(site=True).packed())
        self.sites_cart = sites_cart
    if(self.restraints_manager is not None):
      c = self.restraints_manager.energies_sites(sites_cart = sites_cart,
        compute_gradients=True)
      self.gc = c.gradients
      factor *= self.wc
      if(c.normalization_factor is not None): factor *= c.normalization_factor
    result = None
    if(self.wx is not None):
      result = self.wx * self.gx
    if(self.wc is not None):
      gcw = self.wc * self.gc
      if(result is None): result = gcw
      else: result = result + gcw
    if(factor != 1.0): result *= 1.0 / factor
    #print "norms:", self.gc.norm(), self.gx.norm(), result.norm()
    return result

class gradients_calculator_real_space_simple(object):
  def __init__(self,
               restraints_manager        = None,
               real_space_gradients_delta= 1./4,
               unit_cell                 = None,
               target_map                = None,
               sites_cart                = None,
               wx                        = None,
               wc                        = None,
               update_gradient_threshold = 0):
    adopt_init_args(self, locals())
    assert [self.target_map,         self.wx].count(None) in [0,2]
    assert [self.restraints_manager, self.wc].count(None) in [0,2]
    self.gx, self.gc = 0, 0
    if(self.target_map is not None):
      self.gx = self._compute_gradients()

  def _compute_gradients(self):
    return -1.*maptbx.real_space_gradients_simple(
      unit_cell   = self.unit_cell,
      density_map = self.target_map,
      sites_cart  = self.sites_cart,
      delta       = self.real_space_gradients_delta,
      selection   = flex.bool(self.sites_cart.size(), True))

  def gradients(self, xray_structure, force_update_mask=False):
    factor = 1.0
    sites_cart = xray_structure.sites_cart()
    if(self.target_map is not None):
      max_shift = flex.max(flex.sqrt((self.sites_cart - sites_cart).dot()))
      if(max_shift > self.update_gradient_threshold):
        self.sites_cart = sites_cart
        self.gx = self._compute_gradients()
    if(self.restraints_manager is not None):
      c = self.restraints_manager.energies_sites(sites_cart = sites_cart,
        compute_gradients=True)
      self.gc = c.gradients
      factor *= self.wc
      if(c.normalization_factor is not None): factor *= c.normalization_factor
    result = None
    if(self.wx is not None):
      result = self.wx * self.gx
    if(self.wc is not None):
      gcw = self.wc * self.gc
      if(result is None): result = gcw
      else: result = result + gcw
    if(factor != 1.0): result *= 1.0 / factor
    #print "norms:", self.gc.norm(), self.gx.norm(), result.norm()
    return result

class run(object):
  def __init__(self,
               xray_structure,
               gradients_calculator,
               temperature                      = 300,
               n_steps                          = 200,
               time_step                        = 0.0005,
               initial_velocities_zero_fraction = 0,
               vxyz                             = None,
               n_print                          = 20,
               n_collect                        = 10,
               interleaved_minimization         = False,
               reset_velocities                 = True,
               stop_cm_motion                   = False,
               log                              = None,
               stop_at_diff                     = None,
               random_seed                      = None,
               states_collector                 = None,
               verbose                          = -1):
    adopt_init_args(self, locals())
    assert self.n_print > 0
    assert self.temperature >= 0.0
    assert self.n_steps >= 0
    assert self.time_step >= 0.0
    assert self.log is not None or self.verbose < 1
    self.sites_cart_start = self.xray_structure.sites_cart()
    if(self.states_collector is not None):
      self.states_collector.add(sites_cart = self.sites_cart_start)
    self.k_boltz = boltzmann_constant_akma
    self.current_temperature = 0.0
    self.ekin = 0.0
    self.ekcm = 0.0
    self.timfac = akma_time_as_pico_seconds
    self.atomic_weights = self.xray_structure.atomic_weights()
    if(vxyz is None):
      self.vxyz = flex.vec3_double(self.atomic_weights.size(),(0,0,0))
    else:
      self.vxyz = vxyz
    #
    self.tstep = self.time_step / self.timfac
    self()

  def __call__(self):
    self.center_of_mass_info()
    #print "0:",self.temperature, self.current_temperature
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
    self.current_temperature = kt.temperature
    #print "1:",self.temperature, self.current_temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="restrained dynamics start")
    if(self.reset_velocities):
       self.set_velocities()
       self.center_of_mass_info()
       kt=dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
       self.current_temperature = kt.temperature
       self.ekin = kt.kinetic_energy
       if(self.verbose >= 1):
         self.print_dynamics_stat(text="set velocities")
    if(self.stop_cm_motion):
      self.stop_global_motion()
    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
    self.current_temperature = kt.temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="center of mass motion removed")
    self.velocity_rescaling()
    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
    self.current_temperature = kt.temperature
    #print "2:",self.temperature, self.current_temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="velocities rescaled")
    if(self.verbose >= 1):
      print("integration starts", file=self.log)
    self.verlet_leapfrog_integration()
    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
    self.current_temperature = kt.temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="after final integration step")

  def run_interleaved_minimization(self):
    geo_manager = self.gradients_calculator.restraints_manager.geometry
    sites_cart = self.xray_structure.sites_cart()
    interleaved_lbfgs_minimization(
      restraints_manager = self.gradients_calculator.restraints_manager,
      sites_cart=sites_cart,
      max_iterations=5)
    self.xray_structure.set_sites_cart(sites_cart=sites_cart)
    self.xray_structure.apply_symmetry_sites()

  def set_velocities(self):
    self.vxyz.clear()
    self.vxyz.extend(random_velocities(
      masses=self.atomic_weights,
      target_temperature=self.temperature,
      zero_fraction=self.initial_velocities_zero_fraction,
      seed = self.random_seed))

  def accelerations(self):
    return self.gradients_calculator.gradients(xray_structure=self.xray_structure)

  def center_of_mass_info(self):
    self.rcm = self.xray_structure.center_of_mass()
    result = dynamics.center_of_mass_info(
      self.rcm,
      self.xray_structure.sites_cart(),
      self.vxyz,
      self.atomic_weights)
    self.vcm = flex.vec3_double()
    self.acm = flex.vec3_double()
    self.vcm.append(result.vcm())
    self.acm.append(result.acm())
    self.ekcm = result.ekcm()

  def stop_global_motion(self):
    self.rcm = self.xray_structure.center_of_mass()
    self.vxyz = dynamics.stop_center_of_mass_motion(
      self.rcm,
      self.acm[0],
      self.vcm[0],
      self.xray_structure.sites_cart(),
      self.vxyz,
      self.atomic_weights)

  def velocity_rescaling(self):
    if(self.current_temperature <= 1.e-10):
      factor = 1.0
    else:
      factor = math.sqrt(self.temperature/self.current_temperature)
    self.vxyz = self.vxyz * factor

  def verlet_leapfrog_integration(self):
    # start verlet_leapfrog_integration loop
    for cycle in range(1,self.n_steps+1,1):
      sites_cart = None
      if([self.stop_at_diff,self.states_collector].count(None) != 2):
        sites_cart = self.xray_structure.sites_cart()
      if(self.stop_at_diff is not None):
        dist = flex.mean(flex.sqrt((self.sites_cart_start - sites_cart).dot()))
        if(dist >= self.stop_at_diff): return
      accelerations = self.accelerations()
      print_flag = 0
      switch = math.modf(float(cycle)/self.n_print)[0]
      if((switch==0 or cycle==1 or cycle==self.n_steps) and self.verbose >= 1):
        print_flag = 1
      if(self.states_collector is not None):
        switch2 = math.modf(float(cycle)/self.n_collect)[0]
        if(switch2==0 or cycle==1 or cycle==self.n_steps):
          self.states_collector.add(sites_cart = sites_cart)
      if(print_flag == 1):
        text = "integration step number = %5d"%cycle
        self.center_of_mass_info()
        kt=dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
        self.current_temperature = kt.temperature
        self.ekin = kt.kinetic_energy
        self.print_dynamics_stat(text)
      if(self.stop_cm_motion):
        self.center_of_mass_info()
        self.stop_global_motion()
      # calculate velocities at t+dt/2
      dynamics.vxyz_at_t_plus_dt_over_2(
        self.vxyz, self.atomic_weights, accelerations, self.tstep)
      # calculate the temperature and kinetic energy from new velocities
      kt=dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
      self.current_temperature = kt.temperature
      self.ekin = kt.kinetic_energy
      self.velocity_rescaling()
      if(print_flag == 1 and 0):
        self.center_of_mass_info()
        self.print_dynamics_stat(text)
      # do the verlet_leapfrog_integration to get coordinates at t+dt
      self.xray_structure.set_sites_cart(
        sites_cart=self.xray_structure.sites_cart() + self.vxyz * self.tstep)
      self.xray_structure.apply_symmetry_sites()
      # prevent explosions by doing very quick model geometry regularization
      if(self.interleaved_minimization and cycle==self.n_steps):
        self.run_interleaved_minimization()
      kt=dynamics.kinetic_energy_and_temperature(self.vxyz,self.atomic_weights)
      self.current_temperature = kt.temperature
      self.ekin = kt.kinetic_energy
      if(print_flag == 1 and 0):
        self.center_of_mass_info()
        self.print_dynamics_stat(text)
      self.accelerations()

  def print_dynamics_stat(self, text):
    timfac = akma_time_as_pico_seconds
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    print("| "+text+"-"*(fill_len)+"|", file=self.log)
    print("| kin.energy = %10.3f            " \
      "| information about center of free masses|"%(self.ekin), file=self.log)
    print("| start temperature = %7.3f        " \
      "| position=%8.3f%8.3f%8.3f      |"% (
      self.temperature,self.rcm[0],self.rcm[1],self.rcm[2]), file=self.log)
    print("| curr. temperature = %7.3f        " \
      "| velocity=%8.4f%8.4f%8.4f      |"% (self.current_temperature,
      self.vcm[0][0]/timfac,self.vcm[0][1]/timfac,self.vcm[0][2]/timfac), file=self.log)
    print("| number of integration steps = %4d " \
      "| ang.mom.=%10.2f%10.2f%10.2f|"% (self.n_steps,
      self.acm[0][0]/timfac,self.acm[0][1]/timfac,self.acm[0][2]/timfac), file=self.log)
    print("| time step = %6.4f                 | kin.ener.=%8.3f                     |"% (
      self.time_step,self.ekcm), file=self.log)
    print("|"+"-"*77+"|", file=self.log)
