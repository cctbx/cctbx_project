from mmtbx import dynamics
from mmtbx.dynamics.constants import \
  boltzmann_constant_akma, \
  akma_time_as_pico_seconds
from cctbx import geometry_restraints
from cctbx import xray
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import random
import time
import math
import iotbx.phil

def random_velocities(
      masses,
      target_temperature,
      zero_fraction=0,
      random_gauss=None,
      random_random=None):
  result = flex.vec3_double()
  result.reserve(masses.size())
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
        conservative_pair_proxies,
        sites_cart,
        max_iterations):
    self.conservative_pair_proxies = conservative_pair_proxies
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
    f = 0
    g = flex.vec3_double(sites_cart.size(), (0,0,0))
    for sorted_asu_proxies in [self.conservative_pair_proxies.bond,
                               self.conservative_pair_proxies.angle]:
      if (sorted_asu_proxies is None): continue
      f += geometry_restraints.bond_residual_sum(
        sites_cart=sites_cart,
        sorted_asu_proxies=sorted_asu_proxies,
        gradient_array=g)
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
""")

class cartesian_dynamics(object):
  def __init__(self,
               structure,
               restraints_manager,
               temperature                      = 300,
               n_steps                          = 200,
               time_step                        = 0.0005,
               initial_velocities_zero_fraction = 0,
               vxyz                             = None,
               interleaved_minimization_params  = None,
               n_print                          = 20,
               fmodel                           = None,
               xray_target_weight               = None,
               chem_target_weight               = None,
               shift_update                     = 0.0,
               xray_structure_last_updated      = None,
               xray_gradient                    = None,
               reset_velocities                 = True,
               stop_cm_motion                   = False,
               update_f_calc                    = True,
               time_averaging_data              = None,
               log                              = None,
               stop_at_diff                     = None,
               verbose                          = -1):
    adopt_init_args(self, locals())
    assert self.n_print > 0
    assert self.temperature >= 0.0
    assert self.n_steps >= 0
    assert self.time_step >= 0.0
    assert self.log is not None or self.verbose < 1
    xray.set_scatterer_grad_flags(scatterers = self.structure.scatterers(),
                                  site       = True)
    self.structure_start = self.structure.deep_copy_scatterers()
    self.k_boltz = boltzmann_constant_akma
    self.current_temperature = 0.0
    self.ekin = 0.0
    self.ekcm = 0.0
    self.timfac = akma_time_as_pico_seconds
    self.weights = self.structure.atomic_weights()
    if(vxyz is None):
      self.vxyz = flex.vec3_double(self.weights.size(),(0,0,0))
    else:
      self.vxyz = vxyz
    if(self.time_averaging_data is not None and
       self.time_averaging_data.velocities is not None):
      self.vxyz = self.time_averaging_data.velocities
    if(self.fmodel is not None):
      self.fmodel_copy = self.fmodel.deep_copy()
      self.target_functor = self.fmodel_copy.target_functor()
      assert self.chem_target_weight is not None
      assert self.xray_target_weight is not None
      if(self.xray_gradient is None):
        self.xray_gradient = self.xray_grads()
    #
    imp = self.interleaved_minimization_params
    self.interleaved_minimization_flag = (
      imp is not None and imp.number_of_iterations > 0)
    if (self.interleaved_minimization_flag):
      assert imp.time_step_factor > 0
      self.time_step *= imp.time_step_factor
      if ("bonds" not in imp.restraints):
        raise Sorry(
          'Invalid choice: %s: "bonds" must always be included.'
            % imp.__phil_path__(object_name="restraints"))
      self.interleaved_minimization_angles = "angles" in imp.restraints
    else:
      self.interleaved_minimization_angles = False
    #
    self.tstep = self.time_step / self.timfac
    self.show_gradient_rms = False # XXX debug option
    #
    self()

  def __call__(self):
    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)
    self.current_temperature = kt.temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="restrained dynamics start")
    if(self.reset_velocities):
       self.set_velocities()
       self.center_of_mass_info()
       kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)
       self.current_temperature = kt.temperature
       self.ekin = kt.kinetic_energy
       if(self.verbose >= 1):
         self.print_dynamics_stat(text="set velocities")

    if(self.stop_cm_motion):
      self.stop_global_motion()
    self.center_of_mass_info()

    if(self.time_averaging_data is None):
      kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)
      self.current_temperature = kt.temperature
      self.ekin = kt.kinetic_energy
      if(self.verbose >= 1):
        self.print_dynamics_stat(text="center of mass motion removed")
      self.velocity_rescaling()

    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)
    self.current_temperature = kt.temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="velocities rescaled")

    if(self.verbose >= 1):
      print >> self.log, "integration starts"

    self.verlet_leapfrog_integration()

    self.center_of_mass_info()
    kt = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)
    self.current_temperature = kt.temperature
    self.ekin = kt.kinetic_energy
    if(self.verbose >= 1):
      self.print_dynamics_stat(text="after final integration step")

  def set_velocities(self):
    self.vxyz.clear()
    self.vxyz.extend(random_velocities(
      masses=self.weights,
      target_temperature=self.temperature,
      zero_fraction=self.initial_velocities_zero_fraction))

  def accelerations(self):
    stereochemistry_residuals = self.restraints_manager.energies_sites(
      sites_cart=self.structure.sites_cart(),
      compute_gradients=True)
    result = stereochemistry_residuals.gradients
    d_max = None
    if(self.xray_structure_last_updated is not None and self.shift_update > 0):
      array_of_distances_between_each_atom = \
        flex.sqrt(self.structure.difference_vectors_cart(
           self.xray_structure_last_updated).dot())
      d_max = flex.max(array_of_distances_between_each_atom)
    if(self.fmodel is not None):
      if(d_max is not None):
        if(d_max > self.shift_update):
          self.xray_structure_last_updated = self.structure.deep_copy_scatterers()
          self.xray_gradient = self.xray_grads()
      else:
        self.xray_gradient = self.xray_grads()
      result = self.xray_gradient * self.xray_target_weight \
             + stereochemistry_residuals.gradients * self.chem_target_weight
    factor = 1.0
    if (self.chem_target_weight is not None):
      factor *= self.chem_target_weight
    if (stereochemistry_residuals.normalization_factor is not None):
      factor *= stereochemistry_residuals.normalization_factor
    if (factor != 1.0):
      result *= 1.0 / factor
    if (self.show_gradient_rms and self.xray_gradient is not None):
      gg = stereochemistry_residuals.gradients*(self.chem_target_weight/factor)
      xg = self.xray_gradient*(self.xray_target_weight/factor)
      print "CARDY_GRADIENT_RMS_GEO_XRAY:", \
        flex.mean_sq(gg.as_double())**0.5, \
        flex.mean_sq(xg.as_double())**0.5
    return result

  def xray_grads(self):
    self.fmodel_copy.update_xray_structure(
      xray_structure = self.structure,
      update_f_calc  = True,
      update_f_mask  = False)
    if(self.time_averaging_data is not None):
      self.time_averaging_data.xray_structures.append(
        self.structure.deep_copy_scatterers())
      fcra = self.time_averaging_data.f_calc_running_average
      fmda = self.time_averaging_data.f_model_data_average
      self.time_averaging_data.f_model_data_current = \
        self.fmodel_copy.f_model().data()
      if(fmda is None):
        fmda = self.time_averaging_data.f_model_data_current
      else:
        fmda += self.time_averaging_data.f_model_data_current
      self.time_averaging_data.f_model_data_average = fmda
    if(not self.update_f_calc):
      f_calc = self.fmodel_copy.f_calc()
      if(fcra is None):
        fcra = f_calc
      else:
        a_prime = math.exp(
          -self.time_step/self.time_averaging_data.params.tx)
        f_calc_average_data = fcra.data()
        f_calc_average_data = a_prime * f_calc_average_data + \
                         (1.-a_prime) * f_calc.data()
        fcra = fcra.array(data = f_calc_average_data)
      self.fmodel_copy.update(f_calc = fcra)
      self.time_averaging_data.f_calc_running_average = fcra
      self.time_averaging_data.r_work_running_average = \
        self.fmodel_copy.r_work()
      self.time_averaging_data.r_free_running_average = \
        self.fmodel_copy.r_free()
    sf = self.target_functor(
      compute_gradients=True).gradients_wrt_atomic_parameters(site=True)
    if(self.time_averaging_data is not None):
      self.time_averaging_data.velocities = self.vxyz
    return flex.vec3_double(sf.packed())

  def center_of_mass_info(self):
    self.rcm = self.structure.center_of_mass()
    result = dynamics.center_of_mass_info(
      self.rcm,
      self.structure.sites_cart(),
      self.vxyz,
      self.weights)
    self.vcm = flex.vec3_double()
    self.acm = flex.vec3_double()
    self.vcm.append(result.vcm())
    self.acm.append(result.acm())
    self.ekcm = result.ekcm()

  def stop_global_motion(self):
    self.rcm = self.structure.center_of_mass()
    self.vxyz = dynamics.stop_center_of_mass_motion(
      self.rcm,
      self.acm[0],
      self.vcm[0],
      self.structure.sites_cart(),
      self.vxyz,
      self.weights)

  def velocity_rescaling(self):
    if(self.time_averaging_data is not None):
      self.time_averaging_data.temp.append(self.current_temperature)
    if (self.current_temperature <= 1.e-10):
      factor = 1.0
    else:
      factor = math.sqrt(self.temperature/self.current_temperature)
    self.vxyz = self.vxyz * factor

  def interleaved_minimization(self):
    geo_manager = self.restraints_manager.geometry
    assert geo_manager.shell_sym_tables is not None
    assert len(geo_manager.shell_sym_tables) > 0
    conservative_pair_proxies = self.structure.conservative_pair_proxies(
      bond_sym_table=geo_manager.shell_sym_tables[0],
      conserve_angles=self.interleaved_minimization_angles)
    sites_cart = self.structure.sites_cart()
    interleaved_lbfgs_minimization(
      sites_cart=sites_cart,
      conservative_pair_proxies=conservative_pair_proxies,
      max_iterations=self.interleaved_minimization_params.number_of_iterations)
    self.structure.set_sites_cart(sites_cart=sites_cart)
    self.structure.apply_symmetry_sites()

  def verlet_leapfrog_integration(self):
    # start verlet_leapfrog_integration loop
    for cycle in range(1,self.n_steps+1,1):
      if(self.stop_at_diff is not None):
        diff = flex.mean(self.structure_start.distances(other = self.structure))
        if(diff >= self.stop_at_diff): return
      accelerations = self.accelerations()
      print_flag = 0
      switch = math.modf(float(cycle)/self.n_print)[0]
      if((switch==0 or cycle==1 or cycle==self.n_steps) and self.verbose >= 1):
        print_flag = 1
      if(print_flag == 1):
        text = "integration step number = %5d"%cycle
        self.center_of_mass_info()
        kt = dynamics.kinetic_energy_and_temperature(self.vxyz, self.weights)
        self.current_temperature = kt.temperature
        self.ekin = kt.kinetic_energy
        self.print_dynamics_stat(text)
      if(self.stop_cm_motion):
        self.center_of_mass_info()
        self.stop_global_motion()
      # calculate velocities at t+dt/2
      dynamics.vxyz_at_t_plus_dt_over_2(
        self.vxyz, self.weights, accelerations, self.tstep)
      # calculate the temperature and kinetic energy from new velocities
      kt = dynamics.kinetic_energy_and_temperature(self.vxyz, self.weights)
      self.current_temperature = kt.temperature
      self.ekin = kt.kinetic_energy
      self.velocity_rescaling()
      if(print_flag == 1 and 0):
        self.center_of_mass_info()
        self.print_dynamics_stat(text)
      # do the verlet_leapfrog_integration to get coordinates at t+dt
      self.structure.set_sites_cart(
        sites_cart=self.structure.sites_cart() + self.vxyz * self.tstep)
      self.structure.apply_symmetry_sites()
      if (self.interleaved_minimization_flag):
        self.interleaved_minimization()
      kt = dynamics.kinetic_energy_and_temperature(self.vxyz, self.weights)
      self.current_temperature = kt.temperature
      self.ekin = kt.kinetic_energy
      if(print_flag == 1 and 0):
        self.center_of_mass_info()
        self.print_dynamics_stat(text)
      if(self.time_averaging_data is None):
        self.accelerations()
      else:
        self.time_averaging_data.velocities = self.vxyz

  def print_dynamics_stat(self, text):
    timfac = akma_time_as_pico_seconds
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    print >> self.log, "| "+text+"-"*(fill_len)+"|"
    print >> self.log, "| kin.energy = %10.3f            " \
      "| information about center of free masses|"%(self.ekin)
    print >> self.log, "| start temperature = %7.3f        " \
      "| position=%8.3f%8.3f%8.3f      |"% (
      self.temperature,self.rcm[0],self.rcm[1],self.rcm[2])
    print >> self.log, "| curr. temperature = %7.3f        " \
      "| velocity=%8.4f%8.4f%8.4f      |"% (self.current_temperature,
      self.vcm[0][0]/timfac,self.vcm[0][1]/timfac,self.vcm[0][2]/timfac)
    print >> self.log, "| number of integration steps = %4d " \
      "| ang.mom.=%10.2f%10.2f%10.2f|"% (self.n_steps,
      self.acm[0][0]/timfac,self.acm[0][1]/timfac,self.acm[0][2]/timfac)
    print >> self.log, "| time step = %6.4f                 | kin.ener.=%8.3f                     |"% (
      self.time_step,self.ekcm)
    print >> self.log, "|"+"-"*77+"|"
