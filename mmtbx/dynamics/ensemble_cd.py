from __future__ import absolute_import, division, print_function
from mmtbx import dynamics
from mmtbx.dynamics.constants import \
  boltzmann_constant_akma, \
  akma_time_as_pico_seconds
from cctbx import geometry_restraints
from cctbx import xray
from cctbx.array_family import flex
import scitbx.lbfgs
import scitbx.math
from libtbx.utils import Sorry
from libtbx import adopt_init_args
import random
import math
import sys
import iotbx.phil
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
               temperature                        = 300,
               protein_thermostat                 = False,
               n_steps                            = 200,
               time_step                          = 0.0005,
               initial_velocities_zero_fraction   = 0,
               vxyz                               = None,
               interleaved_minimization_params    = None,
               n_print                            = 20,
               fmodel                             = None,
               xray_target_weight                 = None,
               chem_target_weight                 = None,
               shift_update                       = 0.0,
               xray_structure_last_updated        = None,
               xray_gradient                      = None,
               reset_velocities                   = True,
               stop_cm_motion                     = False,
               update_f_calc                      = True,
               er_data                            = None,
               log                                = None,
               stop_at_diff                       = None,
               verbose                            = -1):
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
    self.ekcm = 0.0
    self.timfac = akma_time_as_pico_seconds
    self.weights = self.structure.atomic_weights()
    if(vxyz is None):
      self.vxyz = flex.vec3_double(self.weights.size(),(0,0,0))
    else:
      self.vxyz = vxyz
    if(self.er_data is not None and
       self.er_data.velocities is not None):
      self.vxyz = self.er_data.velocities

    if self.er_data is not None:
      self.er_data.geo_grad_rms = 0
      self.er_data.xray_grad_rms = 0

    if(self.fmodel is not None):
      if self.er_data is None:
        self.fmodel_copy = self.fmodel.deep_copy()
      else:
        self.fmodel_copy = self.fmodel
        if self.er_data.fix_scale_factor is not None:
          self.fmodel_copy.set_scale_switch = self.er_data.fix_scale_factor
    #
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
          'Invalid choice: %s.restraints: "bonds" must always be included.'
            % imp.__phil_path__())
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
    if self.reset_velocities:
       self.set_velocities()
       self.center_of_mass_info()

    self.non_solvent_vxyz    = self.vxyz.select(~self.er_data.solvent_sel)
    self.non_solvent_weights = self.weights.select(~self.er_data.solvent_sel)
    self.solvent_vxyz        = self.vxyz.select(self.er_data.solvent_sel)
    self.solvent_weights     = self.weights.select(self.er_data.solvent_sel)
        #
    self.non_solvent_kt      = dynamics.kinetic_energy_and_temperature(self.non_solvent_vxyz, self.non_solvent_weights)
    self.solvent_kt          = dynamics.kinetic_energy_and_temperature(self.solvent_vxyz, self.solvent_weights)
    self.kt                  = dynamics.kinetic_energy_and_temperature(self.vxyz,self.weights)

    if self.stop_cm_motion:
      self.stop_global_motion()

    self.velocity_rescaling()

    self.verlet_leapfrog_integration()
    if self.verbose >= 1.0:
      self.print_detailed_dynamics_stats()

  def set_velocities(self):
    self.vxyz.clear()
    if self.er_data is not None:
      seed = self.er_data.seed
    else: seed = None
    self.vxyz.extend(random_velocities(
      masses=self.weights,
      target_temperature=self.temperature,
      zero_fraction=self.initial_velocities_zero_fraction,
      seed = seed))

  def accelerations(self):
    self.stereochemistry_residuals = self.restraints_manager.energies_sites(
      sites_cart=self.structure.sites_cart(),
      compute_gradients=True)

    # Harmonic restraints
    if self.er_data is not None:
      if self.er_data.er_harmonic_restraints_info is not None:
        harmonic_grads = self.restraints_manager.geometry.ta_harmonic_restraints(
                                sites_cart = self.structure.sites_cart(),
                                ta_harmonic_restraint_info = self.er_data.er_harmonic_restraints_info,
                                weight = self.er_data.er_harmonic_restraints_weight,
                                slack = self.er_data.er_harmonic_restraints_slack)
        assert self.stereochemistry_residuals.gradients.size() == harmonic_grads.size()
        self.stereochemistry_residuals.gradients += harmonic_grads
    result = self.stereochemistry_residuals.gradients

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
             + self.stereochemistry_residuals.gradients * self.chem_target_weight

    factor = 1.0
    if (self.chem_target_weight is not None):
      factor *= self.chem_target_weight
    if (self.stereochemistry_residuals.normalization_factor is not None):
      factor *= self.stereochemistry_residuals.normalization_factor


    if (factor != 1.0):
      result *= 1.0 / factor

    #Store RMS non-solvent atom gradients for Xray and Geo
    if self.er_data is not None:
      self.wc = self.chem_target_weight / factor
      self.wx = self.xray_target_weight / factor
      self.gg = self.stereochemistry_residuals.gradients * self.wc
      self.xg = self.xray_gradient * self.wx
      gg_pro = self.gg.select( ~self.er_data.solvent_sel )
      xg_pro = self.xg.select( ~self.er_data.solvent_sel )
      self.er_data.geo_grad_rms  += (flex.mean_sq(gg_pro.as_double())**0.5) / self.n_steps
      self.er_data.xray_grad_rms += (flex.mean_sq(xg_pro.as_double())**0.5) / self.n_steps

    return result

  def xray_grads(self):
    self.fmodel_copy.update_xray_structure(
      xray_structure           = self.structure,
      update_f_calc            = self.update_f_calc,
      update_f_mask            = False)
    sf = self.target_functor(
        compute_gradients=True).gradients_wrt_atomic_parameters(site=True)
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
    if self.protein_thermostat and self.er_data is not None:
      if (self.kt.temperature <= 1.e-10):
        self.v_factor = 1.0
      else:
        self.v_factor = math.sqrt(self.temperature/self.non_solvent_kt.temperature)
    else:
      if (self.kt.temperature <= 1.e-10):
        self.v_factor = 1.0
      else:
        self.v_factor = math.sqrt(self.temperature/self.kt.temperature)

    self.vyz_vscale_remove = self.vxyz * (1.0 - self.v_factor)
    self.kt_vscale_remove = dynamics.kinetic_energy_and_temperature(self.vyz_vscale_remove, self.weights)
    self.vxyz = self.vxyz * self.v_factor

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
    for self.cycle in range(1,self.n_steps+1,1):
      if(self.stop_at_diff is not None):
        diff = flex.mean(self.structure_start.distances(other = self.structure))
        if(diff >= self.stop_at_diff): return
      accelerations = self.accelerations()
      if(self.stop_cm_motion):
        self.center_of_mass_info()
        self.stop_global_motion()

      # calculate velocities at t+dt/2
      dynamics.vxyz_at_t_plus_dt_over_2(
        self.vxyz, self.weights, accelerations, self.tstep)

      # calculate the temperature and kinetic energy from new velocities
      self.non_solvent_vxyz    = self.vxyz.select(~self.er_data.solvent_sel)
      self.non_solvent_weights = self.weights.select(~self.er_data.solvent_sel)
      self.solvent_vxyz        = self.vxyz.select(self.er_data.solvent_sel)
      self.solvent_weights     = self.weights.select(self.er_data.solvent_sel)
      #
      self.kt                  = dynamics.kinetic_energy_and_temperature(self.vxyz, self.weights)
      self.non_solvent_kt      = dynamics.kinetic_energy_and_temperature(self.non_solvent_vxyz, self.non_solvent_weights)
      self.solvent_kt          = dynamics.kinetic_energy_and_temperature(self.solvent_vxyz, self.solvent_weights)
      #Store sys, solvent, nonsolvet temperatures
      if self.er_data is not None:
        self.store_temperatures()
      self.velocity_rescaling()
      if self.verbose >= 1.0:
        print('Scale factor : ', self.v_factor, file=self.log)
        self.vxyz_length_sq = flex.sqrt(self.vxyz.dot())
        print('vxyz_length_sq pst scale', file=self.log)
        self.vxyz_length_sq.min_max_mean().show(out=self.log)
      # do the verlet_leapfrog_integration to get coordinates at t+dt
      self.structure.set_sites_cart(
        sites_cart=self.structure.sites_cart() + self.vxyz * self.tstep)
      self.structure.apply_symmetry_sites()
      if (self.interleaved_minimization_flag):
        self.interleaved_minimization()
      self.kt = dynamics.kinetic_energy_and_temperature(self.vxyz, self.weights)
      if(self.er_data is None):
        self.accelerations()
      else:
        self.er_data.velocities = self.vxyz

  def store_temperatures(self):
    if self.cycle == 1:
      self.er_data.non_solvent_temp = 0
      self.er_data.solvent_temp = 0
      self.er_data.system_temp  = 0
    self.er_data.non_solvent_temp += (self.non_solvent_kt.temperature / self.n_steps)
    self.er_data.solvent_temp += (self.solvent_kt.temperature /self.n_steps)
    self.er_data.system_temp  += (self.kt.temperature /self.n_steps)

    ### Extra dynamics stats
  def print_detailed_dynamics_stats(self):
    # Overall data
    print('\n', file=self.log)
    print('         MC |    Temperature (K)   |   Vscale   | Etot = Ekin + Echem + wxExray', file=self.log)
    print('            |  (sys)  (pro)  (sol) |  Fac  T(K) |   Ekin  Echem     wx  Exray', file=self.log)
    print('  ~E~ {0:5d} | {1:6.1f} {2:6.1f} {3:6.1f} | {4:4.1f} {5:5.1f} | {6:6.1f} {7:6.1f} {8:6.1f} {9:6.1f}'.format(
        self.er_data.macro_cycle,
        self.kt.temperature,
        self.non_solvent_kt.temperature,
        self.solvent_kt.temperature,
        self.v_factor,
        self.kt_vscale_remove.temperature,
        self.kt.kinetic_energy,
        self.stereochemistry_residuals.residual_sum,
        self.xray_target_weight,
        self.target_functor(compute_gradients=False).target_work() * self.fmodel_copy.f_calc_w().data().size(),
        ), file=self.log)
    print('\n', file=self.log)

    # Atomistic histrograms
    # - Kinetic energy
    # - Xray grads
    # - Geo grads
    self.atomic_ke = 0.5 * self.weights * self.vxyz.dot()
    self.atomic_wxray_g = self.xray_gradient * self.xray_target_weight
    self.atomic_wchem_g = self.stereochemistry_residuals.gradients * self.chem_target_weight

    def show_histogram(data,
                       n_slots = 50,
                       out     = None,
                       prefix  = ""):
      if (out is None): out = sys.stdout
      print('\n' + prefix, file=out)

      # Stats
      data_basic_stats = scitbx.math.basic_statistics(data)
      print('\n  Number  : %7.4f ' % (data_basic_stats.n), file=out)
      print('  Min     : %7.4f ' % (data_basic_stats.min), file=out)
      print('  Max     : %7.4f ' % (data_basic_stats.max), file=out)
      print('  Mean    : %7.4f ' % (data_basic_stats.mean), file=out)
      print('  Stdev   : %7.4f ' % (data_basic_stats.biased_standard_deviation), file=out)
      print('  Skew    : %7.4f ' % (data_basic_stats.skew), file=out)
      print('  Sum     : %7.4f ' % (data_basic_stats.sum), file=out)

      # Histo
      histogram = flex.histogram(data    = data,
                                 n_slots = n_slots)
      low_cutoff = histogram.data_min()
      for i,n in enumerate(histogram.slots()):
        high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
        print("%7.3f - %7.3f: %d" % (low_cutoff, high_cutoff, n), file=out)
        low_cutoff = high_cutoff
      out.flush()
      return histogram

    # Select
    for selection_type in ['System', 'Non_solvent', 'Solvent']:
      print('\n\n', file=self.log)
      if selection_type == 'System':
        selection = self.er_data.all_sel
      elif selection_type == 'Non_solvent':
        selection = ~self.er_data.solvent_sel
      elif selection_type == 'Solvent':
        selection = self.er_data.solvent_sel
      else:
        break
      # Data
      for histogram_type in ['Kinetic_energy', 'Xray_grad', 'Chem_grad']:
        if histogram_type == 'Kinetic_energy':
          data = self.atomic_ke.select(selection)
        elif histogram_type == 'Xray_grad':
          data = flex.sqrt(self.atomic_wxray_g.select(selection).dot())
        elif histogram_type == 'Chem_grad':
          data = flex.sqrt(self.atomic_wchem_g.select(selection).dot())
        else:
          break
        # Histrogram
        show_histogram(data    = data,
                       out     = self.log,
                       prefix  = str(self.er_data.macro_cycle) + '_' + selection_type + '_' + histogram_type)
