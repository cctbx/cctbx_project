
from __future__ import division
import mmtbx.solvent.ensemble_ordered_solvent as ensemble_ordered_solvent
from mmtbx.refinement.ensemble_refinement import ensemble_utils
from mmtbx.dynamics import ensemble_cd
import mmtbx.tls.tools as tls_tools
import mmtbx.command_line
import mmtbx.utils
import mmtbx.model
import mmtbx.maps
from iotbx.option_parser import iotbx_option_parser
from iotbx import pdb
import iotbx.phil
import iotbx
from cctbx import geometry_restraints
from cctbx.array_family import flex
from cctbx import miller
from cctbx import adptbx
from cctbx import xray
import scitbx.math
from libtbx.utils import Sorry, user_plus_sys_time, multi_out, show_total_time
from libtbx import adopt_init_args, slots_getstate_setstate
from libtbx.str_utils import format_value, make_header
from libtbx import runtime_utils
from libtbx import easy_mp
import libtbx.load_env
from cStringIO import StringIO
import cPickle
import random
import gzip
import math
import time
import os
import sys

# these supersede the defaults in included scopes
customization_params = iotbx.phil.parse("""
ensemble_refinement.mask.ignore_hydrogens = False
ensemble_refinement.mask.n_radial_shells = 1
ensemble_refinement.mask.radial_shell_width = 1.5
ensemble_refinement.cartesian_dynamics.number_of_steps = 10
ensemble_refinement.ensemble_ordered_solvent.b_iso_min = 0.0
ensemble_refinement.ensemble_ordered_solvent.b_iso_max = 100.0
ensemble_refinement.ensemble_ordered_solvent.find_peaks.map_next_to_model.max_model_peak_dist = 3.0
ensemble_refinement.ensemble_ordered_solvent.find_peaks.map_next_to_model.use_hydrogens = False
""")

# the extra fetch() at the end with the customized parameters gives us the
# actual master phil object
input_phil = mmtbx.command_line.generate_master_phil_with_inputs(
  phil_string="",
  enable_experimental_phases=True,
  as_phil_string=True)
master_params = iotbx.phil.parse(input_phil + """
ensemble_refinement {
  cartesian_dynamics
    .style = menu_item auto_align box
  {
    include scope mmtbx.dynamics.cartesian_dynamics.master_params
    protein_thermostat = True
      .type = bool
      .help = Use protein atoms thermostat
  }
  update_sigmaa_rfree = 0.001
    .type = float
    .help = test function
  ensemble_reduction = True
    .type = bool
    .help = 'Find miminium number of structures to reproduce simulation R-values'
  verbose = -1
    .type = int
  output_file_prefix = None
    .type = str
    .help = 'Prefix for all output files'
# TODO
#  write_mmcif_file = False
#    .type = bool
  gzip_final_model = True
    .type = bool
    .style = hidden
  random_seed = 2679941
    .type = int
    .help = 'Random seed'
  nproc = 1
    .type = int
    .short_caption = Number of processors
    .style = renderer:draw_nproc_widget
  tx = None
    .type = float
    .help = 'Relaxation time (ps)'
    .short_caption = Relaxation time (ps)
  equilibrium_n_tx = 10
    .type = int
    .help = 'Length of equilibration period, n times tx'
    .short_caption = Length of equilibration period
  acquisition_block_n_tx = 2
    .type = int
    .help = 'Length of acquisition block, n times tx'
    .short_caption = Length of acquisition block
  number_of_acquisition_periods = 10
    .type = int
    .help = 'Number of acquisition periods'
  pdb_stored_per_block = 50
    .type = int
    .help = 'Number of model coordinates stored per acquisition block'
    .short_caption = Models stored per acquisition block
  wxray_coupled_tbath = True
    .type = bool
    .help = 'Use temperature control wxray'
    .short_caption = Use temperature control X-ray weight
  wxray_coupled_tbath_offset = 5.0
    .type = float
    .help = 'Temperature offset, increasing offset increases wxray'
    .short_caption = Temperature ofset
  wxray = 1.0
    .type = float
    .help = 'Multiplier for xray weighting; used if wxray_coupled_tbath = Flase'
    .short_caption = X-ray weight
  tls_group_selections = None
    .type = atom_selection
    .multiple = True
    .help = 'TLS groups to use for TLS fitting (TLS details in PDB header not used)'
    .style = use_list
  ptls = 0.80
    .type = floats
    .optional = False
    .help = 'The fraction of atoms to include in TLS fitting'
    .short_caption = Fraction of atoms to include in TLS fitting
    .style = bold
  max_ptls_cycles = 25
    .type = int
    .help = 'Maximum cycles to use in TLS fitting; TLS will stop prior to this if convergence is reached'
    .short_caption = Max. number of cycles of TLS fitting
  isotropic_b_factor_model = False
    .type = bool
    .help = 'Use isotropic B-factor model instead of TLS'
    .short_caption = Use isotropic B-factor model
  pwilson = 0.8
    .type = float
    .help = 'Scale factor for isotropic b-factor model: all atoms = Bwilson * pwilson'
    .short_caption = Scale factor for isotropic b-factor model
  set_occupancies = False
    .type = bool
    .help = 'Set all atoms aoccupancy to 1.0'
    .short_caption = Reset occupancies to 1.0
  target_name = *ml mlhl ls_wunit_k1_fixed ls_wunit_k1
    .type = choice
    .short_caption = Refinement target
    .help = 'Choices for refinement target'
  remove_alt_conf_from_input_pdb = True
    .type = bool
    .help = 'Removes any alternative conformations if present in input PDB model'
  scale_wrt_n_calc_start = False
    .type = bool
    .help = 'Scale <Ncalc> to starting Ncalc'
    .short_caption = Scale <Ncalc> to starting Ncalc
  show_wilson_plot = False
    .type = bool
    .help = 'Print Wilson plot during simulation'
  output_running_kinetic_energy_in_occupancy_column = False
    .type = bool
    .help = 'Output PDB file contains running kinetic energy in place of occupancy'
  ordered_solvent_update = True
    .type = bool
    .help = 'Ordered water molecules automatically updated every nth macro cycle'
  ordered_solvent_update_cycle = 25
    .type = int
    .help = 'Number of macro-cycles per ordered solvent update'
    .short_caption = Solvent update cycles
  harmonic_restraints
    .style = menu_item auto_align box
  {
    selections = None
      .type = atom_selection
      .help = 'Atoms selections to apply harmonic restraints'
    weight = 0.001
      .type = float
      .help = 'Harmonic restraints weight'
    slack = 1.0
      .type = float
      .help = 'Harmonic restraints slack distance'
  }
  electron_density_maps
    .style = menu_item
  {
    apply_default_maps = True
      .type = bool
    include scope mmtbx.maps.map_coeff_params_str
  }
  mask
    .short_caption = Bulk solvent mask
    .style = menu_item auto_align box
  {
    include scope mmtbx.masks.mask_master_params
  }
  ensemble_ordered_solvent
    .style = menu_item auto_align box
  {
    diff_map_cutoff = 2.5
      .type = float
    e_map_cutoff_keep = 0.5
      .type = float
    e_map_cutoff_find = 0.5
      .type = float
    tolerance = 0.9
      .type = float
    ordered_solvent_map_to_model = True
      .type = bool
    include scope mmtbx.solvent.ordered_solvent.output_params_str
    primary_map_type = mFo-DFmodel
      .type=str
    primary_map_cutoff = 3.0
      .type=float
    secondary_map_type = 2mFo-DFmodel
      .type=str
    secondary_map_cutoff_keep = 1.0
      .type=float
    secondary_map_cutoff_find = 1.0
      .type=float
    include scope mmtbx.solvent.ordered_solvent.h_bond_params_str
    include scope mmtbx.solvent.ordered_solvent.adp_occ_params_str
    refine_occupancies = False
      .type = bool
      .help = Refine solvent occupancies.
      .expert_level = 1
    add_hydrogens = False
      .type = bool
      .help = Adds hydrogens to water molecules (except those on special positions)
    refilter = True
      .type = bool
    temperature = 300
      .type = float
      .help = Target temperature for random velocity assignment
    seed = 343534534
      .type = int
      .help = Fixes the random seed for velocity assignment
    preserved_solvent_minimum_distance = 7.0
      .type = float
    find_peaks {
      include scope mmtbx.find_peaks.master_params
    }
  }
  include scope mmtbx.geometry_restraints.external.external_energy_params_str
}
refinement.geometry_restraints.edits
  .short_caption = Edit geometry restraints
  .style = menu_item box auto_align
{
  include scope mmtbx.monomer_library.pdb_interpretation.geometry_restraints_edits_str
}
gui
  .help = Phenix GUI parameters, not used in command-line program
{
  include scope libtbx.phil.interface.tracking_params
  output_dir = None
    .type = path
    .short_caption = Output directory
    .style = output_dir
}
extra_restraints_file = None
  .type = path
  .short_caption = Custom geometry restraints
  .help = File containing custom geometry restraints, using the same format \
    as phenix.refine.  On the command line this can be specified directly \
    as a command-line argument, but this parameter is used by the Phenix GUI.
  .style = input_file file_type:phil
""", process_includes=True).fetch(source=customization_params)

class er_pickle(object):
  def __init__(self,
               pickle_object,
               pickle_filename):
    cPickle.dump(pickle_object, gzip.open(pickle_filename, 'wb'))

class ensemble_refinement_data(object):
  def __init__(self, f_calc_running                      = None,
                     f_calc_data_total                   = None,
                     f_calc_data_current                 = None,
                     f_mask_running                      = None,
                     f_mask_current                      = None,
                     f_mask_total                        = None,
                     total_SF_cntr                       = 0,
                     total_SF_cntr_mask                  = 0,
                     fix_scale_factor                    = None,
                     non_solvent_temp                    = None,
                     solvent_temp                        = None,
                     system_temp                         = None,
                     xray_structures                     = [],
                     pdb_hierarchys                      = [],
                     xray_structures_diff_map            = [],
                     seed                                = None,
                     velocities                          = None,
                     ke_protein_running                  = None,
                     ke_pdb                              = [],
                     geo_grad_rms                        = None,
                     xray_grad_rms                       = None,
                     solvent_sel                         = None,
                     all_sel                             = None,
                     er_harmonic_restraints_info         = None,
                     er_harmonic_restraints_weight       = 0.001,
                     er_harmonic_restraints_slack        = 1.0,
                     macro_cycle                         = None,
                     ):
    adopt_init_args(self, locals())

class er_tls_manager(object):
  def __init__(self, tls_selection_strings_no_sol       = None,
                     tls_selection_strings_no_sol_no_hd = None,
                     tls_selections_with_sol            = None,
                     tls_selections_no_sol              = None,
                     tls_selections_no_sol_no_hd        = None,
                     tls_operators                      = None):
    adopt_init_args(self, locals())

class run_ensemble_refinement(object):
  def __init__(self, fmodel,
                     model,
                     log,
                     raw_data,
                     raw_flags,
                     params,
                     ptls,
                     run_number=None) :
    adopt_init_args(self, locals())
#    self.params = params.extract().ensemble_refinement

    if self.params.target_name in ['ml', 'mlhl'] :
      self.fix_scale = False
    else:
      self.fix_scale = True
    if not self.params.wxray_coupled_tbath:
      self.params.wxray_coupled_tbath_offset = 0.0
    self.wxray = self.params.wxray
    self.params.ensemble_ordered_solvent.temperature = self.params.cartesian_dynamics.temperature
    self.ensemble_utils = ensemble_utils.manager(ensemble_obj = self)
    self.xray_gradient = None
    self.fc_running_ave = None
    self.macro_cycle = 1
    self.sf_model_ave = None
    self.fmodel_total_block_list = []
    self.reset_velocities = True
    self.cmremove = True
    self.cdp = self.params.cartesian_dynamics
    self.bsp = mmtbx.bulk_solvent.bulk_solvent_and_scaling.master_params.extract()
    if (self.params.target_name == 'mlhl') :
      self.bsp.target = 'ml'
    else :
      self.bsp.target = self.params.target_name
    if self.params.tx == None:
      print >> log, "\nAutomatically set Tx (parameter not defined)"
      print >> log, "Tx          :  2(1/dmin)**2"
      self.params.tx = round(2.0 * ((1.0/self.fmodel.f_obs().d_min())**2),1)
      print >> log, 'Dmin        : ', self.fmodel.f_obs().d_min()
      print >> log, 'Set Tx      : ', self.params.tx
    self.n_mc_per_tx = self.params.tx / (self.cdp.time_step * self.cdp.number_of_steps)

    # Set simulation length
    make_header("Simulation length:", out = self.log)
    print >> log, "Number of time steps per macro cycle    : ", self.cdp.number_of_steps
    print >> log, "Tx                                      : ", self.params.tx
    print >> log, "Number macro cycles per Tx period       : ", self.n_mc_per_tx
    self.equilibrium_macro_cycles = int(self.n_mc_per_tx * self.params.equilibrium_n_tx)
    self.acquisition_block_macro_cycles = int(self.n_mc_per_tx * self.params.acquisition_block_n_tx)
    self.total_macro_cycles = int(self.equilibrium_macro_cycles \
                            + (self.acquisition_block_macro_cycles * self.params.number_of_acquisition_periods))
    #
    print >> log, "\nEquilibration"
    print >> log, "    Number Tx periods    : ", self.params.equilibrium_n_tx
    print >> log, "    Number macro cycles  : ", self.equilibrium_macro_cycles
    print >> log, "    Time (ps)            : ", self.equilibrium_macro_cycles \
                                                  * self.cdp.number_of_steps * self.cdp.time_step
    #
    print >> log, "\nAcquisition block"
    print >> log, "    Number Tx periods    : ",  self.params.acquisition_block_n_tx
    print >> log, "    Number macro cycles  : ",  self.acquisition_block_macro_cycles
    print >> log, "    Time (ps)            : ",  self.acquisition_block_macro_cycles \
                                                  * self.cdp.number_of_steps\
                                                  * self.cdp.time_step
    #
    print >> log, "\nSimulation total"
    print >> log, "    Number Tx periods    : ", self.params.equilibrium_n_tx\
                                                + (self.params.number_of_acquisition_periods\
                                                   * self.params.acquisition_block_n_tx)
    print >> log, "    Number macro cycles  : ", self.total_macro_cycles
    self.total_time = self.total_macro_cycles\
                        * self.cdp.number_of_steps\
                        * self.cdp.time_step
    print >> log, "    Time (ps)            : ", self.total_time
    print >> log, "    Total = Equilibration + nAcquisition"
    # Store block
    self.block_store_cycle_cntr = 0
    self.block_store_cycle = \
        range(self.acquisition_block_macro_cycles + self.equilibrium_macro_cycles,
              self.acquisition_block_macro_cycles + self.total_macro_cycles,
              self.acquisition_block_macro_cycles
              )
    # Store pdb
    self.pdb_store_cycle = max(int(self.acquisition_block_macro_cycles \
                         / self.params.pdb_stored_per_block), 1)

    #Setup ensemble_refinement_data_object
    self.er_data = ensemble_refinement_data()
    #Setup fmodels for running average   = refinement target
    #                  total average     = final model
    #                  current model     = model at time point n
    self.fmodel_running = self.fmodel
    self.fmodel_total = None
    self.fmodel_current = None
    self.tls_manager = None
    self.er_data.seed = self.params.random_seed
    self.run_time_stats_dict = {}

    #Dummy miller array
    self.copy_ma = self.fmodel_running.f_masks()[0].array(data = self.fmodel_running.f_masks()[0].data()*0).deep_copy()
    #
    self.fmodel_running.xray_structure = self.model.xray_structure
    assert self.fmodel_running.xray_structure is self.model.xray_structure
    self.pdb_hierarchy = self.model.pdb_hierarchy

    #Atom selections
    self.atom_selections()

    self.model.show_geometry_statistics(message   = "Starting model",
                                        ignore_hd = True,
                                        out       = self.log)

    self.setup_bulk_solvent_and_scale()

    # Wilson plot for input model
    if self.params.show_wilson_plot:
      self.wilson_plot(miller_data = self.fmodel_running.f_obs().data(), header = "Fobs input")
      self.wilson_plot(miller_data = self.fmodel_running.f_calc().data(), header = "Fcalc input")
      self.wilson_plot(miller_data = self.fmodel_running.f_model().data(), header = "Fmodel input")
      self.wilson_plot(miller_data = self.fmodel_running.f_model_scaled_with_k1().data(), header = "Fmodel_k1")

    self.fmodel_running.info(
      free_reflections_per_bin = 100,
      max_number_of_bins       = 999).show_rfactors_targets_in_bins(out = self.log)

    if self.params.target_name in ['ml', 'mlhl'] :
      #Must be called before reseting ADPs
      if self.params.scale_wrt_n_calc_start:
        make_header("Calculate Ncalc and restrain to scale kn", out = self.log)
        self.fmodel_running.n_obs_n_calc(update_nobs_ncalc = True)
        n_obs  = self.fmodel_running.n_obs
        n_calc = self.fmodel_running.n_calc
        self.scale_n1_reference = self.scale_helper(target    = n_calc,
                                                    reference = n_obs
                                                    )
        self.scale_n1_target    = self.scale_n1_reference
        self.scale_n1_current   = self.scale_n1_reference
        self.n_calc_reference = self.fmodel_running.n_calc.deep_copy()
        self.n_mc_per_ncalc_update = max(1, int(self.n_mc_per_tx / 10) )
        print >> self.log, "Number macro cycles per tx     : {0:5.0f}".format(self.n_mc_per_tx)
        print >> self.log, "Number macro cycles per update : {0:5.0f}".format(self.n_mc_per_ncalc_update)
        #
        self.fixed_k1_from_start = self.fmodel_running.scale_k1()
        self.target_k1 = self.fmodel_running.scale_k1()
        self.update_normalisation_factors()
      else:
        make_header("Calculate and fix scale of Ncalc", out = self.log)
        self.fmodel_running.n_obs_n_calc(update_nobs_ncalc = True)
        print >> self.log, "Fix Ncalc scale          : True"
        print >> self.log, "Sum current Ncalc        : {0:5.3f}".format(sum(self.fmodel_running.n_calc))

      # XXX test
      self.wilson_plot(miller_data = self.fmodel_running.n_obs, header = "Nobs reference")
      self.wilson_plot(miller_data = self.fmodel_running.n_calc, header = "Ncalc reference")
      print >> self.log, "|"+"-"*77+"|\n"

    #Set ADP model
    self.tls_manager = er_tls_manager()
    self.setup_tls_selections(tls_group_selection_strings = self.params.tls_group_selections)
    self.fit_tls(input_model = self.model)
    self.assign_solvent_tls_groups()

    #Set occupancies to 1.0
    if self.params.set_occupancies:
      make_header("Set occupancies to 1.0", out = self.log)
      self.model.xray_structure.set_occupancies(
        value      = 1.0)
      self.model.show_occupancy_statistics(out = self.log)

    # Wilson plot for start model
    if self.params.show_wilson_plot:
      self.wilson_plot(miller_data = self.fmodel_running.f_obs().data(), header = "Fobs start")
      self.wilson_plot(miller_data = self.fmodel_running.f_calc().data(), header = "Fcalc start")
      self.wilson_plot(miller_data = self.fmodel_running.f_model().data(), header = "Fmodel start")
      self.wilson_plot(miller_data = self.fmodel_running.f_model_scaled_with_k1().data(), header = "Fmodel_k1 start")

    #Initiates running average SFs
    self.er_data.f_calc_running = self.fmodel_running.f_calc().data().deep_copy()
    #self.fc_running_ave = self.fmodel_running.f_calc()
    self.fc_running_ave = self.fmodel_running.f_calc().deep_copy()

    #Initial sigmaa array, required for ML target function
    #Set eobs and ecalc normalization factors in Fmodel, required for ML
    if self.params.target_name in ['ml', 'mlhl'] :
      self.sigmaa_array = self.fmodel_running.sigmaa().sigmaa().data()
      self.best_r_free = self.fmodel_running.r_free()
      self.fmodel_running.set_sigmaa = self.sigmaa_array

    #Harmonic restraints
    if self.params.harmonic_restraints.selections is not None:
      self.add_harmonic_restraints()

############################## START Simulation ################################
    make_header("Start simulation", out = self.log)
    while self.macro_cycle <= self.total_macro_cycles:
      self.er_data.macro_cycle = self.macro_cycle
      self.time = self.cdp.time_step * self.cdp.number_of_steps * self.macro_cycle
      #XXX Debug
      if False and self.macro_cycle % 10==0:
        print >> self.log, "Sys temp  : ", self.er_data.system_temp
        print >> self.log, "Xray grad : ", self.er_data.xray_grad_rms
        print >> self.log, "Geo grad  : ", self.er_data.geo_grad_rms
        print >> self.log, "Wx        : ", self.wxray

      if self.fmodel_running.target_name in ['ml', 'mlhl'] :
        if self.macro_cycle < self.equilibrium_macro_cycles:
          if self.params.scale_wrt_n_calc_start and self.macro_cycle%self.n_mc_per_ncalc_update == 0:
            self.update_normalisation_factors()
          elif self.macro_cycle%int(self.n_mc_per_tx)==0:
            self.update_normalisation_factors()

      # Ordered Solvent Update
      if self.params.ordered_solvent_update \
          and (self.macro_cycle == 1\
          or self.macro_cycle%self.params.ordered_solvent_update_cycle == 0):
        self.ordered_solvent_update()

      xrs_previous = self.model.xray_structure.deep_copy_scatterers()
      assert self.fmodel_running.xray_structure is self.model.xray_structure

      if self.cdp.verbose >= 1:
        if self.macro_cycle == 1 or self.macro_cycle%100 == 0:
          cdp_verbose = 1
        else:
          cdp_verbose = -1
      else:
        cdp_verbose = -1

      cd_manager = ensemble_cd.cartesian_dynamics(
        structure                   = self.model.xray_structure,
        restraints_manager          = self.model.restraints_manager,
        temperature                 = self.cdp.temperature - self.params.wxray_coupled_tbath_offset,
        protein_thermostat          = self.cdp.protein_thermostat,
        n_steps                     = self.cdp.number_of_steps,
        n_print                     = self.cdp.n_print,
        time_step                   = self.cdp.time_step,
        initial_velocities_zero_fraction = self.cdp.initial_velocities_zero_fraction,
        fmodel                      = self.fmodel_running,
        xray_target_weight          = self.wxray,
        chem_target_weight          = 1.0,
        xray_structure_last_updated = None,
        shift_update                = 0.0,
        xray_gradient               = self.xray_gradient,
        reset_velocities            = self.reset_velocities,
        stop_cm_motion              = self.cmremove,
        update_f_calc               = False,
        er_data                     = self.er_data,
        verbose                     = cdp_verbose,
        log                         = self.log)

      self.reset_velocities = False
      self.cmremove = False

      #Calc rolling average KE energy
      self.kinetic_energy_running_average()
      #Show KE stats
      if self.params.verbose > 0 and self.macro_cycle % 500 == 0:
        self.ensemble_utils.kinetic_energy_stats()

      #Update Fmodel
      self.fmodel_running.update_xray_structure(
        xray_structure      = self.model.xray_structure,
        update_f_calc       = True,
        update_f_mask       = True,
        force_update_f_mask = True)

      #Save current Fmask
      self.er_data.f_mask_current = self.fmodel_running.f_masks()[0].data().deep_copy()

      #Save current Fcalc
      self.er_data.f_calc_data_current = self.fmodel_running.f_calc().data().deep_copy()

      #Total Fmask calculation
      if self.er_data.f_mask_total is None:
        self.er_data.f_mask_total = self.fmodel_running.f_masks()[0].data().deep_copy()
        self.er_data.total_SF_cntr_mask = 1
      else:
        self.er_data.f_mask_total += self.fmodel_running.f_masks()[0].data().deep_copy()
        self.er_data.total_SF_cntr_mask += 1

      #Total Fcalc calculation
      if self.er_data.f_calc_data_total is None:
        self.er_data.f_calc_data_total = self.fmodel_running.f_calc().data().deep_copy()
        self.er_data.total_SF_cntr = 1
      else:
        self.er_data.f_calc_data_total += self.fmodel_running.f_calc().data().deep_copy()
        self.er_data.total_SF_cntr += 1

      #Running average Fcalc calculation
      if self.params.tx > 0:
        self.a_prime = math.exp(-(self.cdp.time_step * self.cdp.number_of_steps)/self.params.tx)
      else:
        self.a_prime = 0

      self.er_data.f_calc_running \
        = (self.a_prime * self.er_data.f_calc_running) + ((1-self.a_prime) * self.fmodel_running.f_calc().data().deep_copy())
      self.fc_running_ave = self.fc_running_ave.array(data = self.er_data.f_calc_running)

      #Update running average Fmask
      if self.macro_cycle == 1:
        self.er_data.f_mask_running = self.fmodel_running.f_masks()[0].data().deep_copy()
      else:
        self.er_data.f_mask_running \
          = (self.a_prime * self.er_data.f_mask_running) + ((1-self.a_prime) * self.fmodel_running.f_masks()[0].data())
      self.running_f_mask_update = self.copy_ma.array(data = self.er_data.f_mask_running).deep_copy()

      #Update runnning average Fcalc and Fmask
      self.fmodel_running.update(f_calc = self.fc_running_ave,
                                 f_mask = self.running_f_mask_update)

      #Update total average Fcalc
      total_f_mask_update \
          = self.copy_ma.array(data = self.er_data.f_mask_total / self.er_data.total_SF_cntr_mask).deep_copy()


      if self.fmodel_total == None:
        self.fmodel_total = self.fmodel_running.deep_copy()
        self.fmodel_total.update(
          f_calc = self.copy_ma.array(data = self.er_data.f_calc_data_total / self.er_data.total_SF_cntr ),
          f_mask = total_f_mask_update)

        if(self.er_data.fix_scale_factor is not None):
          self.fmodel_total.set_scale_switch = self.er_data.fix_scale_factor
      else:
        self.fmodel_total.update(
          f_calc = self.copy_ma.array(data = self.er_data.f_calc_data_total / self.er_data.total_SF_cntr),
          f_mask = total_f_mask_update)

      #Update current time-step Fcalc
      current_f_mask_update = self.copy_ma.array(data = self.er_data.f_mask_current)

      if self.fmodel_current == None:
        self.fmodel_current = self.fmodel_running.deep_copy()
        self.fmodel_current.update(
          f_calc = self.copy_ma.array(data = self.er_data.f_calc_data_current),
          f_mask = current_f_mask_update)
        if(self.er_data.fix_scale_factor is not None):
          self.fmodel_current.set_scale_switch = self.er_data.fix_scale_factor
      else:
        self.fmodel_current.update(
          f_calc = self.copy_ma.array(data = self.er_data.f_calc_data_current),
          f_mask = current_f_mask_update)

      #ML params update
      if self.params.target_name in ['ml', 'mlhl'] :
        if self.macro_cycle < self.equilibrium_macro_cycles:
          if self.fmodel_running.r_free() < (self.best_r_free - self.params.update_sigmaa_rfree):
            self.update_sigmaa()

      # XXX wilson plot
      if self.params.show_wilson_plot and self.macro_cycle%int(self.n_mc_per_tx) == 0:
        self.wilson_plot(miller_data = self.fmodel_running.f_obs().data(), header = "Fobs running")
        self.wilson_plot(miller_data = self.fmodel_running.f_calc().data(), header = "Fcalc running")
        self.wilson_plot(miller_data = self.fmodel_running.f_model().data(), header = "Fmodel running")
        self.wilson_plot(miller_data = self.fmodel_running.f_model_scaled_with_k1().data(), header = "Fmodel_k1 running")

      # Wxray coupled to temperature bath
      if self.params.wxray_coupled_tbath:
        if self.macro_cycle < 5:
          self.wxray        = 2.5
        elif self.macro_cycle < self.equilibrium_macro_cycles:
          if self.params.tx == 0:
            a_prime_wx = 0
          else:
            wx_tx = min(self.time, self.params.tx)
            a_prime_wx = math.exp(-(self.cdp.time_step * self.cdp.number_of_steps)/wx_tx)
          wxray_t = self.wxray * max(0.01, self.cdp.temperature / self.er_data.non_solvent_temp)
          self.wxray = (a_prime_wx * self.wxray) + ((1-a_prime_wx) * wxray_t)

      #Store current structure, current KE
      if self.macro_cycle % self.pdb_store_cycle == 0 \
           and self.macro_cycle >= self.equilibrium_macro_cycles:
        self.er_data.xray_structures.append(self.model.xray_structure.deep_copy_scatterers())
        self.er_data.pdb_hierarchys.append(self.model.pdb_hierarchy().deep_copy())
        if self.er_data.ke_protein_running is None:
          self.er_data.ke_pdb.append(flex.double(self.model.xray_structure.sites_cart().size(), 0.0) )
        else:
          ke_expanded = flex.double(self.model.xray_structure.sites_cart().size(), 0.0)
          ke_expanded.set_selected(~self.model.solvent_selection(),
                                   self.er_data.ke_protein_running)
          self.er_data.ke_pdb.append(ke_expanded)

      #Current structural deviation vs starting structure and previous macro-cycle structure
      if xrs_previous.distances(other = self.model.xray_structure).min_max_mean().mean > 1.0:
        print >> self.log, "\n\nWARNING:"
        print >> self.log, "Macro cycle too long, max atomic deviation w.r.t. previous cycle"
        print >> self.log, "greater than 1.0A"
        print >> self.log, "Reduce params.cartesian_dynamics.number_of_steps"
        print >> self.log, "Max deviation : {0:1.3f}"\
          .format(xrs_previous.distances(other = self.model.xray_structure).min_max_mean().mean)

      if self.fmodel_running.r_work() > 0.75:
        raise Sorry("Simulation aborted, running Rfree > 75%")

      #Print run time stats
      if self.macro_cycle == 1 or self.macro_cycle%50 == 0:
        print >> self.log, "\n________________________________________________________________________________"
        print >> self.log, "    MC        Time     |  Current  |  Rolling  |   Total   | Temp |  Grad Wxray "
        print >> self.log, "          (ps)     (%) |   Rw   Rf |   Rw   Rf |   Rw   Rf |  (K) |   X/G       "
      print >> self.log, \
          "~{0:5d} {1:7.2f} {2:7.2f} | {3:4.1f} {4:4.1f} | {5:4.1f} {6:4.1f} | {7:4.1f} {8:4.1f} | {9:4.0f} | {10:5.2f} {11:5.2f}"\
          .format(self.macro_cycle,
                  self.time,
                  100 * self.time / self.total_time,
                  100*self.fmodel_current.r_work(),
                  100*self.fmodel_current.r_free(),
                  100*self.fmodel_running.r_work(),
                  100*self.fmodel_running.r_free(),
                  100*self.fmodel_total.r_work(),
                  100*self.fmodel_total.r_free(),
                  self.er_data.non_solvent_temp,
                  self.er_data.xray_grad_rms / self.er_data.geo_grad_rms,
                  self.wxray)

      if self.params.verbose > 0:
        if self.macro_cycle == 1\
            or self.macro_cycle%100 == 0\
            or self.macro_cycle == self.total_macro_cycles:
          self.print_fmodels_scale_and_solvent_stats()

      if self.params.number_of_acquisition_periods > 1:
        if self.macro_cycle in self.block_store_cycle:
          self.save_multiple_fmodel()

      #End of equilibration period, reset total structure factors, atomic cords, kinetic energies
      if self.macro_cycle == self.equilibrium_macro_cycles:
        self.reset_totals()
      #
      assert self.model.xray_structure is cd_manager.structure
      assert self.fmodel_running.xray_structure is cd_manager.structure
      if self.fix_scale == True:
        assert self.fmodel_running.scale_k1() == self.er_data.fix_scale_factor
      self.macro_cycle +=1

############################## END Simulation ##################################

    self.macro_cycle = self.total_macro_cycles
    #Find optimum section of acquisition period
    if self.params.number_of_acquisition_periods > 1:
      self.optimise_multiple_fmodel()
    else:
      self.fmodel_total.set_scale_switch = 0
      self.fmodel_total.update_solvent_and_scale(
                            verbose       = self.params.verbose,
                            out           = self.log,
                            params        = self.bsp,
                            optimize_mask = False)

    #Minimize number of ensemble models
    if self.params.ensemble_reduction:
      self.ensemble_utils.ensemble_reduction()

    #Optimise fmodel_total k, b_aniso, k_sol, b_sol
    self.fmodel_total.set_scale_switch = 0
    self.print_fmodels_scale_and_solvent_stats()
    self.fmodel_total.update_solvent_and_scale(
                          verbose       = self.params.verbose,
                          out           = self.log,
                          params        = self.bsp,
                          optimize_mask = False
                          )
    self.print_fmodels_scale_and_solvent_stats()
    print >> self.log, "FINAL Rwork = %6.4f Rfree = %6.4f Rf/Rw = %6.4f"\
        %(self.fmodel_total.r_work(),
          self.fmodel_total.r_free(),
          self.fmodel_total.r_free() / self.fmodel_total.r_work()
          )
    print >> self.log, "Final Twork = %6.4f Tfree = %6.4f Tf/Tw = %6.4f"\
        %(self.fmodel_total.target_w(),
          self.fmodel_total.target_t(),
          self.fmodel_total.target_t() / self.fmodel_total.target_w()
          )
    info = self.fmodel_total.info(free_reflections_per_bin = 100,
                                  max_number_of_bins       = 999
                                  )
    info.show_remark_3(out = self.log)
    info.show_rfactors_targets_in_bins(out = self.log)

    # Final wilson plot
    if self.params.show_wilson_plot:
      self.wilson_plot(miller_data = self.fmodel_total.f_obs().data(), header = "Fobs final")
      self.wilson_plot(miller_data = self.fmodel_total.f_calc().data(), header = "Fcalc final")
      self.wilson_plot(miller_data = self.fmodel_total.f_model().data(), header = "Fmodel final")
      self.wilson_plot(miller_data = self.fmodel_total.f_model_scaled_with_k1().data(), header = "Fmodel_k1 final")

    self.write_output_files(run_number=run_number)

############################## END ER ##########################################

  def write_output_files (self, run_number=None) :
    #PDB output
    prefix = self.params.output_file_prefix
    if (run_number is not None) :
      prefix += "_%g" % run_number
    pdb_out = prefix + ".pdb"
    cif_out = prefix + ".cif"
    if (self.params.gzip_final_model) :
      pdb_out += ".gz"
      self.write_ensemble_pdb(out = gzip.open(pdb_out, 'wb'))
      # TODO
      if False :#(self.params.write_cif_file) :
        self.write_ensemble_mmcif(out=gzip.open(cif_out, 'wb'))
    else :
      self.write_ensemble_pdb(out = open(pdb_out, 'wb'))
      # TODO
      if False :#(self.params.write_cif_file) :
        self.write_ensemble_mmcif(out=open(cif_out, 'wb'))
    self.pdb_file = pdb_out
    # Map output
    assert (self.fmodel_total is not None)
    self.mtz_file = write_mtz_file(
      fmodel_total=self.fmodel_total,
      raw_data=self.raw_data,
      raw_flags=self.raw_flags,
      prefix=prefix,
      params=self.params)

  def show_overall(self, message = "", fmodel_running = True):
    if fmodel_running:
      message = "Running: " + message
      self.fmodel_running.info().show_rfactors_targets_scales_overall(header = message, out = self.log)
    else:
      message = "Total: " + message
      self.fmodel_total.info().show_rfactors_targets_scales_overall(header = message, out = self.log)

  def add_harmonic_restraints(self):
    make_header("Add specific harmonic restraints", out = self.log)
    # ensures all solvent atoms are at the end prior to applying harmonic restraints
    self.ordered_solvent_update()
    all_chain_proxies = self.generate_all_chain_proxies(model = self.model)
    hr_selections = mmtbx.utils.get_atom_selections(
        all_chain_proxies = all_chain_proxies,
        selection_strings = self.params.harmonic_restraints.selections,
        xray_structure    = self.model.xray_structure)
    pdb_atoms = self.pdb_hierarchy().atoms()
    print >> self.log, "\nAdd atomic harmonic restraints:"
    restraint_info = []
    for i_seq in hr_selections[0]:
      atom_info = pdb_atoms[i_seq].fetch_labels()
      print >> self.log, '    {0} {1} {2} {3} {4}     '.format(
                                   atom_info.name,
                                   atom_info.i_seq+1,
                                   atom_info.resseq,
                                   atom_info.resname,
                                   atom_info.chain_id,
                                   )
      restraint_info.append((i_seq, pdb_atoms[i_seq].xyz))
    self.er_data.er_harmonic_restraints_info = restraint_info
    self.er_data.er_harmonic_restraints_weight = self.params.harmonic_restraints.weight
    self.er_data.er_harmonic_restraints_slack  = self.params.harmonic_restraints.slack
    print >> self.log, "\n|"+"-"*77+"|\n"

  def setup_bulk_solvent_and_scale(self):
    make_header("Setup bulk solvent and scale", out = self.log)
    self.show_overall(message = "pre solvent and scale")
    #
    self.fmodel_running.update_solvent_and_scale(
        params        = self.bsp,
        verbose       = self.params.verbose,
        out           = self.log,
        optimize_mask = True)

    #Fixes scale factor for rolling average #ESSENTIAL for LSQ
    if self.fix_scale == True:
      self.er_data.fix_scale_factor = self.fmodel_running.scale_k1()
      self.fmodel_running.set_scale_switch = self.er_data.fix_scale_factor
    self.show_overall(message = "post solvent and scale")

  def scale_helper(self, reference, target):
    return flex.sum(reference * target) / flex.sum(flex.pow2(target))

  def update_normalisation_factors(self):
    if self.params.scale_wrt_n_calc_start:
      # Adaptive scaling
      # Ncalc_start / Ncalc_current
      make_header("Update Ncalc and restrain to Ncalc ref", out = self.log)
      # Get N_calc current, compare with reference
      n_obs, n_calc =\
        self.fmodel_running.n_obs_n_calc(update_nobs_ncalc = False)
      ref_div_current = self.n_calc_reference / n_calc

      n_calc_coeff    = 1.0-math.exp(-self.n_mc_per_ncalc_update/self.n_mc_per_tx)
      n_calc_scaled   = ref_div_current * n_calc_coeff
      n_calc_update   = (self.fmodel_running.n_calc * (1.0-n_calc_coeff) ) + (self.fmodel_running.n_calc * ref_div_current * n_calc_coeff)

      # Update with scaled array
      self.fmodel_running.n_calc = n_calc_update

    else:
      # Normalise to reference Sum(Ncalc)
      make_header("Update and renormalise Ncalc array", out = self.log)
      eobs_norm_factor, ecalc_norm_factor =\
        self.fmodel_running.n_obs_n_calc(update_nobs_ncalc = False)
      self.scale_n1_current = self.scale_helper(target    = ecalc_norm_factor,
                                                reference = eobs_norm_factor
                                                )
      print >> self.log, "Kn current               : {0:5.3f}".format(self.scale_n1_current)
      ecalc_k = sum(self.fmodel_running.n_calc) / sum(ecalc_norm_factor)
      ecalc_k_alt = flex.sum(self.fmodel_running.n_calc * ecalc_norm_factor) / flex.sum(flex.pow2(ecalc_norm_factor) )
      print >> self.log, "Sum current Ncalc        : {0:5.3f}".format(sum(self.fmodel_running.n_calc) )
      print >> self.log, "Sum updated Ncalc        : {0:5.3f}".format(sum(ecalc_norm_factor) )
      print >> self.log, "Rescaling factor         : {0:5.3f}".format(ecalc_k)
      print >> self.log, "Rescaling factor alt     : {0:5.3f}".format(ecalc_k_alt)
      ecalc_norm_factor = ecalc_k * ecalc_norm_factor
      self.fmodel_running.n_calc = ecalc_norm_factor
    print >> self.log, "|"+"-"*77+"|\n"

  def update_sigmaa(self):
    make_header("Update sigmaa", out = self.log)
    if self.params.verbose > 0:
      print >> self.log, "Previous best Rfree      : ", self.best_r_free
      print >> self.log, "Current       Rfree      : ", self.fmodel_running.r_free()
      self.print_ml_stats()
      print >> self.log, "  Update sigmaa"
    self.sigmaa_array = self.fmodel_running.sigmaa().sigmaa().data()
    self.fmodel_running.set_sigmaa = self.sigmaa_array
    if self.params.verbose > 0:
      self.print_ml_stats()
    self.best_r_free = self.fmodel_running.r_free()
    print >> self.log, "|"+"-"*77+"|\n"

  def setup_tls_selections(self, tls_group_selection_strings):
    make_header("Generating TLS selections from input parameters (not including solvent)", out = self.log)
    model_no_solvent = self.model.deep_copy()
    model_no_solvent = model_no_solvent.remove_solvent()
    all_chain_proxies = self.generate_all_chain_proxies(model = model_no_solvent)

    if len(tls_group_selection_strings) < 1:
      print >> self.log, '\nNo TLS groups supplied - automatic setup'
      # Get chain information
      chains_info = []
      for chain in model_no_solvent.pdb_hierarchy().chains():
        count_h = 0
        for atom in chain.atoms():
          if atom.element_is_hydrogen(): count_h+=1
        chain_id_non_h = (chain.id, chain.atoms_size() - count_h)
        chains_info.append(chain_id_non_h)
      # Check all chains > 63 heavy atoms for TLS fitting
      chains_size = flex.int(zip(*chains_info)[1])
      chains_size_ok = flex.bool(chains_size > 63)
      if sum(chains_size) < 63:
        print >> self.log, '\nStructure contains less than 63 atoms (non H/D, non solvent)'
        print >> self.log, '\nUnable to perform TLS fitting, will use isotropic B-factor model'
      elif chains_size_ok.count(False) == 0:
        print >> self.log, '\nTLS selections:'
        print >> self.log, 'Chain, number atoms (non H/D)'
        for chain in chains_info:
          tls_group_selection_strings.append('chain ' + chain[0])
          print >> self.log, chain[0], chain[1]
      else:
        print >> self.log, '\nFollowing chains contain less than 63 atoms (non H/D):'
        tls_group_selection_strings.append('chain ')
        for chain in chains_info:
          tls_group_selection_strings[0] += (chain[0] + ' or chain ')
          if chain[1] < 63:
            print >> self.log, chain[0], chain[1]
        print >> self.log, 'Combining all chains to single TLS group'
        print >> self.log, 'WARNING: this may not be the optimum tls groupings to use'
        print >> self.log, 'TLS selections:'
        tls_group_selection_strings[0] = tls_group_selection_strings[0][0:-10]
        print >> self.log, tls_group_selection_strings[0]
    #
    tls_no_sol_selections =  mmtbx.utils.get_atom_selections(
        all_chain_proxies = all_chain_proxies,
        selection_strings = tls_group_selection_strings,
        xray_structure    = model_no_solvent.xray_structure)
    #
    tls_no_hd_selection_strings = []
    for selection_string in tls_group_selection_strings:
      no_hd_string = '(' + selection_string + ') and not (element H or element D)'
      tls_no_hd_selection_strings.append(no_hd_string)

    tls_no_sol_no_hd_selections = mmtbx.utils.get_atom_selections(
        all_chain_proxies = all_chain_proxies,
        selection_strings = tls_no_hd_selection_strings,
        xray_structure    = model_no_solvent.xray_structure)

    #
    assert self.tls_manager is not None
    self.tls_manager.tls_selection_strings_no_sol       = tls_group_selection_strings
    self.tls_manager.tls_selection_strings_no_sol_no_hd = tls_no_hd_selection_strings
    self.tls_manager.tls_selections_no_sol              = tls_no_sol_selections
    self.tls_manager.tls_selections_no_sol_no_hd        = tls_no_sol_no_hd_selections
    self.tls_manager.tls_operators = mmtbx.tls.tools.generate_tlsos(
        selections     = self.tls_manager.tls_selections_no_sol,
        xray_structure = model_no_solvent.xray_structure,
        value          = 0.0)

    self.model.tls_groups = mmtbx.tls.tools.tls_groups(
        selection_strings = self.tls_manager.tls_selection_strings_no_sol,
        tlsos             = self.tls_manager.tls_operators)

  def generate_all_chain_proxies(self, model = None):
    if model == None:
      model = self.model
    raw_records = [pdb.format_cryst1_record(crystal_symmetry=self.model.xray_structure.crystal_symmetry())]
    pdb_hierarchy = model.pdb_hierarchy
    raw_records.extend(pdb_hierarchy().as_pdb_string().splitlines())
    pip = model.processed_pdb_files_srv.pdb_interpretation_params
    pip.clash_guard.nonbonded_distance_threshold = -1.0
    pip.clash_guard.max_number_of_distances_below_threshold = 100000000
    pip.clash_guard.max_fraction_of_distances_below_threshold = 1.0
    pip.proceed_with_excessive_length_bonds=True
    model.processed_pdb_files_srv.pdb_interpretation_params.\
        clash_guard.nonbonded_distance_threshold=None
    processed_pdb_file, pdb_inp = model.processed_pdb_files_srv.\
      process_pdb_files(raw_records = raw_records)
    return processed_pdb_file.all_chain_proxies

  def fit_tls(self, input_model, verbose = False):
    make_header("Fit TLS from reference model", out = self.log)
    model_copy = input_model.deep_copy()
    model_copy = model_copy.remove_solvent()
    print >> self.log, 'Reference model :'
    model_copy.show_adp_statistics(padded = True, out = self.log)
    start_xrs = model_copy.xray_structure.deep_copy_scatterers()
    start_xrs.convert_to_isotropic()
    start_biso = start_xrs.scatterers().extract_u_iso()/adptbx.b_as_u(1)
    model_copy.xray_structure.convert_to_anisotropic()
    tls_selection_no_sol_hd            = self.tls_manager.tls_selections_no_sol_no_hd
    tls_selection_no_sol_hd_exclusions = self.tls_manager.tls_selections_no_sol_no_hd
    pre_fitted_mean = 999999.99
    #
    use_isotropic = False
    for group in self.tls_manager.tls_selections_no_sol_no_hd:
      if group.size() < 63:
        self.params.isotropic_b_factor_model = True
      elif self.ptls * group.size() < 63:
        self.ptls = 64.0 / group.size()
        print >> self.log, '\nAutomatically increasing pTLS to : {0:5.3f}'.format(self.ptls)
    if self.params.isotropic_b_factor_model:
      print >> self.log, '\nModel contains less than 63 non-solvent, non-H/D atoms'
      print >> self.log, 'Insufficient to fit TLS model, using isotropic model'
      iso_b  = self.fmodel_running.wilson_b() * self.params.pwilson
      episq = 8.0*(math.pi**2)
      print >> self.log, 'Isotropic translation (B) : {0:5.3f}'.format(iso_b)
      print >> self.log, '  = Wilson b-factor * pwilson'
      iso_u = iso_b / episq
      print >> self.log, 'Isotropic translation (U) : {0:5.3f}'.format(iso_u)
      fit_tlsos = []
      for tls_group in self.tls_manager.tls_operators:
        tls_t_new = (iso_u,
                     iso_u,
                     iso_u,
                     0.0,
                     0.0,
                     0.0)
        tls_l_new = (0.0,0.0,0.0,0.0,0.0,0.0)
        tls_s_new = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        fit_tlsos.append(tls_tools.tlso(t      = tls_t_new,
                                      l      = tls_l_new,
                                      s      = tls_s_new,
                                      origin = tls_group.origin))
      self.tls_manager.tls_operators = fit_tlsos
      for tls_group in self.tls_manager.tls_operators:
        mmtbx.tls.tools.show_tls_one_group(tlso = tls_group,
                                           out  = self.log)

    else:
      for fit_cycle in xrange(self.params.max_ptls_cycles):
        fit_tlsos = mmtbx.tls.tools.generate_tlsos(
          selections     = tls_selection_no_sol_hd_exclusions,
          xray_structure = model_copy.xray_structure,
          value          = 0.0)
        print >> self.log, '\nFitting cycle : ', fit_cycle+1
        for rt,rl,rs in [[1,0,1],[1,1,1],[0,1,1],
                         [1,0,0],[0,1,0],[0,0,1],[1,1,1],
                         [0,0,1]]*10:
          fit_tlsos = mmtbx.tls.tools.tls_from_uanisos(
            xray_structure               = model_copy.xray_structure,
            selections                   = tls_selection_no_sol_hd_exclusions,
            tlsos_initial                = fit_tlsos,
            number_of_macro_cycles       = 10,
            max_iterations               = 100,
            refine_T                     = rt,
            refine_L                     = rl,
            refine_S                     = rs,
            enforce_positive_definite_TL = True,
            verbose                      = -1,
            out                          = self.log)
        fitted_tls_xrs = model_copy.xray_structure.deep_copy_scatterers()
        us_tls = mmtbx.tls.tools.u_cart_from_tls(
               sites_cart = fitted_tls_xrs.sites_cart(),
               selections = self.tls_manager.tls_selections_no_sol,
               tlsos      = fit_tlsos)
        fitted_tls_xrs.set_u_cart(us_tls)
        fitted_tls_xrs.convert_to_isotropic()
        fitted_biso = fitted_tls_xrs.scatterers().extract_u_iso()/adptbx.b_as_u(1)
        mmtbx.tls.tools.show_tls(tlsos = fit_tlsos, out = self.log)
        #For testing
        if verbose:
          pdb_hierarchy = model_copy.pdb_hierarchy
          pdb_atoms = pdb_hierarchy().atoms()
          not_h_selection = pdb_hierarchy().atom_selection_cache().selection('not element H')
          ca_selection = pdb_hierarchy().atom_selection_cache().selection('name ca')
          print >> self.log, '\nCA atoms (Name/res number/res name/chain/atom number/ref biso/fit biso:'
          for i_seq, ca in enumerate(ca_selection):
            if ca:
              atom_info = pdb_atoms[i_seq].fetch_labels()
              print >> self.log, atom_info.name, atom_info.resseq, atom_info.resname, atom_info.chain_id, " | ", i_seq, start_biso[i_seq], fitted_biso[i_seq]

        delta_ref_fit = start_biso - fitted_biso
        hd_selection = model_copy.xray_structure.hd_selection()
        delta_ref_fit_no_h = delta_ref_fit.select(~hd_selection)
        delta_ref_fit_no_h_basic_stats = scitbx.math.basic_statistics(delta_ref_fit_no_h )
        start_biso_no_hd = start_biso.select(~hd_selection)
        fitted_biso_no_hd = fitted_biso.select(~hd_selection)

        if verbose:
          print >> self.log, 'pTLS                                    : ', self.ptls

        sorted_delta_ref_fit_no_h = sorted(delta_ref_fit_no_h)
        percentile_cutoff = sorted_delta_ref_fit_no_h[int(len(sorted_delta_ref_fit_no_h) * self.ptls)-1]
        if verbose:
          print >> self.log, 'Cutoff (<)                              : ', percentile_cutoff

        print >> self.log, 'Number of atoms (non HD)                : ', delta_ref_fit_no_h.size()
        delta_ref_fit_no_h_include = flex.bool(delta_ref_fit_no_h < percentile_cutoff)
        print >> self.log, 'Number of atoms (non HD) used in fit    : ', delta_ref_fit_no_h_include.count(True)
        print >> self.log, 'Percentage (non HD) used in fit         : {0:5.3f}'.format(delta_ref_fit_no_h_include.count(True) / delta_ref_fit_no_h.size())

        # Convergence test
        if fitted_biso_no_hd.min_max_mean().mean == pre_fitted_mean:
          break
        else:
          pre_fitted_mean = fitted_biso_no_hd.min_max_mean().mean

        # N.B. map on to full array including hydrogens for i_seqs
        include_array = flex.bool(delta_ref_fit < percentile_cutoff)
        #
        include_i_seq = []
        assert delta_ref_fit.size() == model_copy.xray_structure.sites_cart().size()
        assert include_array.size() == model_copy.xray_structure.sites_cart().size()
        for i_seq, include_flag in enumerate(include_array):
          if include_flag and not hd_selection[i_seq]:
            include_i_seq.append(i_seq)
        tls_selection_no_sol_hd_exclusions = []
        for group in xrange(len(tls_selection_no_sol_hd)):
          new_group = flex.size_t()
          for x in tls_selection_no_sol_hd[group]:
            if x in include_i_seq:
              new_group.append(x)
          if len(new_group) < 63:
            raise Sorry("Number atoms in TLS too small; increase size of group or reduce cut-off")
          if verbose:
            print >> self.log, 'TLS group ', group+1, ' number atoms ', len(new_group)
          tls_selection_no_sol_hd_exclusions.append(new_group)

    print >> self.log, '\nFinal non-solvent b-factor model'
    model_copy.xray_structure.convert_to_anisotropic()

    us_tls = mmtbx.tls.tools.u_cart_from_tls(
             sites_cart = model_copy.xray_structure.sites_cart(),
             selections = self.tls_manager.tls_selections_no_sol_no_hd,
             tlsos      = fit_tlsos)
    model_copy.xray_structure.set_u_cart(us_tls)
    model_copy.show_adp_statistics(padded = True, out = self.log)
    del model_copy

    #Update TLS params
    self.model.tls_groups.tlsos = fit_tlsos
    self.tls_manager.tls_operators = fit_tlsos
    self.assign_solvent_tls_groups()

  def tls_parameters_update(self):
    self.model.xray_structure.convert_to_anisotropic()
    #Apply TLS w.r.t. to atomic position
    selections = self.tls_manager.tls_selections_with_sol
    us_tls = mmtbx.tls.tools.u_cart_from_tls(
               sites_cart = self.model.xray_structure.sites_cart(),
               selections = selections,
               tlsos      = self.tls_manager.tls_operators)
    for selection in selections:
      self.model.xray_structure.set_u_cart(us_tls, selection = selection)
    self.fmodel_running.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = False,
      update_f_mask  = False)

  def assign_solvent_tls_groups(self):
    self.model.xray_structure.convert_to_anisotropic(selection =  self.model.solvent_selection())
    self.fmodel_running.update_xray_structure(
      xray_structure  = self.model.xray_structure,
      update_f_calc   = False,
      update_f_mask   = False)
    #
    self.tls_manager.tls_selections_with_sol = []
    for grp in self.tls_manager.tls_selections_no_sol:
      self.tls_manager.tls_selections_with_sol.append(grp.deep_copy())
    #
    if len(self.tls_manager.tls_selections_with_sol) == 1:
      pdb_atoms     = self.pdb_hierarchy().atoms()
      hoh_selection = self.model.solvent_selection()
      for n, atom in enumerate(pdb_atoms):
        if hoh_selection[n]:
          self.tls_manager.tls_selections_with_sol[0].append(n)
    else:
      model             = self.model.deep_copy()
      solvent_selection = model.solvent_selection()
      solvent_xrs       = model.xray_structure.select(solvent_selection)
      model             = model.remove_solvent()
      closest_distances = model.xray_structure.closest_distances(
                              sites_frac      = solvent_xrs.sites_frac(),
                              use_selection   = ~model.xray_structure.hd_selection(),
                              distance_cutoff = 10.0)
      assert len(solvent_xrs.sites_cart()) == len(closest_distances.i_seqs)
      number_non_solvent_atoms = model.xray_structure.sites_cart().size()
      for n, i_seq in enumerate(closest_distances.i_seqs):
        for grp in self.tls_manager.tls_selections_with_sol:
          if i_seq in grp:
            grp.append(n+number_non_solvent_atoms)
            break
    #
    self.tls_parameters_update()

  def kinetic_energy_running_average(self):
    #Kinetic energy
    atomic_weights = self.model.xray_structure.atomic_weights()
    ke = 0.5 * atomic_weights * self.er_data.velocities.dot()
    #Select non-solvent atoms
    ke = ke.select(~self.model.solvent_selection() )
    if self.er_data.ke_protein_running == None:
      self.er_data.ke_protein_running = ke
    else:
      self.er_data.ke_protein_running\
        = (self.a_prime * self.er_data.ke_protein_running) + ( (1-self.a_prime) * ke)

  def ordered_solvent_update(self):
    ensemble_ordered_solvent_manager = ensemble_ordered_solvent.manager(
        model             = self.model,
        fmodel            = self.fmodel_running,
        verbose           = self.params.verbose,
        params            = self.params.ensemble_ordered_solvent,
        velocities        = self.er_data.velocities,
        log               = self.log)
    self.model = ensemble_ordered_solvent_manager.model
    self.er_data.velocities = ensemble_ordered_solvent_manager.velocities
    self.fmodel_running.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = False,
      update_f_mask  = False)
    assert self.fmodel_running.xray_structure is self.model.xray_structure
    self.xray_gradient = None
    #Update atom selections
    self.pdb_hierarchy = self.model.pdb_hierarchy
    self.atom_selections()
    #Reset solvent tls groups
    if self.tls_manager is not None:
      self.assign_solvent_tls_groups()

  def reset_totals(self):
    make_header("Reseting structure ensemble and total Fmodel", out = self.log)
    self.er_data.xray_structures = []
    self.er_data.xray_structures_diff_map = []
    self.er_data.pdb_hierarchys = []
    self.er_data.ke_pdb = []
    self.er_data.f_calc_data_total = None
    self.er_data.total_SF_cntr = 0
    self.er_data.f_mask_total = None
    self.er_data.total_SF_cntr_mask = 0

  #Generates list of atom selections (needed to pass to CD)
  def atom_selections(self):
    self.er_data.all_sel     = flex.bool(self.model.xray_structure.sites_cart().size(), True)
    self.er_data.solvent_sel = self.model.solvent_selection()

  def save_multiple_fmodel(self):
    make_header("Saving fmodel block", out = self.log)
    #Stores fcalc, fmask, xray structure, pdb hierarchys
    print >> self.log, '{0:<23}: {1:>8} {2:>8} {3:>8} {4:>8}'.format('','MC','Block','Rwork','Rfree')
    print >> self.log, "{0:<23}: {1:8d} {2:8d} {3:8.3f} {4:8.3f}".format(
        'Fmodel block info',
        self.macro_cycle,
        self.block_store_cycle_cntr+1,
        100 * self.fmodel_total.r_work(),
        100 * self.fmodel_total.r_free() )
    fcalc_block  = self.er_data.f_calc_data_total / self.er_data.total_SF_cntr
    fmask_block  = self.er_data.f_mask_total / self.er_data.total_SF_cntr_mask
    xrs_block    = self.er_data.xray_structures
    xrs_dm_block = self.er_data.xray_structures_diff_map
    pdb_h_block  = self.er_data.pdb_hierarchys
    ke_pdb_block = self.er_data.ke_pdb

    block_info = (fcalc_block,
                  fmask_block,
                  xrs_block,
                  xrs_dm_block,
                  pdb_h_block,
                  ke_pdb_block)

    self.block_store_cycle_cntr+1
    if self.block_store_cycle_cntr+1 == 1:
      self.block_temp_file_list = []
    prefix = self.params.output_file_prefix
    if (self.run_number is not None) :
      prefix += "_%g" % self.run_number
    filename = str(self.block_store_cycle_cntr+1)+'_block_'+prefix+'_TEMP.pZ'
    self.block_temp_file_list.append(filename)
    er_pickle(pickle_object = block_info, pickle_filename = filename)
    self.block_store_cycle_cntr += 1
    if self.macro_cycle != self.total_macro_cycles:
      self.reset_totals()

  def optimise_multiple_fmodel(self):
    make_header("Block selection by Rwork", out = self.log)
    best_r_work = None

    # Load all temp files
    self.fmodel_total_block_list = []
    for filename in self.block_temp_file_list:
      block_info = cPickle.load(gzip.open(filename,'rb'))
      self.fmodel_total_block_list.append(block_info)
      os.remove(filename)

    self.fmodel_total.set_scale_switch = 0
    print >> self.log, '  {0:>17} {1:>8} {2:>8}'\
      .format('Block range','Rwork','Rfree','k1')
    for x in xrange(len(self.fmodel_total_block_list)):
      x2 = x+1
      y = len(self.fmodel_total_block_list)
      while y > x:
        fcalc_tot = self.fmodel_total_block_list[x][0].deep_copy()
        fmask_tot = self.fmodel_total_block_list[x][1].deep_copy()
        cntr      = 1
        while x2 < y:
          fcalc_tot += self.fmodel_total_block_list[x2][0].deep_copy()
          fmask_tot += self.fmodel_total_block_list[x2][1].deep_copy()
          x2     += 1
          cntr   += 1
        self.fmodel_total.update(
          f_calc = self.copy_ma.array(data = (fcalc_tot / cntr)),
          f_mask = self.copy_ma.array(data = (fmask_tot / cntr)) )
        self.fmodel_total.update_solvent_and_scale(
          verbose = self.params.verbose,
          params = self.bsp,
          out = self.log,
          optimize_mask = False)
        print >> self.log, "  {0:8d} {1:8d} {2:8.3f} {3:8.3f}"\
          .format(x+1,
                  y,
                  self.fmodel_total.r_work(),
                  self.fmodel_total.r_free(),
                  self.fmodel_total.scale_k1()
                  )
        if best_r_work == None:
          best_r_work = self.fmodel_total.r_work()
          best_r_work_block = [x,y]
          best_r_work_fcalc = (fcalc_tot / cntr)
          best_r_work_fmask = (fmask_tot / cntr)
        elif self.fmodel_total.r_work() < (best_r_work - 0.01):
          best_r_work = self.fmodel_total.r_work()
          best_r_work_block = [x,y]
          best_r_work_fcalc = (fcalc_tot / cntr)
          best_r_work_fmask = (fmask_tot / cntr)
        x2 = x+1
        y -= 1
    self.fmodel_total.update(
      f_calc = self.copy_ma.array(data = best_r_work_fcalc),
      f_mask = self.copy_ma.array(data = best_r_work_fmask) )
    self.fmodel_total.update_solvent_and_scale(
          verbose       = self.params.verbose,
          params        = self.bsp,
          out           = self.log,
          optimize_mask = False)

    print >> self.log, "\nOptimium block :"
    print >> self.log, "  {0:8d} {1:8d} {2:8.3f} {3:8.3f} {4:8.3f} {5:8.3f}"\
      .format(best_r_work_block[0]+1,
              best_r_work_block[1],
              self.fmodel_total.r_work(),
              self.fmodel_total.r_free(),
              self.fmodel_total.scale_k1(),
              self.fmodel_total.fmodel_kbu().k_sols()[0],
              self.fmodel_total.fmodel_kbu().b_sol() )
    #Update self.er_data.xray_structures and self.er_data.pdb_hierarchys to correspond to optimum fmodel_total
    self.er_data.xray_structures = []
    self.er_data.xray_structures_diff_map =[]
    self.er_data.pdb_hierarchys  = []
    self.er_data.ke_pdb          = []
    for x in xrange(len(self.fmodel_total_block_list)):
      if x >= best_r_work_block[0] and x < best_r_work_block[1]:
        print  >> self.log, "Block | Number of models in block : ", x+1, " | ", len(self.fmodel_total_block_list[x][2])
        self.er_data.xray_structures.extend(self.fmodel_total_block_list[x][2])
        self.er_data.xray_structures_diff_map.extend(self.fmodel_total_block_list[x][3])
        self.er_data.pdb_hierarchys.extend(self.fmodel_total_block_list[x][4])
        self.er_data.ke_pdb.extend(self.fmodel_total_block_list[x][5])
    assert len(self.er_data.xray_structures) == len(self.er_data.pdb_hierarchys)
    assert len(self.er_data.xray_structures) == len(self.er_data.ke_pdb)
    print >> self.log, "Number of models for PBD          : ", len(self.er_data.xray_structures)
    print >> self.log, "|"+"-"*77+"|\n"

  def print_fmodels_scale_and_solvent_stats(self):
    make_header("Fmodel statistics | macrocycle: "+str(self.macro_cycle),
      out = self.log)
    print >> self.log, '{0:<23}: {1:>8} {2:>8} {3:>8} {4:>8}'.format('','MC',
      'k1','Bsol','ksol')
    if self.fmodel_current is not None:
      print >> self.log, "{0:<23}: {1:8d} {2:8.3f} {3:8.3f} {4:8.3f}"\
        .format('Fmodel current',
                self.macro_cycle,
                self.fmodel_current.scale_k1(),
                self.fmodel_current.fmodel_kbu().b_sol(),
                self.fmodel_current.fmodel_kbu().k_sols()[0],
                )
    if self.fmodel_running is not None:
      print >> self.log, "{0:<23}: {1:8d} {2:8.3f} {3:8.3f} {4:8.3f}"\
        .format('Fmodel running',
                self.macro_cycle,
                self.fmodel_running.scale_k1(),
                self.fmodel_running.fmodel_kbu().b_sol(),
                self.fmodel_running.fmodel_kbu().k_sols()[0] )
    if self.fmodel_total is not None:
      print >> self.log, "{0:<23}: {1:8d} {2:8.3f} {3:8.3f} {4:8.3f}"\
        .format('Fmodel_Total',
                self.macro_cycle,
                self.fmodel_total.scale_k1(),
                self.fmodel_total.fmodel_kbu().b_sol(),
                self.fmodel_total.fmodel_kbu().k_sols()[0] )
    if self.fmodel_current is not None:
      print >> self.log, "Fmodel current bcart   : {0:14.2f} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f} {5:5.2f}".format(*self.fmodel_current.fmodel_kbu().b_cart())
    if self.fmodel_running is not None:
      print >> self.log, "Fmodel running bcart   : {0:14.2f} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f} {5:5.2f}".format(*self.fmodel_running.fmodel_kbu().b_cart())
    if self.fmodel_total  is not None:
      print >> self.log, "Fmodel total bcart     : {0:14.2f} {1:5.2f} {2:5.2f} {3:5.2f} {4:5.2f} {5:5.2f}".format(*self.fmodel_total.fmodel_kbu().b_cart())
    print >> self.log, "|"+"-"*77+"|\n"

  def write_diff_map_ensemble(self, out):
    crystal_symmetry = self.er_data.xray_structures_diff_map[0].crystal_symmetry()
    print >> out, pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
    print >> out, pdb.format_scale_records(unit_cell = crystal_symmetry.unit_cell())
    for n,xrs in enumerate(self.er_data.xray_structures_diff_map):
      print >> out, "MODEL %8d"%(n+1)
      print >> out, xrs.as_pdb_file()
      print >> out, "ENDMDL"
    print >> out, "END"

  def update_single_hierarchy (self, i_model) :
    xrs = self.er_data.xray_structures[i_model]
    scatterers = xrs.scatterers()
    sites_cart = xrs.sites_cart()
    u_isos = xrs.extract_u_iso_or_u_equiv()
    occupancies = scatterers.extract_occupancies()
    u_carts = scatterers.extract_u_cart_plus_u_iso(xrs.unit_cell())
    scat_types = scatterers.extract_scattering_types()
    i_model_pdb_hierarchy = self.er_data.pdb_hierarchys[i_model]
    pdb_atoms = i_model_pdb_hierarchy.atoms()
    i_model_ke = self.er_data.ke_pdb[i_model]
    for j_seq, atom in enumerate(pdb_atoms):
      if j_seq < len(sites_cart):
        atom.xyz = sites_cart[j_seq]
        if self.params.output_running_kinetic_energy_in_occupancy_column:
          #XXX * 0.1 to fit in occ col
          atom.occ = 0.1 * i_model_ke[j_seq]
        else:
          atom.occ = 1.0 / len(self.er_data.xray_structures)
        atom.b = adptbx.u_as_b(u_isos[j_seq])
        e = scat_types[j_seq]
        if (len(e) > 1 and "+-0123456789".find(e[1]) >= 0):
          atom.element = "%2s" % e[:1]
          atom.charge = "%-2s" % e[1:]
        elif (len(e) > 2):
          atom.element = "%2s" % e[:2]
          atom.charge = "%-2s" % e[2:]
        else:
          atom.element = "%2s" % e
          atom.charge = "  "
        if (scatterers[j_seq].flags.use_u_aniso()):
          atom.uij = u_carts[j_seq]
        elif(False):
          atom.uij = self.u_cart
        else:
          atom.uij = (-1,-1,-1,-1,-1,-1)
    return i_model_pdb_hierarchy

  def write_ensemble_pdb(self, out):
    crystal_symmetry = self.er_data.xray_structures[0].crystal_symmetry()
    pr = "REMARK   3"
    print >> out, pr
    print >> out,  "REMARK   3 TIME-AVERAGED ENSEMBLE REFINEMENT."
    from phenix import phenix_info # FIXME ???
    ver, tag = phenix_info.version_and_release_tag(f = out)
    if(ver is None):
      prog = "   PROGRAM     : PHENIX (phenix.ensemble_refinement)"
    else:
      if(tag is not None):
        ver = ver+"_"+tag
      prog = "   PROGRAM     : PHENIX (phenix.ensemble_refinement: %s)"%ver
    print >> out,pr+prog
    authors = phenix_info.phenix_developers_last
    l = pr+"   AUTHORS     :"
    j = 0
    i = j
    n = len(l) + 1
    while (j != len(authors)):
      a = len(authors[j]) + 1
      if (n+a > 79):
        print >> out,l, ",".join(authors[i:j]) + ","
        l = pr+"               :"
        i = j
        n = len(l) + 1
      n += a
      j += 1
    if (i != j):
      print >> out,l, ",".join(authors[i:j])
    fmodel_info = self.fmodel_total.info()
    fmodel_info.show_remark_3(out = out)
#    model_stats = mmtbx.model_statistics.model(model     = self.model,
#                                               ignore_hd = False)
#    # set mode_stats.geometry to None as refers to final structure NOT ensemble
#    model_stats.geometry = None
#    model_stats.show(out = out, pdb_deposition =True)
    # get mean geometry stats for ensemble
    self.final_geometry_pdb_string = self.ensemble_utils.ensemble_mean_geometry_stats(
        restraints_manager       = self.model.restraints_manager,
        xray_structure           = self.model.xray_structure,
        ensemble_xray_structures = self.er_data.xray_structures,
        ignore_hd                = True,
        verbose                  = False,
        out                      = self.log,
        return_pdb_string        = True)
    print >> out, self.final_geometry_pdb_string
    print >> out, pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry)
    print >> out, pdb.format_scale_records(unit_cell = crystal_symmetry.unit_cell())
    atoms_reset_serial = True
    #
    cntr = 0
    assert len(self.er_data.ke_pdb) == len(self.er_data.xray_structures)
    assert len(self.er_data.pdb_hierarchys) == len(self.er_data.xray_structures)
    for i_model, xrs in enumerate(self.er_data.xray_structures):
      cntr += 1
      print >> out, "MODEL %8d"%cntr
      i_model_pdb_hierarchy = self.update_single_hierarchy(i_model)
      if (atoms_reset_serial):
        atoms_reset_serial_first_value = 1
      else:
        atoms_reset_serial_first_value = None
      out.write(i_model_pdb_hierarchy.as_pdb_string(
        append_end=False,
        atoms_reset_serial_first_value=atoms_reset_serial_first_value))
      #
      print >> out, "ENDMDL"
    print >> out, "END"

  def print_ml_stats(self):
    if self.fmodel_running.set_sigmaa is not None:
      self.run_time_stats_dict.update({'Sigma_a':self.fmodel_running.set_sigmaa})
    if self.params.target_name in ['ml', 'mlhl'] :
      self.run_time_stats_dict.update({'Alpha':self.fmodel_running.alpha_beta()[0].data()})
      self.run_time_stats_dict.update({'Beta':self.fmodel_running.alpha_beta()[1].data()})
    if self.fmodel_running.n_obs is not None:
      self.run_time_stats_dict.update({'Eobs(fixed)':self.fmodel_running.n_obs})
    if self.fmodel_running.n_calc is not None:
      self.run_time_stats_dict.update({'Ecalc(fixed)':self.fmodel_running.n_calc})

    make_header("ML statistics", out = self.log)
    print >> self.log, '  {0:<23}: {1:>12} {2:>12} {3:>12}'.format('','min','max','mean')
    for key in sorted(self.run_time_stats_dict.keys()):
      info = self.run_time_stats_dict[key].min_max_mean()
      print >> self.log, '  {0:<23}: {1:12.3f} {2:12.3f} {3:12.3f}'.format(
        key,
        info.min,
        info.max,
        info.mean)
    print >> self.log, "|"+"-"*77+"|\n"

  def wilson_plot(self, miller_data, header = ""):
    make_header("Wilson Plot " + header, out = self.log)
    self.copy_ma = self.copy_ma.array(data = flex.abs(miller_data))
    reflections_per_bin = min(250, self.copy_ma.data().size())
    self.copy_ma.setup_binner(reflections_per_bin = reflections_per_bin)
    wilson_plot = self.copy_ma.wilson_plot(use_binning=True)
    wilson_plot.show(f = self.log)

################################################################################

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels, log):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags().data().count(True)*1./fmodel.r_free_flags().data().size()
  print >> log, "Data statistics"
  try: f_obs_labels = f_obs_labels[:f_obs_labels.index(",")]
  except ValueError: pass
  result = " \n    ".join([
    "data_label                          : %s" % f_obs_labels,
    "high_resolution                     : "+format_value("%-5.2f",info.d_min),
    "low_resolution                      : "+format_value("%-6.2f",info.d_max),
    "completeness_in_range               : " + \
      format_value("%-6.2f",info.completeness_in_range),
    "completeness(d_min-inf)             : " + \
      format_value("%-6.2f",info.completeness_d_min_inf),
    "completeness(6A-inf)                : " + \
      format_value("%-6.2f",info.completeness_6_inf),
    "wilson_b                            : " + \
      format_value("%-6.1f",fmodel.wilson_b()),
    "number_of_reflections               : " + \
      format_value("%-8d",  info.number_of_reflections),
    "test_set_size                       : " + \
      format_value("%-8.4f",flags_pc),
    "test_flag_value                     : " + \
      format_value("%-d",   test_flag_value),
    "number_of_Fobs_outliers             : " + format_value("%-8d",  n_outl),
    "anomalous_flag                      : " + \
      format_value("%-6s",  fmodel.f_obs().anomalous_flag())])
  print >> log, "   ", result

def show_model_vs_data(fmodel, log):
  d_max, d_min = fmodel.f_obs().d_max_min()
  flags_pc = fmodel.r_free_flags().data().count(True)*100./\
    fmodel.r_free_flags().data().size()
  if(flags_pc == 0): r_free = None
  else: r_free = fmodel.r_free()
  k_sol = format_value("%-5.2f",fmodel.fmodel_kbu().k_sols()[0])
  b_sol = format_value("%-7.2f",fmodel.fmodel_kbu().b_sol())
  b_cart = " ".join([("%8.2f"%v).strip() for v in fmodel.fmodel_kbu().b_cart()])
  print >> log, "Model vs data statistics"
  result = " \n    ".join([
    "r_work(re-computed)                 : " + \
      format_value("%-6.4f",fmodel.r_work()),
    "r_free(re-computed)                 : " + format_value("%-6.4f",r_free),
    "scale_k1                            : " + \
      format_value("%-6.4f",fmodel.scale_k1()),
    "bulk_solvent_(k_sol,b_sol)          : %s%s" % (k_sol, b_sol),
    "overall_anisotropic_scale_(b_cart)  : " + format_value("%-s",b_cart)])
  print >> log, "   ", result

def write_mtz_file (fmodel_total, raw_data, raw_flags, prefix, params) :
  assert (fmodel_total is not None)
  class labels_decorator:
    def __init__(self, amplitudes_label, phases_label):
      self._amplitudes = amplitudes_label
      self._phases = phases_label
    def amplitudes(self):
      return self._amplitudes
    def phases(self, root_label, anomalous_sign=None):
      assert anomalous_sign is None or not anomalous_sign
      return self._phases
  xray_suffix = "_xray"
  f_obs_label = "F-obs"
  i_obs_label = "I-obs"
  flag_label = "R-free-flags"
  if (raw_data.is_xray_intensity_array()) :
    column_root_label = i_obs_label
  else:
    column_root_label = f_obs_label
  mtz_dataset_original = raw_data.as_mtz_dataset(
    column_root_label=column_root_label)
  mtz_dataset_original.add_miller_array(
    miller_array = raw_flags,
    column_root_label=flag_label)
  mtz_dataset_original.set_name("Original-experimental-data")
  new_dataset = mtz_dataset_original.mtz_crystal().add_dataset(
    name = "Experimental-data-used-in-refinement", wavelength=1)
  new_dataset.add_miller_array(
    miller_array = fmodel_total.f_obs(),
    column_root_label="F-obs-filtered"+xray_suffix)
  another_dataset = new_dataset.mtz_crystal().add_dataset(
    name = "Model-structure-factors-(bulk-solvent-and-all-scales-included)",
    wavelength=1)
  another_dataset.add_miller_array(
    miller_array = fmodel_total.f_model_scaled_with_k1_composite_work_free(),
    column_root_label="F-model"+xray_suffix)
  yet_another_dataset = another_dataset.mtz_crystal().add_dataset(
    name = "Fourier-map-coefficients", wavelength=1)
  cmo = mmtbx.maps.compute_map_coefficients(
    fmodel = fmodel_total,
    params = params.electron_density_maps.map_coefficients)
  for ma in cmo.mtz_dataset.mtz_object().as_miller_arrays():
    labels=ma.info().labels
    ld = labels_decorator(amplitudes_label=labels[0], phases_label=labels[1])
    yet_another_dataset.add_miller_array(
      miller_array      = ma,
      column_root_label = labels[0],
      label_decorator   = ld)
  yet_another_dataset.mtz_object().write(
    file_name=prefix+".mtz")
  return prefix + ".mtz"

#-----------------------------------------------------------------------
def run(args, command_name = "phenix.ensemble_refinement", out=None,
    validate=False, replace_stderr=True):
  if(len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="%s reflection_file pdb_file [options]" % command_name,
    description='Example: %s data.mtz model.pdb' % command_name
  ).enable_dry_run().enable_show_defaults()).process(args=args)
  if (out is None) :
    out = sys.stdout
  if(command_line.expert_level is not None) :
    master_params.show(
      expert_level=command_line.expert_level,
      attributes_level=command_line.attributes_level,
      out=out)
    return
  inputs = mmtbx.command_line.load_model_and_data(
    args=command_line.args,
    master_phil=master_params,
    out=out,
    create_fmodel=False,
    process_pdb_file=False)
  working_phil = inputs.working_phil
  params = working_phil.extract()
  if (params.extra_restraints_file is not None) :
    # XXX this is a revolting hack...
    print >> out, "Processing custom geometry restraints in file:"
    print >> out, "  %s" % params.extra_restraints_file
    restraints_phil = iotbx.phil.parse(file_name=params.extra_restraints_file)
    cleanup_phil = iotbx.phil.parse("extra_restraints_file=None")
    working_phil = master_params.fetch(
      sources=[working_phil, restraints_phil, cleanup_phil])
    params = working_phil.extract()
  er_params = params.ensemble_refinement

  if er_params.electron_density_maps.apply_default_maps != False\
    or len(er_params.electron_density_maps.map_coefficients) == 0:
    maps_par = libtbx.env.find_in_repositories(
      relative_path=\
        "cctbx_project/mmtbx/refinement/ensemble_refinement/maps.params",
      test=os.path.isfile)
    maps_par_phil = iotbx.phil.parse(file_name=maps_par)
    working_params = mmtbx.refinement.ensemble_refinement.master_params.fetch(
                        sources = [working_phil]+[maps_par_phil])
    er_params = working_params.extract().ensemble_refinement

  if er_params.output_file_prefix == None:
    er_params.output_file_prefix = os.path.splitext(
      os.path.basename(inputs.pdb_file_names[0]))[0] + "_ensemble"
  log = multi_out()
  log.register(label="stdout", file_object=out)
  log.register(
    label="log_buffer",
    file_object=StringIO(),
    atexit_send_to=None)
  if (replace_stderr) :
    sys.stderr = log
  log_file = open(er_params.output_file_prefix+'.log', "w")
  log.replace_stringio(
      old_label="log_buffer",
      new_label="log",
      new_file_object=log_file)
  timer = user_plus_sys_time()
  mmtbx.utils.print_programs_start_header(log=log, text=command_name)
  make_header("Ensemble refinement parameters", out = log)
  working_phil.show(out = log)
  make_header("Model and data statistics", out = log)
  print >> log, "Data file                               : %s" % \
    format_value("%5s", os.path.basename(params.input.xray_data.file_name))
  print >> log, "Model file                              : %s \n" % \
    (format_value("%5s",os.path.basename(inputs.pdb_file_names[0])))
  print >> log, "\nTLS MUST BE IN ATOM RECORDS OF INPUT PDB\n"
  f_obs = inputs.f_obs
  number_of_reflections = f_obs.indices().size()

  r_free_flags = inputs.r_free_flags
  raw_flags = inputs.raw_flags
  raw_data = inputs.raw_data

  print >> log, "\nPDB file name : ", inputs.pdb_file_names[0]

  # Process PDB file
  cif_objects = inputs.cif_objects
  pdb_file = inputs.pdb_file_names[0]
  pdb_ip = mmtbx.monomer_library.pdb_interpretation.master_params.extract()
  pdb_ip.clash_guard.nonbonded_distance_threshold = -1.0
  pdb_ip.clash_guard.max_number_of_distances_below_threshold = 100000000
  pdb_ip.clash_guard.max_fraction_of_distances_below_threshold = 1.0
  pdb_ip.proceed_with_excessive_length_bonds=True
  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
    cif_objects               = cif_objects,
    pdb_interpretation_params = pdb_ip,
    crystal_symmetry          = inputs.xray_structure,
    log                       = log)
  processed_pdb_file, pdb_inp = \
    processed_pdb_files_srv.process_pdb_files(pdb_file_names = [pdb_file])

  # Remove alternative conformations if present
  hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  if er_params.remove_alt_conf_from_input_pdb:
    atoms_size_pre = hierarchy.atoms().size()
    for model in hierarchy.models() :
      for chain in model.chains() :
        for residue_group in chain.residue_groups() :
          atom_groups = residue_group.atom_groups()
          assert (len(atom_groups) > 0)
          for atom_group in atom_groups :
            if (not atom_group.altloc in ["", "A"]) :
              residue_group.remove_atom_group(atom_group=atom_group)
            else :
              atom_group.altloc = ""
          if (len(residue_group.atom_groups()) == 0) :
            chain.remove_residue_group(residue_group=residue_group)
        if (len(chain.residue_groups()) == 0) :
          model.remove_chain(chain=chain)
    atoms = hierarchy.atoms()
    new_occ = flex.double(atoms.size(), 1.0)
    atoms.set_occ(new_occ)
    atoms_size_post = hierarchy.atoms().size()
    if atoms_size_pre != atoms_size_post:
      pdb_file_removed_alt_confs = pdb_file[0:-4]+'_removed_alt_confs.pdb'
      print >> log, "\nRemoving alternative conformations"
      print >> log, "All occupancies reset to 1.0"
      print >> log, "New PDB : ", pdb_file_removed_alt_confs, "\n"
      hierarchy.write_pdb_file(file_name        = pdb_file_removed_alt_confs,
                               crystal_symmetry = pdb_inp.crystal_symmetry())
      processed_pdb_file, pdb_inp = \
      processed_pdb_files_srv.process_pdb_files(
        pdb_file_names = [pdb_file_removed_alt_confs])

  d_min = f_obs.d_min()
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = "n_gaussian",
    d_min               = d_min,
    log                = log)
  if(len(xsfppf.xray_structures) > 1):
    raise Sorry("Multiple models not supported.")
  xray_structure = xsfppf.xray_structures[0].deep_copy_scatterers()

  # TODO Amber hooks
  amber_structs = use_amber = None
  if hasattr(params.ensemble_refinement, "amber") :
    use_amber = params.ensemble_refinement.amber.use_amber
    if (use_amber) :
      amber_params = params.ensemble_refinement.amber
      import amber_adaptbx
      make_header("Initializing AMBER", out=log)
      print >> log, "  topology: %s" % amber_params.topology_file_name
      amber_structs = amber_adaptbx.get_amber_structs(
        prmtop=amber_params.topology_file_name,
        ambcrd=amber_params.coordinate_file_name)
      if (params.hydrogens.refine.lower() == "auto") :
      #  params.hydrogens.refine = "individual"
        params.hydrogens.force_riding_adp = True

  # Geometry manager
  sctr_keys = \
         xray_structure.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys

  geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      plain_pairs_radius           = 5,
      params_edits                 = \
        params.refinement.geometry_restraints.edits,
      params_remove                = None,
      hydrogen_bond_proxies        = None,
      hydrogen_bond_params         = None,
      custom_nonbonded_exclusions  = None,
      external_energy_function     = None,
      assume_hydrogens_all_missing = not has_hd)

  restraints_manager = mmtbx.restraints.manager(
      geometry      = geometry,
      normalization = True,
      use_amber     = use_amber,
      amber_structs = amber_structs)

  # Refinement flags
  class rf:
    def __init__(self, size):
      self.individual_sites     = True
      self.individual_adp       = False
      self.sites_individual     = flex.bool(size, True)
      self.sites_torsion_angles = None
      self.torsion_angles       = None
      self.adp_individual_iso   = None
      self.adp_individual_aniso = None
    def inflate(self, **keywords): pass

  refinement_flags = rf(size = xray_structure.scatterers().size())

  # Model
  model = mmtbx.model.manager(
    processed_pdb_files_srv = processed_pdb_files_srv,
    refinement_flags = refinement_flags,
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = hierarchy,
    tls_groups = None,
    anomalous_scatterer_groups = None,
    log = log)

  # Geometry file
  xray_structure = model.xray_structure
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  model.restraints_manager.geometry.show_sorted(
    sites_cart=sites_cart,
    site_labels=site_labels,
    f=open(er_params.output_file_prefix+'.geo','w') )

  print >> log, "Unit cell                               :", f_obs.unit_cell()
  print >> log, "Space group                             :", \
    f_obs.crystal_symmetry().space_group_info().symbol_and_number()
  print >> log, "Number of symmetry operators            :", \
    f_obs.crystal_symmetry().space_group_info().type().group().order_z()
  print >> log, "Unit cell volume                        : %-15.4f" % \
    f_obs.unit_cell().volume()
  f_obs_labels = f_obs.info().label_string()

  if (command_line.options.dry_run) :
    return None

  fmodel = mmtbx.utils.fmodel_simple(
    f_obs                      = f_obs,
    xray_structures            = [model.xray_structure],
    scattering_table           = "n_gaussian",
    r_free_flags               = r_free_flags,
    target_name                = er_params.target_name,
    bulk_solvent_and_scaling   = False,
    bss_params                 = None,
    mask_params                = None,
    twin_laws                  = None,
    skip_twin_detection        = True,
    twin_switch_tolerance      = 2.0,
    outliers_rejection         = True,
    bulk_solvent_correction    = True,
    anisotropic_scaling        = True,
    log                        = log)
  hl_coeffs = inputs.hl_coeffs
  if (hl_coeffs is not None) and (params.input.use_experimental_phases) :
    print >> log, "Using MLHL target with experimental phases"
    er_params.target_name = "mlhl"
    hl_coeffs = hl_coeffs.common_set(other=fmodel.f_obs())
  else :
    hl_coeffs = None
  # XXX is this intentional?
  fmodel = mmtbx.f_model.manager(
    mask_params                  = er_params.mask,
    xray_structure               = model.xray_structure,
    f_obs                        = fmodel.f_obs(),
    r_free_flags                 = fmodel.r_free_flags(),
    target_name                  = er_params.target_name,
    abcd                         = hl_coeffs)
  hd_sel = model.xray_structure.hd_selection()
  model.xray_structure.set_occupancies(
        value     = 1.0,
        selection = hd_sel)
  model.show_occupancy_statistics(out = log)

  fmodel.update_xray_structure(
    xray_structure      = model.xray_structure,
    update_f_calc       = True,
    update_f_mask       = False,
    force_update_f_mask = False)

  n_outl = f_obs.data().size() - fmodel.f_obs().data().size()
  show_data(fmodel          = fmodel,
            n_outl          = n_outl,
            test_flag_value = inputs.test_flag_value,
            f_obs_labels    = f_obs_labels,
            log             = log)
  show_model_vs_data(fmodel = fmodel,
                     log    = log)

  best_trial = None
  if (len(er_params.ptls) == 1) :
    best_trial = run_wrapper(
      fmodel               = fmodel,
      model                = model,
      er_params            = er_params,
      raw_data             = raw_data,
      raw_flags            = raw_flags,
      log                  = log).__call__(
        ptls=er_params.ptls[0],
        buffer_output=False,
        append_ptls=False)
  else :
    driver = run_wrapper(
      fmodel               = fmodel,
      model                = model,
      er_params            = er_params,
      raw_data             = raw_data,
      raw_flags            = raw_flags,
      log                  = log)
    trials = []
    if (er_params.nproc in [1, None]) or (sys.platform == "win32") :
      for ptls in er_params.ptls :
        make_header("Running with pTLS = %g" % ptls, out=log)
        trial_result = driver(ptls, buffer_output=False, write_log=False)
        assert (trial_result is not None)
        trials.append(trial_result)
    else :
      trials = easy_mp.pool_map(
        fixed_func=driver,
        args=er_params.ptls,
        processes=er_params.nproc)
    assert (not None in trials)
    best_trial = min(trials, key=lambda t: t.r_free)
    best_trial.save_final(er_params.output_file_prefix)

  show_total_time(out = log)
  return result(
    best_trial=best_trial,
    prefix=er_params.output_file_prefix,
    validate=validate,
    log=log)

class run_wrapper (object) :
  def __init__ (self, model, fmodel, raw_data, raw_flags, er_params, log) :
    adopt_init_args(self, locals())

  def __call__ (self, ptls, buffer_output=True, write_log=True,
      append_ptls=True) :
    out = self.log
    log_out = None
    if (buffer_output) :
      out = StringIO()
    run_number = None
    if (append_ptls) :
      run_number = ptls
    ensemble_refinement = run_ensemble_refinement(
      fmodel               = self.fmodel.deep_copy(),
      model                = self.model.deep_copy(),
      params               = self.er_params,
      raw_data             = self.raw_data,
      raw_flags            = self.raw_flags,
      run_number           = run_number,
      ptls                 = ptls,
      log                  = out)
    if (buffer_output) :
      log_out = out.getvalue()
      if (write_log):
        log_file_name = self.er_params.output_file_prefix + '_ptls-' + \
          str(ptls) + '.log'
        log_file = open(log_file_name, 'w')
        log_file.write(log_out)
    return trial(
      ptls=ptls,
      r_work=ensemble_refinement.fmodel_total.r_work(),
      r_free=ensemble_refinement.fmodel_total.r_free(),
      pdb_file=ensemble_refinement.pdb_file,
      mtz_file=ensemble_refinement.mtz_file,
      log_out=log_out,
      number_of_models=len(ensemble_refinement.er_data.xray_structures))

class trial (slots_getstate_setstate) :
  __slots__ = ["r_work", "r_free", "pdb_file", "mtz_file", "number_of_models",
               "log_out", "ptls"]
  def __init__ (self, **kwds) :
    kwds = dict(kwds)
    for name in self.__slots__ :
      setattr(self, name, kwds[name])

  def save_final (self, prefix) :
    pdb_out = prefix + ".pdb"
    if (self.pdb_file.endswith(".gz")) :
      pdb_out += ".gz"
    os.rename(self.pdb_file, pdb_out)
    os.rename(self.mtz_file, prefix + ".mtz")
    self.pdb_file = pdb_out
    self.mtz_file = prefix + ".mtz"

########################################################################
# Phenix GUI hooks
class result (slots_getstate_setstate) :
  __slots__ = [
    "directory", "r_work", "r_free",
    "number_of_models", "pdb_file", "mtz_file","validation",
  ]
  def __init__ (self,
      best_trial,
      prefix,
      log,
      validate=False) :
    for attr in ["r_work", "r_free", "number_of_models", "pdb_file",
                 "mtz_file"] :
      setattr(self, attr, getattr(best_trial, attr))
    self.directory = os.getcwd()
    self.validation = None
    if (validate) :
      from mmtbx.command_line import validation_summary
      self.validation = validation_summary.run(
        args=[self.pdb_file],
        out=log)
      assert (type(self.validation).__name__ == 'ensemble')

  def get_result_files (self, output_dir=None) :
    if (output_dir is None) :
      output_dir = self.directory
    return (os.path.join(self.directory, self.pdb_file),
            os.path.join(self.directory, self.mtz_file))

  def finish_job (self) :
    pdb_file, mtz_file = self.get_result_files()
    return (
      [(pdb_file, "Final ensemble"),
       (mtz_file, "Map coefficients")],
      [("R-work", "%.4f" % self.r_work),
       ("R-free", "%.4f" % self.r_free),
       ("Models", str(self.number_of_models)),]
    )

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    os.mkdir(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args),
      out=sys.stdout,
      validate=True)

def validate_params (params) :
  if (params.ensemble_refinement.ptls is None) :
    raise Sorry("You must specify a fraction of atoms to use for TLS fitting.")
  elif (len(params.input.xray_data.labels[0].split(",")) > 2) :
    raise Sorry("Anomalous data are not allowed in this program.")
  return mmtbx.command_line.validate_input_params(params)
