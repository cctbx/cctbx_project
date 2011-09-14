from libtbx import adopt_init_args
from cctbx.array_family import flex
import iotbx.phil
import mmtbx.utils
from mmtbx import find_peaks
from cctbx import xray
import random
from mmtbx.dynamics.constants import boltzmann_constant_akma

master_params_str = """\
  tolerance = 1.0
    .type = float
  reset_all = False
    .type = bool
    .help = Removes all water atoms prior to re-picking using mFobs-DFmodel and 2mFo-DFmodel
  output_residue_name = HOH
    .type=str
    .input_size = 50
  output_chain_id = S
    .type=str
    .input_size = 50
  output_atom_name = O
    .type=str
    .input_size = 50
  scattering_type = O
    .type=str
    .help = Defines scattering factors for newly added waters
    .expert_level=2
    .input_size = 50
  primary_map_type = mFo-DFmodel
    .type=str
  primary_map_cutoff = 3.0
    .type=float
  secondary_map_type = 2mFo-DFmodel
    .type=str
  secondary_map_cutoff = 3.0
    .type=float
  h_bond_min_mac = 1.8
    .type = float
    .short_caption = H-bond minimum for solvent-model
    .expert_level = 1
  h_bond_min_sol = 1.8
    .type = float
    .short_caption = H-bond minimum for solvent-solvent
    .expert_level = 1
  h_bond_max = 3.2
    .type = float
    .short_caption = Maximum H-bond length
    .expert_level = 1
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
    .short_caption = New solvent ADP type
    .expert_level = 1
  b_iso_min = 0.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
    .short_caption = Minimum B-factor
    .expert_level = 1
  b_iso_max = 100.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
    .short_caption = Maximum B-factor
    .expert_level = 1
  anisotropy_min = 0.1
    .type = float
    .help = For solvent refined as anisotropic: remove is less than this value
    .short_caption = Minimum anisotropic B-factor
    .expert_level = 1
  b_iso = None
    .type=float
    .help = Initial B-factor value for newly added water
    .short_caption = Initial B-factor value
    .expert_level = 1
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
    .expert_level = 1
  occupancy_min = 0.1
    .type=float
    .help = Minimum occupancy value, waters with smaller value will be rejected
    .short_caption = Minimum occupancy
  occupancy_max = 1.0
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
    .short_caption = Maximum occupancy
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
    .short_caption = Initial occupancy value
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
  preserve_existing_solvent = True
    .type = bool

  find_peaks {
    use_sigma_scaled_maps = True
      .type=bool
      .help = Default is sigma scaled map, map in absolute scale is used \
              otherwise.
    resolution_factor = 1./4.
      .type=float
    map_next_to_model
      .expert_level=2
    {
      min_model_peak_dist = 1.8
        .type=float
      max_model_peak_dist = 6.0
        .type=float
      min_peak_peak_dist = 1.8
        .type=float
      use_hydrogens = True
        .type = bool
    }
    max_number_of_peaks = None
      .type=int
      .expert_level=1
    peak_search
      .expert_level=1
    {
      peak_search_level = 1
        .type=int
      max_peaks = 0
        .type=int
        .short_caption=Maximum peaks
      interpolate = True
        .type=bool
      min_distance_sym_equiv = None
        .type=float
        .short_caption=Minimum distance between symmetry-equivalent atoms
      general_positions_only = False
        .type=bool
      min_cross_distance = 1.8
        .type=float
        .short_caption=Minimum cross distance
      min_cubicle_edge = 5
        .type=float
        .short_caption=Minimum edge length of cubicles used for \
          fast neighbor search
        .expert_level=2
    }
  }
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

class manager(object):
  def __init__(self, fmodel,
                     model,
                     params = None,
                     velocities = None,
                     log = None):
    ###### Methodology
    # 1) copy exisiting solvent x,y,z,v
    # 2) find peaks in fo-fc map, map back to exisiting solvent x,y,z
    #       if distance peak:solvent < 1.0Ang keep x,y,z,v
    # 3) Repick fo-fc map with preserved waters in place, generate new x,y,z,v
    #       #N.B. repicking is done to prevent clashes
    # 4) Assert all solvent atoms have peak < x sigma in 2fo-fc map
    adopt_init_args(self, locals())
    mmtbx.utils.print_header("Time-averaging ordered solvent update", out = self.log)
    if(self.params is None): self.params = master_params().extract()
    self.fpp = self.params.find_peaks
    assert self.fmodel.xray_structure is self.model.xray_structure
    self.move_solvent_to_the_end_of_atom_list()
    self.show(message = 'Number of waters in current model')
    if(not self.is_water_last()):
      raise RuntimeError("Water picking failed: solvent must be last.")

    #copy existing solvent information for later comparison
    self.existing_solvent_xrs_selection = self.model.solvent_selection()
    self.existing_solvent_xrs = self.model.xray_structure.select(self.existing_solvent_xrs_selection)
    if self.velocities is not None:
      self.existing_solvent_velocities = self.velocities.select(self.existing_solvent_xrs_selection)
      self.velocities = self.velocities.select(~self.existing_solvent_xrs_selection)
    self.remove_solvent()
    assert self.fmodel.xray_structure is self.model.xray_structure
    #
    assert self.params.primary_map_type is not None
    ##########################################################################
    #New method where prexisting solvent atoms near picked peaks are kept
    print >> self.log, "\nCycle 1 - Evaluate Existing Solvent Atoms vs Primary Map and Secondary Map"
    #Calculate if existing atoms has a significant peak within 1.0A of either map
    #N.B. consider perfectly modelled water... - no diff map peak, significant 2Fo-Fc peak
    self.solvent_atoms_near_pick_selection = flex.bool(self.existing_solvent_xrs.scatterers().size(),False)
    #Solvent Pick from 1st Map
    self.sites = None
    self.heights = None

    #override distances for finding peaks
    store_min_peak_dist = self.fpp.map_next_to_model.min_model_peak_dist
    store_max_peak_dist = self.fpp.map_next_to_model.max_model_peak_dist
    self.fpp.map_next_to_model.min_model_peak_dist = 1.8
    self.fpp.map_next_to_model.max_model_peak_dist = 10.0

    peaks = self.find_peaks(
      map_type     = self.params.primary_map_type,
      map_cutoff   = self.params.primary_map_cutoff,
      use_kick_map = False).peaks_mapped()
    if(peaks is not None):
      self.sites, self.heights = peaks.sites, peaks.heights
    #Solvent Pick from 2nd Map (optional)
    if self.params.secondary_map_type is not None:
      peaks = self.find_peaks(
        map_type     = self.params.secondary_map_type,
        map_cutoff   = self.params.secondary_map_cutoff,
        use_kick_map = False).peaks_mapped()
      if(peaks is not None and self.sites is not None):
        self.sites.extend(peaks.sites)
        self.heights.extend(peaks.heights)
      elif (peaks is not None and self.sites is None):
        self.sites, self.heights = peaks.sites, peaks.heights
    if self.sites is not None:
      self.compare_peaks_with_positions()
    #Select solvent atoms with significant neighbouring peak (within 1.0A)...
    solvent_near_peak = self.existing_solvent_xrs.select(self.solvent_atoms_near_pick_selection)
    #Return solvent atoms near peaks
    self.sites = solvent_near_peak.sites_frac()
    self.add_new_solvent()
    #Book keeping for velocities (keep record of preserved waters)
    if self.velocities is not None:
      solvent_velocities_near_peaks = self.existing_solvent_velocities.select(self.solvent_atoms_near_pick_selection)
      self.velocities.extend(solvent_velocities_near_peaks)
    self.show(message = 'Number of preserved waters')
    assert self.fmodel.xray_structure is self.model.xray_structure
    ###########################################################################

    if(not self.is_water_last()):
      raise RuntimeError("Water picking failed: solvent must be last.")
    print >> self.log, "\nCycle 2 - Picking New Solvent Atoms from Fo-Fc Map and 2Fo-Fc Map"
    # Peak present in Fo-Fc and 2Fo-Fc map
    self.sites = None
    self.heights = None
    self.fpp.map_next_to_model.min_model_peak_dist = store_min_peak_dist
    self.fpp.map_next_to_model.max_model_peak_dist = store_max_peak_dist
    peaks_fo_fc = self.find_peaks(
      map_type     = self.params.primary_map_type,
      map_cutoff   = self.params.primary_map_cutoff,
      use_kick_map = False).peaks_mapped()
    if peaks_fo_fc.sites is not None:
      new_scatterers = flex.xray_scatterer(
        peaks_fo_fc.sites.size(),
        xray.scatterer(occupancy       = self.params.occupancy,
        b                              = self.params.b_iso,
        scattering_type                = self.params.scattering_type))
      new_scatterers.set_sites(peaks_fo_fc.sites)
      new_solvent_xray_structure = xray.structure(
        special_position_settings = self.model.xray_structure,
        scatterers                = new_scatterers)

      peaks_2fo_fc = self.find_peaks(
          map_type     = self.params.secondary_map_type,
          map_cutoff   = self.params.secondary_map_cutoff,
          use_kick_map = False).peaks_mapped()

      atom_nearest_to_peak = new_solvent_xray_structure.closest_distances(sites_frac = peaks_2fo_fc.sites, distance_cutoff = self.params.tolerance, use_selection = None)
      new_solvent_atoms_near_pick_selection = flex.bool(new_solvent_xray_structure.scatterers().size(),False)

      for x in atom_nearest_to_peak.i_seqs:
        if x > 0:
          new_solvent_atoms_near_pick_selection[x] = True
      print >> self.log, "Number of additional waters             : ", new_solvent_atoms_near_pick_selection.count(True)
      new_solvent_near_peak = new_solvent_xray_structure.select(new_solvent_atoms_near_pick_selection)
      #Return solvent atoms near peaks
      self.sites = new_solvent_near_peak.sites_frac()
      self.add_new_solvent()
      self.show(message = 'Total number of waters')
      assert self.fmodel.xray_structure is self.model.xray_structure

    if self.velocities is not None:
      self.model.xray_structure.scatterers().size()
      self.new_solvent_atom_velocities = flex.vec3_double((self.model.xray_structure.scatterers().size() - self.velocities.size()),[0,0,0])
      self.randomize_new_velocities()
      self.velocities.extend(self.new_solvent_atom_velocities)
    print >> self.log, "|"+"-"*77+"|\n"

  def randomize_new_velocities(self):
    if self.params.seed is not None:
      random.seed(self.params.seed)
    random_gauss = random.gauss
    random_random = random.random
    kt = boltzmann_constant_akma * self.params.temperature
    mass_oxygen = 16.0
    sigma = (kt / mass_oxygen)**0.5
    random_random()
    for x in xrange(self.new_solvent_atom_velocities.size()):
      self.new_solvent_atom_velocities[x] = [random_gauss(0, sigma) for i in (1,2,3)]

  def compare_peaks_with_positions(self):
    atom_nearest_to_peak = self.existing_solvent_xrs.closest_distances(sites_frac = self.sites, distance_cutoff = self.params.tolerance, use_selection = None)
    for x in atom_nearest_to_peak.i_seqs:
      if x > 0:
        self.solvent_atoms_near_pick_selection[x] = True

  def move_solvent_to_the_end_of_atom_list(self):
    solsel = flex.bool(self.model.solvent_selection().count(False), False)
    solsel.extend(flex.bool(self.model.solvent_selection().count(True), True))
    xrs_sol =  self.model.xray_structure.select(self.model.solvent_selection())
    if(xrs_sol.hd_selection().count(True) == 0):
      self.reset_solvent(
        solvent_selection      = solsel,
        solvent_xray_structure = xrs_sol)
    self.model.renumber_water()
    self.fmodel.xray_structure = self.model.xray_structure

  def remove_solvent(self):
    self.model = self.model.remove_solvent()
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = False)

  def reset_solvent(self, solvent_selection, solvent_xray_structure):
    assert solvent_selection.count(True) == \
      solvent_xray_structure.scatterers().size()
    self.model = self.model.remove_solvent()
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)

  def is_water_last(self):
    result = True
    sol_sel = self.model.solvent_selection()
    i_sol_sel = sol_sel.iselection()
    i_mac_sel = (~sol_sel).iselection()
    if(i_sol_sel.size() > 0 and i_mac_sel.size() > 0):
      if(flex.min_default(i_sol_sel,0)-flex.max_default(i_mac_sel,0) != 1):
        result = False
    return result

  def find_peaks(self, map_type, map_cutoff, use_kick_map=False):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = False)
    #N.B. essential that 'use_all_data    = False' so only working reflections are used
    return find_peaks.manager(fmodel          = self.fmodel,
                              map_type        = map_type,
                              map_cutoff      = map_cutoff,
                              params          = self.fpp,
                              use_kick_map    = False,
                              kick_map_params = None,
                              use_all_data    = False,
                              log             = self.log)

  def add_new_solvent(self):
    b_solv = self.params.b_iso
    # Isotropic
    new_scatterers = flex.xray_scatterer(
      self.sites.size(),
      xray.scatterer(occupancy       = self.params.occupancy,
      b                              = b_solv,
      scattering_type                = self.params.scattering_type))
    new_scatterers.set_sites(self.sites)
    solvent_xray_structure = xray.structure(
      special_position_settings = self.model.xray_structure,
      scatterers                = new_scatterers)
    xrs_sol = self.model.xray_structure.select(self.model.solvent_selection())
    xrs_mac = self.model.xray_structure.select(~self.model.solvent_selection())
    xrs_sol = xrs_sol.concatenate(other = solvent_xray_structure)
    sol_sel = flex.bool(xrs_mac.scatterers().size(), False)
    sol_sel.extend( flex.bool(xrs_sol.scatterers().size(), True) )
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = False)

  def show(self, message = None):
    if message is not None:
      spacer = " " * (40 - len(message) )
      print >> self.log, message+spacer+": ", self.model.number_of_ordered_solvent_molecules()
    else:
      print >> self.log, "Number of waters                        : ", self.model.number_of_ordered_solvent_molecules()
