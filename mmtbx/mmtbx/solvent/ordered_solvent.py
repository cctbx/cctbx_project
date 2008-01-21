from iotbx import pdb
from cctbx.array_family import flex
from cctbx import xray
import math,sys
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx import masks
from mmtbx import bulk_solvent
import mmtbx.restraints
from cctbx import crystal
from mmtbx.max_lik import max_like_non_uniform
from libtbx import adopt_init_args
import mmtbx.f_model
from libtbx import introspection
from cctbx import maptbx
from cctbx import adptbx
from cctbx import crystal
import libtbx.load_env
import iotbx.xplor.map
from scitbx import matrix
from libtbx.test_utils import approx_equal
import iotbx.phil
from libtbx.str_utils import format_value
from mmtbx import find_peaks

master_params = iotbx.phil.parse("""\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
  mode = *auto filter_only every_macro_cycle
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
  output_residue_name = HOH
    .type=str
  output_chain_id = None
    .type=str
  output_atom_name = O
    .type=str
  b_iso_min = 1.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
  b_iso_max = 60.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
  anisotropy_min = 0.1
    .type = float
    .help = For solvent refined as anisotropic: remove is less than this value
  b_iso = None
    .type=float
    .help = Initial B-factor value for newly added water
  scattering_type = O
    .type=str
    .help = Defines scattering factors for newly added waters
    .expert_level=2
  occupancy_min = 0.1
    .type=float
    .help = Minimum occupancy value, waters with smaller value will be rejected
  occupancy_max = 1.2
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
  primary_map_type = mFobs-DFmodel
    .type=str
  primary_map_cutoff = 3.0
    .type=float
  secondary_map_type = 2mFobs-DFmodel
    .type=str
  secondary_map_cutoff = 1.5
    .type=float
  use_h_bond_rejection_criteria = True
    .type = bool
  h_bond_min = 1.8
    .type = float
  h_bond_max = 3.2
    .type = float
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
""")

class water_ids(object):
  residue_names = ["HOH","SOL","SOLV","WAT","DOD","TIP3"]
  atom_names = ["O","OH2","H","H1","H2","D","D1","D2"]
  element_types = ["O","H","D","", " "]


class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     params = master_params.extract(),
                     find_peaks_params = None,
                     log    = None):
    adopt_init_args(self, locals())
    if(self.params is None): self.params = master_params.extract()
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    if(self.params.mode == "filter_only"): self.filter_only = True
    else: self.filter_only = False
    if(self.log is None): self.log = sys.stdout
    assert self.model.xray_structure == self.fmodel.xray_structure
    self.sites = None
    self.heights = None
    if(self.find_peaks_params.max_number_of_peaks is None):
      if(self.model.solvent_selection().count(False) > 0):
        self.find_peaks_params.max_number_of_peaks = \
          self.model.solvent_selection().count(False)
      else:
        self.find_peaks_params.max_number_of_peaks = \
          self.model.xray_structure.scatterers().size()
    self.move_solvent_to_the_end_of_atom_list()
    self.assert_water_is_last()
    self.show(message = "Start model:")
    self.filter_solvent()
    if(not self.filter_only):
      assert self.params.primary_map_type is not None
      peaks = self.find_peaks(
        map_type   = self.params.primary_map_type,
        map_cutoff = self.params.primary_map_cutoff).peaks_mapped()
      self.sites, self.heights = peaks.sites, peaks.heights
      self.add_new_solvent()
      self.filter_solvent()
      if(self.params.secondary_map_type is not None):
        self.find_peaks_2fofc()
        self.show(message = "2Fo-Fc map selection:")
    #
    if(not self.filter_only and self.params.refine_adp and
       self.model.refinement_flags.individual_adp and
       self.model.solvent_selection().count(True) > 0):
      self.fmodels.update_xray_structure(
         xray_structure = self.model.xray_structure,
         update_f_calc  = True,
         update_f_mask  = True)
      from mmtbx.refinement import minimization
      import scitbx.lbfgs
      if(self.params.new_solvent == "anisotropic"):
        selection_aniso = flex.bool(
          self.model.refinement_flags.adp_individual_aniso.size(), False)
        for sel in self.model.refinement_flags.adp_individual_aniso:
          selection_aniso = selection_aniso | sel
        selection_aniso.set_selected(~self.model.solvent_selection(), False)
        self.model.set_refine_individual_adp(selection_aniso = selection_aniso)
      else:
        selection_iso = self.model.refinement_flags.adp_individual_iso.deep_copy()
        selection_iso.set_selected(~self.model.solvent_selection(), False)
        assert selection_iso.count(True) == self.model.solvent_selection().count(True)
        self.model.set_refine_individual_adp(selection_iso = selection_iso)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = 25)
      minimized = minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = fmodels,
        model                    = model,
        refine_adp               = True,
        lbfgs_termination_params = lbfgs_termination_params)

  def move_solvent_to_the_end_of_atom_list(self):
    solsel = flex.bool(self.model.solvent_selection().count(False), False)
    solsel.extend(flex.bool(self.model.solvent_selection().count(True), True))
    xrs_sol =  self.model.xray_structure.select(self.model.solvent_selection())
    scat_types = xrs_sol.scatterers().extract_scattering_types()
    if(scat_types.count('H')+scat_types.count('D')==0):
      self.reset_solvent(
        solvent_selection      = solsel,
        solvent_xray_structure = xrs_sol)
    else: raise RuntimeError("Water contains H: not implemented.")

  def assert_water_is_last(self):
    sol_sel = self.model.solvent_selection()
    i_sol_sel = sol_sel.iselection()
    i_mac_sel = (~sol_sel).iselection()
    if(i_sol_sel.size() > 0 and i_mac_sel.size() > 0):
      if(not flex.max_default(i_sol_sel,0)-flex.min_default(i_mac_sel,0) != 1):
        raise RuntimeError("Water picking failed: solvent must be last.")

  def filter_solvent(self):
    sol_sel = self.model.solvent_selection()
    xrs_mac_h = self.model.xray_structure.select(~sol_sel)
    xrs_sol_h = self.model.xray_structure.select(sol_sel)
    hd_sol = self.model.xray_structure.hd_selection().select(sol_sel)
    hd_mac = self.model.xray_structure.hd_selection().select(~sol_sel)
    xrs_sol = xrs_sol_h.select(~hd_sol)
    xrs_mac = xrs_mac_h.select(~hd_mac)
    selection = xrs_sol.all_selection()
    scat_sol = xrs_sol.scatterers()
    occ_sol = scat_sol.extract_occupancies()
    b_isos_sol = scat_sol.extract_u_iso_or_u_equiv(
      self.model.xray_structure.unit_cell()) * math.pi**2*8
    anisotropy_sol = scat_sol.anisotropy(unit_cell =
      self.model.xray_structure.unit_cell())
    result = xrs_mac.closest_distances(sites_frac = xrs_sol.sites_frac(),
      distance_cutoff =
        self.find_peaks_params.map_next_to_model.max_model_peak_dist)
    selection &= b_isos_sol >= self.params.b_iso_min
    selection &= b_isos_sol <= self.params.b_iso_max
    selection &= occ_sol >= self.params.occupancy_min
    selection &= occ_sol <= self.params.occupancy_max
    selection &= anisotropy_sol > 0.1
    selection &= result.smallest_distances<= \
      self.find_peaks_params.map_next_to_model.max_model_peak_dist
    selection &= result.smallest_distances>= \
      self.find_peaks_params.map_next_to_model.min_model_peak_dist
    ########
    if(self.params.use_h_bond_rejection_criteria):
      aal = []
      for i_seq, sel in enumerate(self.model.xray_structure.hd_selection()):
        if(not sel):
          aal.append(self.model.atom_attributes_list[i_seq])
      sfs = xrs_sol.sites_frac()
      for i, i_seq in enumerate(result.i_seqs):
        if(selection[i]):
          if(result.smallest_distances[i] <= self.params.h_bond_max and
             result.smallest_distances[i] >= self.params.h_bond_min):
            assert aal[xrs_mac.scatterers().size()+i].element.strip() == 'O'
            if(aal[i_seq].element.strip() not in ['O','N']):
              selection[i] = False
          else:
            dist = 9999.
            for j, j_seq in enumerate(result.i_seqs):
              if(i != j):
                d = xrs_sol.unit_cell().distance(sfs[i], sfs[j])
                if(d < dist):
                  dist = d
            if(dist>= self.params.h_bond_max or dist<= self.params.h_bond_min):
              selection[i] = False
    new_ss = flex.bool(self.model.solvent_selection().count(False), True)
    new_ss.extend(selection)
    self.model = self.model.select(new_ss)

  def reset_solvent(self, solvent_selection, solvent_xray_structure):
    assert solvent_selection.count(True) == \
      solvent_xray_structure.scatterers().size()
    self.model = self.model.remove_solvent()
    self.model.add_solvent(
      solvent_selection      = solvent_selection,
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)

  def show(self, message):
    print >> self.log, message
    xrs_mac = self.model.xray_structure.select(~self.model.solvent_selection())
    xrs_sol = self.model.xray_structure.select(self.model.solvent_selection())
    scat = xrs_sol.scatterers()
    occ = scat.extract_occupancies()
    b_isos = scat.extract_u_iso_or_u_equiv(
      self.model.xray_structure.unit_cell()) * math.pi**2*8
    smallest_distances = xrs_mac.closest_distances(
      sites_frac      = xrs_sol.sites_frac(),
      distance_cutoff = self.find_peaks_params.map_next_to_model.\
        max_model_peak_dist).smallest_distances
    number = format_value("%-7d",scat.size())
    b_min  = format_value("%-7.2f", flex.min_default( b_isos, None))
    b_max  = format_value("%-7.2f", flex.max_default( b_isos, None))
    b_ave  = format_value("%-7.2f", flex.mean_default(b_isos, None))
    bl_min = format_value("%-7.2f", self.params.b_iso_min).strip()
    bl_max = format_value("%-7.2f", self.params.b_iso_max).strip()
    o_min  = format_value("%-7.2f", flex.min_default(occ, None))
    o_max  = format_value("%-7.2f", flex.max_default(occ, None))
    ol_min = format_value("%-7.2f", self.params.occupancy_min).strip()
    ol_max = format_value("%-7.2f", self.params.occupancy_max).strip()
    d_min  = format_value("%-7.2f", flex.min_default(smallest_distances, None))
    d_max  = format_value("%-7.2f", flex.max_default(smallest_distances, None))
    dl_min = format_value("%-7.2f",
      self.find_peaks_params.map_next_to_model.min_model_peak_dist).strip()
    dl_max = format_value("%-7.2f",
      self.find_peaks_params.map_next_to_model.max_model_peak_dist).strip()
    print >> self.log,"  number           = %s"%number
    print >> self.log,"  b_iso_min        = %s (limit = %s)"%(b_min, bl_min)
    print >> self.log,"  b_iso_max        = %s (limit = %s)"%(b_max, bl_max)
    print >> self.log,"  b_iso_mean       = %s             "%(b_ave)
    print >> self.log,"  occupancy_min    = %s (limit = %s)"%(o_min, ol_min)
    print >> self.log,"  occupancy_max    = %s (limit = %s)"%(o_max, ol_max)
    print >> self.log,"  dist_sol_mol_min = %s (limit = %s)"%(d_min, dl_min)
    print >> self.log,"  dist_sol_mol_max = %s (limit = %s)"%(d_max, dl_max)

  def find_peaks(self, map_type, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    return find_peaks.manager(fmodel     = self.fmodel,
                              map_type   = map_type,
                              map_cutoff = map_cutoff,
                              params     = self.find_peaks_params,
                              log        = self.log)

  def find_peaks_2fofc(self):
    peaks = self.find_peaks(
      map_type   = self.params.secondary_map_type,
      map_cutoff = self.params.secondary_map_cutoff).peaks()
    sites_2nd, heights_2nd = peaks.sites, peaks.heights
    step= self.fmodel.f_obs.d_min()*self.find_peaks_params.resolution_factor
    if(step < 0.3): step = 0.3 # XXX
    zz = self.model.xray_structure.select(self.model.solvent_selection())
    result = zz.closest_distances(sites_frac = sites_2nd,
      distance_cutoff = 6.0) # XXX
    smallest_distances = result.smallest_distances
    selection = (smallest_distances <= step) & (smallest_distances >= 0)
    cs = self.model.xray_structure.crystal_symmetry()
    sp = crystal.special_position_settings(cs)
    scatterers = flex.xray_scatterer()
    for site in result.sites_frac:
      scatterers.append(xray.scatterer("o", site))
    xrs_2nd = xray.structure(sp, scatterers)
    smallest_distances = xrs_2nd.closest_distances(sites_frac =
      zz.sites_frac(), distance_cutoff = 6.0).smallest_distances # XXX
    selection = (smallest_distances <= step) & (smallest_distances >= 0)
    new_ss = flex.bool(self.model.solvent_selection().count(False), True)
    new_ss.extend(selection)
    self.model = self.model.select(new_ss)

  def add_new_solvent(self):
    if(self.params.b_iso is None):
      b = self.model.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
      b_solv = flex.mean_default(b, None)
      if(b_solv is not None and b_solv < self.params.b_iso_min or
         b_solv > self.params.b_iso_max):
        b_solv = (self.params.b_iso_min + self.params.b_iso_max) / 2.
    else:
      b_solv = self.params.b_iso
    if(self.params.new_solvent == "isotropic"):
      new_scatterers = flex.xray_scatterer(
        self.sites.size(),
        xray.scatterer(occupancy       = self.params.occupancy,
                       b               = b_solv,
                       scattering_type = self.params.scattering_type))
    elif(self.params.new_solvent == "anisotropic"):
      u_star = adptbx.u_iso_as_u_star(self.model.xray_structure.unit_cell(),
        adptbx.b_as_u(b_solv))
      new_scatterers = flex.xray_scatterer(
        self.sites.size(),
        xray.scatterer(
          occupancy       = self.params.occupancy,
          u               = u_star,
          scattering_type = self.params.scattering_type))
    else: raise RuntimeError
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
      solvent_selection      = sol_sel,
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)

def show_histogram(data,
                   n_slots,
                   out=None,
                   prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix
    histogram = flex.histogram(data    = data,
                               n_slots = n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> out, "%7.3f - %7.3f: %d" % (low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    out.flush()
    return histogram
