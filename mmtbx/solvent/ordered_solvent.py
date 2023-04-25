from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import xray
import math,sys
from libtbx import adopt_init_args
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.phil
from mmtbx import find_peaks
from mmtbx.refinement import minimization
import scitbx.lbfgs
import mmtbx.utils
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from six.moves import range
from libtbx.utils import user_plus_sys_time
from libtbx.utils import null_out

output_params_str = """
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
"""

h_bond_params_str = """
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
"""

adp_occ_params_str = """
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
    .short_caption = New solvent ADP type
    .expert_level = 1
  b_iso_min = 1.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
    .short_caption = Minimum B-factor
    .expert_level = 1
  b_iso_max = 80.0
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
"""

master_params_str = """\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
    .short_caption = Minimum resolution
  mode = *second_half filter_only every_macro_cycle every_macro_cycle_after_first
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
    .short_caption = Mode
  n_cycles = 3
    .type = int
    .short_caption = Number of cycles
  %s
  primary_map_type = mFobs-DFmodel
    .type=str
    .help = Map used to identify candidate water sites - by default this is \
      the standard difference map.
  primary_map_cutoff = 3.0
    .type=float
    .short_caption = Primary map cutoff (sigma)
  secondary_map_and_map_cc_filter
    .short_caption = Secondary map filter
    .style = auto_align box
  {
    cc_map_1_type = "Fmodel"
      .type = str
      .short_caption = Model map type for CC calculation
    cc_map_2_type = 2mFo-DFmodel
      .type = str
      .short_caption = Experimental map type for CC calculation
    poor_cc_threshold = 0.7
      .type = float
      .short_caption = Minimum map correlation
    poor_map_value_threshold = 0.7
      .type = float
      .short_caption = Minimum map value (sigma)
  }
  %s
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
    .expert_level = 1
  %s
  filter_at_start = True
    .type = bool
    .expert_level = 1
  ignore_final_filtering_step = False
    .type = bool
    .expert_level=2
  correct_drifted_waters = True
    .type = bool
    .expert_level=2
""" % (output_params_str, h_bond_params_str, adp_occ_params_str)

# XXX
# XXX Consolidate with phenix.douse and may be more
# XXX
class msg_accumulator(object):
  def __init__(self, log=None):
    self.log = log
    self.messages = []
    self.prefix=""

  def set_prefix(self, prefix):
    self.prefix=prefix

  def add(self, msg):
    msg = "%s%s"%(self.prefix, msg)
    self.messages.append(msg)
    if(self.log is not None): print(msg, file=self.log)

  def show(self, log=None):
    assert [self.log, log].count(None) == 1
    if(log is not None):
      assert self.log is None
    else:
      assert self.log is not None
      log = self.log
    for msg in self.messages:
      print(msg, file=log)

def master_params():
  return iotbx.phil.parse(master_params_str)

class maps(object):
  def __init__(self, fmodel, map_1_type, map_2_type, grid_step=0.6, radius=2.0):
    self.radius = radius
    self.e_map = fmodel.electron_density_map()
    self.crystal_symmetry = fmodel.xray_structure.crystal_symmetry()
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = self.crystal_symmetry.unit_cell(),
      space_group_info = self.crystal_symmetry.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = grid_step)
    self.map_1 = self._get_real_map(map_type = map_1_type)
    self.map_2 = self._get_real_map(map_type = map_2_type)

  def _get_real_map(self, map_type):
    coeffs = self.e_map.map_coefficients(
      map_type     = map_type,
      fill_missing = False,
      isotropize   = True)
    fft_map = coeffs.fft_map(crystal_gridding = self.crystal_gridding)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

  def score_atom(self, atom, min_cc, min_value):
    site_cart = atom.xyz
    uc = self.crystal_symmetry.unit_cell()
    site_frac = uc.fractionalize(site_cart)
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = uc,
      fft_n_real = self.map_1.focus(),
      fft_m_real = self.map_1.all(),
      sites_cart = flex.vec3_double([site_cart]),
      site_radii = flex.double([self.radius]))
    cc = flex.linear_correlation(
      x=self.map_1.select(sel),
      y=self.map_2.select(sel)).coefficient()
    value_2 = self.map_2.eight_point_interpolation(site_frac)
    result = cc > min_cc and value_2 > min_value*atom.occ
    return result

class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     is_neutron_scat_table,
                     params = master_params().extract(),
                     find_peaks_params = None,
                     log = sys.stdout):
    adopt_init_args(self, locals())
    self.ma         = msg_accumulator(log = self.log)
    self.total_time = 0
    self._maps      = None
    self._peaks     = None
    self.n_water    = None
    self.model_size_init = self.model.size()
    #
    self._call(msg="Start")
    self._call(msg="Compute maps",     func=self._get_maps)
    self._call(msg="Filter",           func=self._filter_solvent)
    self._call(msg="Find peaks",       func=self._find_peaks)
    self._call(msg="Add new water",    func=self._add_new_solvent)
    self._call(msg="Refine new water", func=self._refine)
    self._call(msg="Compute maps",     func=self._get_maps)
    self._call(msg="Filter",           func=self._filter_solvent)
    self._call(msg="Correct drifted",  func=self._correct_drifted_waters)

  def _call(self, msg, func = None):
    timer = user_plus_sys_time()
    self.ma.add(msg)
    self._assert_same_model()
    if(func is not None): func()
    self._get_and_set_n_water()
    rs="r_work=%6.4f r_free=%6.4f"%(self.fmodel.r_work(), self.fmodel.r_free())
    nw="n_water=%3d"%(self.n_water)
    t = timer.elapsed()
    self.total_time += t
    self.ma.add("  %s | %s | time (s): %s (total time: %s)"%(rs, nw,
      ("%8.3f"%t).strip(), ("%8.3f"%self.total_time).strip()))

  def _get_and_set_n_water(self):
    self.n_water = self.model.solvent_selection().count(True)

  def _assert_same_model(self):
    #assert self.model.get_xray_structure() is self.fmodel.xray_structure
    mmtbx.utils.assert_xray_structures_equal(
      x1 = self.model.get_xray_structure(),
      x2 = self.fmodel.xray_structure)
    self.model.is_same_model(other=self.model)

  def _get_maps(self):
    self._maps = maps(
      fmodel     = self.fmodel,
      map_1_type = self.params.secondary_map_and_map_cc_filter.cc_map_1_type,
      map_2_type = self.params.secondary_map_and_map_cc_filter.cc_map_2_type)

  def _filter_solvent(self):
    sol_sel   = self.model.solvent_selection()
    hd_sel    = self.model.get_hd_selection()
    n_sol_start = self.n_water
    mfp = self.params.secondary_map_and_map_cc_filter
    # Select by distance
    distance_iselection = mmtbx.utils.select_water_by_distance(
      sites_frac_all      = self.model.get_xray_structure().sites_frac(),
      element_symbols_all = self.model.get_xray_structure().scattering_types(),
      water_selection_o   = sol_sel.set_selected(hd_sel, False).iselection(),
      dist_max            = self.params.h_bond_max,
      dist_min_mac        = self.params.h_bond_min_mac,
      dist_min_sol        = self.params.h_bond_min_sol,
      unit_cell           = self.model.get_xray_structure().unit_cell())
    distance_selection = flex.bool(self.model.size(), distance_iselection)
    selection = distance_selection.set_selected(~sol_sel, True)
    # Select by attributes
    get_class = iotbx.pdb.common_residue_names_get_class
    scatterers = self.model.get_xray_structure().scatterers()
    occ = scatterers.extract_occupancies()
    b_isos = scatterers.extract_u_iso_or_u_equiv(
      self.model.get_xray_structure().unit_cell()) * math.pi**2*8
    anisotropy = scatterers.anisotropy(unit_cell =
      self.model.get_xray_structure().unit_cell())
    for m in self.model.get_hierarchy().models():
      for c in m.chains():
        for conf in c.conformers():
          for r in conf.residues():
            if(not get_class(name=r.resname) == "common_water"): continue
            i_seqs = r.atoms().extract_i_seq()
            keep = True
            for atom in r.atoms():
              if(atom.element_is_hydrogen()): continue
              i_seq = atom.i_seq
              if(not distance_selection[i_seq]): keep = False
              if(atom.occ > self.params.occupancy_max or
                 atom.occ < self.params.occupancy_min): keep = False
              assert approx_equal(atom.occ, occ[i_seq], 1.e-3)
              if(anisotropy[i_seq] < self.params.anisotropy_min): keep = False
              if(b_isos[i_seq] < self.params.b_iso_min or
                 b_isos[i_seq] > self.params.b_iso_max): keep = False
              good_map = self._maps.score_atom(
                atom      = atom,
                min_cc    = mfp.poor_cc_threshold,
                min_value = mfp.poor_map_value_threshold)
              if(not good_map): keep = False
              #
              # XXX It does not work for 1f8t. WHY? XXX
              #
              #o = maptbx.sphericity_by_heuristics(
              #  map_data    = self._maps.map_2,
              #  unit_cell   = self.model.get_xray_structure().unit_cell(),
              #  center_cart = atom.xyz,
              #  radius      = 1.0)
              #if(o.ccs[0] < 0.1): keep = False
              #
              if(not keep):
                selection = selection.set_selected(i_seqs, False)
    # Apply selection
    self.model  = self.model.select(selection)
    n_sol_final = self.model.solvent_selection().count(True)
    if(n_sol_final != n_sol_start):
      self.fmodel.update_xray_structure(
        xray_structure = self.model.get_xray_structure(),
        update_f_calc  = True,
        update_f_mask  = True)

  def _refine(self):
    if(self.params.mode == "filter_only"): return
    if(self.model.size() == self.model_size_init or self.n_water==0): return
    for i in range(self.params.n_cycles):
      self.refine_adp()
      self.refine_occupancies()

  def _find_peaks(self):
    if(self.params.mode == "filter_only"): return
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    self.find_peaks_params.max_number_of_peaks=self.model.get_number_of_atoms()
    assert self.params.primary_map_type is not None
    self._peaks = find_peaks.manager(
      fmodel     = self.fmodel,
      map_type   = self.params.primary_map_type,
      map_cutoff = self.params.primary_map_cutoff,
      params     = self.find_peaks_params,
      log        = null_out()).peaks_mapped()

  def _add_new_solvent(self):
    if(self._peaks is None): return
    sites, heights = self._peaks.sites, self._peaks.heights
    if(sites.size()==0): return
    if(self.params.mode == "filter_only"): return
    if(self.params.b_iso is None):
      b_solv = min(max(
        self.params.b_iso_min, flex.mean(self.model.get_b_iso())),
          self.params.b_iso_max)
    else:
      b_solv = self.params.b_iso
    if(self.params.new_solvent == "isotropic"):
      new_scatterers = flex.xray_scatterer(
        sites.size(),
        xray.scatterer(occupancy       = self.params.occupancy,
                       b               = b_solv,
                       scattering_type = self.params.scattering_type))
    elif(self.params.new_solvent == "anisotropic"):
      u_star = adptbx.u_iso_as_u_star(
        self.model.crystal_symmetry().unit_cell(), adptbx.b_as_u(b_solv))
      new_scatterers = flex.xray_scatterer(
        sites.size(),
        xray.scatterer(
          occupancy       = self.params.occupancy,
          u               = u_star,
          scattering_type = self.params.scattering_type))
    else: raise RuntimeError
    new_scatterers.set_sites(sites)
    solvent_xray_structure = xray.structure(
      special_position_settings = self.model.get_xray_structure(),
      scatterers                = new_scatterers)
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True,
      update_f_mask  = True)

  def refine_adp(self):
    if(self.params.refine_adp #and
       #self.model.refinement_flags.individual_adp and
       ):
      hd_sel     = self.model.get_hd_selection()
      not_hd_sel = ~hd_sel
      sol_sel    = self.model.solvent_selection()
      not_sol_sel= ~sol_sel
      selection_aniso = self.model.get_xray_structure().use_u_aniso().deep_copy()
      if(self.params.new_solvent == "anisotropic"):
        selection_aniso.set_selected(sol_sel, True)
      selection_iso   = self.model.get_xray_structure().use_u_iso().deep_copy()
      selection_aniso.set_selected(not_sol_sel, False)
      selection_iso  .set_selected(not_sol_sel, False)
      if(not self.is_neutron_scat_table):
        selection_aniso.set_selected(hd_sel, False)
        selection_iso.set_selected(hd_sel, False)
      selection_aniso.set_selected(selection_iso, False)
      selection_iso.set_selected(selection_aniso, False)
      self.model.set_refine_individual_adp(
        selection_aniso = selection_aniso, selection_iso = selection_iso)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = 50)
      minimized = minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = self.fmodels,
        model                    = self.model,
        is_neutron_scat_table    = self.is_neutron_scat_table,
        refine_adp               = True,
        lbfgs_termination_params = lbfgs_termination_params)
      self.model.adopt_xray_structure(
        xray_structure = self.fmodel.xray_structure)
      print("  ADP (water only), final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def refine_occupancies(self):
    if(self.params.refine_occupancies #and
       #self.model.refinement_flags.occupancies and
       ):
      self.fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
        state = False)
      i_selection = self.model.solvent_selection().iselection()
      self.fmodels.fmodel_xray().xray_structure.scatterers(
        ).flags_set_grad_occupancy(iselection = i_selection)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = 50)
      minimized = mmtbx.refinement.minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = self.fmodels,
        model                    = self.model,
        is_neutron_scat_table    = self.is_neutron_scat_table,
        lbfgs_termination_params = lbfgs_termination_params)
      #
      self.model.adopt_xray_structure(
        xray_structure = self.fmodel.xray_structure)
      print("  occ (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def _correct_drifted_waters(self):
    if(self.params.mode != "filter_only"): return
    if(not self.params.correct_drifted_waters): return
    map_cutoff = self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold/2
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    peaks = find_peaks.manager(
      fmodel         = self.fmodel,
      map_type       = "2mFobs-DFmodel",
      map_cutoff     = map_cutoff,
      params         = find_peaks_params_drifted,
      log            = null_out()).peaks_mapped()
    if(peaks is not None and self.fmodel.r_work() > 0.01):
      sites_frac, heights = peaks.sites, peaks.heights
      model_sites_frac = self.model.get_xray_structure().sites_frac()
      solvent_selection = self.model.solvent_selection()
      mmtbx.utils.correct_drifted_waters(
        sites_frac_all   = model_sites_frac,
        sites_frac_peaks = sites_frac,
        water_selection  = solvent_selection,
        unit_cell        = self.model.get_xray_structure().unit_cell())
      self.model.get_xray_structure().set_sites_frac(sites_frac = model_sites_frac)
      self.fmodel.update_xray_structure(
        xray_structure = self.model.get_xray_structure(),
        update_f_calc  = True,
        update_f_mask  = True)
