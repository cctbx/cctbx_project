from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import xray
import math,sys
from libtbx import adopt_init_args
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.phil
from mmtbx import find_peaks
import mmtbx.utils
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from six.moves import range
from libtbx.utils import user_plus_sys_time
from libtbx.utils import null_out
from mmtbx.solvent import map_to_water
from libtbx import group_args
import string
import libtbx.log

def get_unique_altloc(exclude):
  for l in string.ascii_uppercase:
    if not l in exclude:
      return l

def get_unique_altloc2(available, exclude):
  l = None
  for l in available:
    if l in [' ', '']: continue
    if not l in exclude:
      return l
  return l

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
  dist_min = 1.8
    .type = float
    .short_caption = Min distance between water and any atom
    .expert_level = 1
  dist_max = 3.2
    .type = float
    .short_caption = Max distance between water and any atom
    .expert_level = 1
  dist_min_altloc = 0.5
    .type = float
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
  occupancy_min = 0.2
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
  mask_atoms_selection = protein and (name CA or name CB or name N or name C or name O)
    .type = str
    .help = Mask macromolecule atoms in peak picking map
  include_altlocs = False
    .type = bool
    .help = Search for water with altlocs
  n_cycles = 1
    .type = int
    .short_caption = Number of cycles
  %s
  primary_map_type = mFobs-DFmodel
    .type=str
    .help = Map used to identify candidate water sites - by default this is \
      the standard difference map.
  primary_map_cutoff = 3.
    .type=float
    .short_caption = Primary map cutoff (sigma)
  secondary_map_and_map_cc_filter
    .short_caption = Secondary map filter
    .style = auto_align box
  {
    cc_map_1_type = "Fmodel"
      .type = str
      .short_caption = Model map type for CC calculation
    cc_map_2_type = 2mFobs-DFmodel
      .type = str
      .short_caption = Experimental map type for CC calculation
    poor_cc_threshold = 0.70
      .type = float
      .short_caption = Minimum map correlation
    poor_map_value_threshold = 1.0
      .type = float
      .short_caption = Minimum map value (sigma)
  }
  %s
  refine_oat = False
    .type = bool
    .help = Q & B coarse grid search.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_occupancies = True
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

def master_params():
  return iotbx.phil.parse(master_params_str)

class maps(object):
  def __init__(self,
               fmodel,
               model,
               mask_atoms_selection,
               difference_map_type,
               model_map_type,
               data_map_type,
               grid_step=0.6,
               radius=2.0):
    self.radius = radius
    self.fmodel = fmodel
    self.model  = model
    self.e_map = fmodel.electron_density_map()
    self.crystal_symmetry = fmodel.xray_structure.crystal_symmetry()
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell        = self.crystal_symmetry.unit_cell(),
      space_group_info = self.crystal_symmetry.space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      step             = grid_step)
    self.difference_map = self._get_real_map(map_type = difference_map_type)
    self.model_map      = self._get_real_map(map_type = model_map_type)
    self.data_map       = self._get_real_map(map_type = data_map_type)
    #self._estimate_diff_map_cutoff()
    # Compute mask in P1 to mask out desired regions
    if mask_atoms_selection is not None:
      bb_sel = self.model.selection(string=mask_atoms_selection)
      if bb_sel.count(True)>0:
        xrs = fmodel.xray_structure.select(bb_sel)
        #
        # Both ways to compute mask should be the same, but they are slighly not,
        # expectedly.
        #
        mask_p1 = mmtbx.masks.mask_from_xray_structure(
          xray_structure           = xrs,
          p1                       = True,
          for_structure_factors    = True,
          solvent_radius           = 0,
          shrink_truncation_radius = 0,
          atom_radius              = 1.2,
          n_real                   = self.crystal_gridding.n_real(),
          in_asu                   = False).mask_data
        maptbx.unpad_in_place(map=mask_p1)
        sel0 = mask_p1 > 0.1
        mask_p1 = mask_p1.set_selected(sel0,  1)
        mask_p1 = mask_p1.set_selected(~sel0, 0)

        #xrs = xrs.expand_to_p1()
        #import boost_adaptbx.boost.python as bp
        #cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
        #radii = flex.double(xrs.scatterers().size(), 1.2)
        #mask_p1 = cctbx_maptbx_ext.mask(
        #  sites_frac                  = xrs.sites_frac(),
        #  unit_cell                   = xrs.unit_cell(),
        #  n_real                      = self.crystal_gridding.n_real(),
        #  mask_value_inside_molecule  = 0,
        #  mask_value_outside_molecule = 1,
        #  radii                       = radii,
        #  wrapping                    = True)

        self.difference_map = self.difference_map * mask_p1

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
      fft_n_real = self.model_map.focus(),
      fft_m_real = self.model_map.all(),
      sites_cart = flex.vec3_double([site_cart]),
      site_radii = flex.double([self.radius]))
    cc = flex.linear_correlation(
      x=self.model_map.select(sel),
      y=self.data_map.select(sel)).coefficient()
    value_2 = self.data_map.eight_point_interpolation(site_frac)
    diff_map_val = self.difference_map.eight_point_interpolation(site_frac)
    result = (cc > min_cc and value_2 > min_value*atom.occ) and \
             not (diff_map_val < -3)
    return group_args(result = result, cc=cc, value_2=value_2)

  def _estimate_diff_map_cutoff(self):
    scatterers = self.fmodel.xray_structure.scatterers()
    sel = flex.random_bool(scatterers.size(), 0.1)
    result = flex.double()
    cntr=0
    for s, sc in zip(sel, scatterers):
      if not s: continue
      print(self.fmodel.r_work())
      occ = sc.occupancy
      sc.occupancy=0
      self.fmodel.update_xray_structure(update_f_calc=True)

      coeffs = self.fmodel.electron_density_map().map_coefficients(
        map_type     = "mFobs-DFmodel",
        fill_missing = False,
        isotropize   = True)
      fft_map = coeffs.fft_map(crystal_gridding = self.crystal_gridding)
      fft_map.apply_sigma_scaling()
      map_data = fft_map.real_map_unpadded()
      mv = map_data.tricubic_interpolation(sc.site)
      result.append(mv)

      print(self.fmodel.r_work(), mv)
      sc.occupancy=occ
      self.fmodel.update_xray_structure(update_f_calc=True)
      print(self.fmodel.r_work())
      print()
      if cntr>100: break
      cntr+=1
    mean = flex.mean(result)
    if mean/2 > 3: cutoff = 3

def fix_altlocs_and_filter(model, dist_min=1.8, fix_only=False):
  present_altlocs = list(
    model.get_hierarchy().get_conformer_indices().index_altloc_mapping.keys())
  sps = model.crystal_symmetry().special_position_settings()
  get_class = iotbx.pdb.common_residue_names_get_class
  only_model = model.get_hierarchy().only_model()
  eps = 1.e-3
  dist_min = dist_min-eps
  sites_cart = model.get_sites_cart()
  atoms = only_model.atoms()
  remove_sel = flex.size_t()
  for agi in only_model.atom_groups(): # loop over water
    if(not get_class(name=agi.resname) == "common_water"): continue
    for ai in agi.atoms():
      if ai.element_is_hydrogen(): continue # skip water H
      #
      selection_around_ai = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_min, i_seq=ai.i_seq)
      if len(selection_around_ai) == 0: continue
      #
      altlocs_inside = []
      for j in selection_around_ai:
        altlocs_inside.append( atoms[j].parent().altloc )
      #
      skip = False
      for j in selection_around_ai:
        aj = atoms[j]
        if aj.element_is_hydrogen(): continue
        agj = aj.parent()
        if(agj.altloc in [' ', '']):
          if(get_class(name=agj.resname) != "common_water"):
            remove_sel.extend(agi.atoms().extract_i_seq())
            skip=True
          else:
            new_altloc = get_unique_altloc2(
              available = present_altlocs,
              exclude   = altlocs_inside+[agi.altloc])
            if new_altloc is not None: agj.altloc = new_altloc
            else:
              remove_sel.extend(agi.atoms().extract_i_seq())
              skip = True
      if skip: continue
      #
      altlocs_inside = []
      for j in selection_around_ai:
        altlocs_inside.append( atoms[j].parent().altloc )
      #
      if agi.altloc in altlocs_inside or agi.altloc in [' ', '']:
        #print(ai.i_seq, selection_around_ai, altlocs_inside, [agi.altloc])
        #new_altloc = get_unique_altloc(exclude=altlocs_inside+[agi.altloc])
        new_altloc = get_unique_altloc2(
              available = present_altlocs,
              exclude   = altlocs_inside+[agi.altloc])
        if new_altloc is not None:
          agi.altloc = new_altloc
        else:
          remove_sel.extend(agi.atoms().extract_i_seq())
  #
  if remove_sel.size() > 0 and not fix_only:
    remove_sel = ~flex.bool(model.size(), remove_sel)
    model = model.select(selection = remove_sel)
  return model

def get_sphere_selection(
      sites_cart_all, special_position_settings, radius, i_seq):
  sel = flex.bool(sites_cart_all.size(), False)
  sel[i_seq] = True
  selection_around_i_seq = special_position_settings.pair_generator(
    sites_cart      = sites_cart_all,
    distance_cutoff = radius
      ).neighbors_of(primary_selection = sel).iselection()
  selection_around_i_seq = list(selection_around_i_seq)
  if len(selection_around_i_seq) > 0:
    selection_around_i_seq.remove(i_seq)
  return selection_around_i_seq

def filter_by_distance(model, fix_altlocs_and_filter_was_run, dist_min=1.8,
                       dist_max=3.2):
  interaction_selection = model.selection(
    map_to_water.selection_string_interaction)
  sps = model.crystal_symmetry().special_position_settings()
  get_class = iotbx.pdb.common_residue_names_get_class
  only_model = model.get_hierarchy().only_model()
  sites_cart = model.get_sites_cart()
  atoms = only_model.atoms()
  remove_sel = flex.size_t()
  for agi in only_model.atom_groups(): # loop over water
    if(not get_class(name=agi.resname) == "common_water"): continue
    for ai in agi.atoms():
      if ai.element_is_hydrogen(): continue # skip water H
      # Get selections
      selection_around_ai_min = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_min, i_seq=ai.i_seq)
      selection_around_ai_max = get_sphere_selection(
        sites_cart_all=sites_cart, special_position_settings=sps,
        radius=dist_max, i_seq=ai.i_seq)
      selection_shell = list( set(selection_around_ai_max) -
                              set(selection_around_ai_min) )
      #
      altloc_i = ai.parent().altloc.strip()
      # Make sure anything inside smaller sphere are altlocs
      for j in selection_around_ai_min:
        aj = atoms[j]
        altloc_j = aj.parent().altloc
        if fix_altlocs_and_filter_was_run:
          if altloc_i =="" or altloc_i==" " or altloc_j =="" or altloc_j==" " or altloc_i==altloc_j:
            remove_sel.extend(agi.atoms().extract_i_seq())

      # Check water inside shell dist_min < dist < dist_max
      found = False
      for j in selection_shell:
        if not interaction_selection[j]: continue
        aj = atoms[j]
        if aj.element_is_hydrogen(): continue
        agj = aj.parent()
        altloc_j = agj.altloc.strip()
        if altloc_i=="" or altloc_j=="":        found = True
        if altloc_i==altloc_j and altloc_i!="": found = True
        if get_class(name=agj.resname) == "common_water":
          if altloc_i!="" and altloc_j!="": found = True
      if not found:
        remove_sel.extend(agi.atoms().extract_i_seq())
  #
  if remove_sel.size() > 0:
    remove_sel = ~flex.bool(model.size(), remove_sel)
    model = model.select(selection = remove_sel)
  return model

class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     is_neutron_scat_table,
                     params = master_params().extract(),
                     find_peaks_params = None,
                     log = sys.stdout):
    adopt_init_args(self, locals())

    # XXX Rationalize this:

    self.find_peaks_params.map_next_to_model.min_peak_peak_dist=self.params.dist_max
    if self.params.include_altlocs:
      self.find_peaks_params.peak_search.min_cross_distance=0.5
      self.find_peaks_params.map_next_to_model.min_model_peak_dist=0.5
      self.find_peaks_params.map_next_to_model.min_peak_peak_dist=0.5

    self.ma         = libtbx.log.manager(log = self.log)
    self.total_time = 0
    self._maps      = None
    self._peaks     = None
    self.n_water    = None
    self.model_size_init = self.model.size()
    self.new_solvent_selection = None
    #
    self._call(msg="Start",            func=None)
    self._call(msg="Filter (dist)",    func=self._filter_dist_fix_altlocs)
    self._call(msg="Filter (q & B)",   func=self._filter_q_b)
    self._call(msg="Compute maps",     func=self._get_maps)
    self._call(msg="Filter (map)",     func=self._filter_map)
    self._call(msg="Find peaks",       func=self._find_peaks)
    self._call(msg="Add new water",    func=self._add_new_solvent)
    self._call(msg="Refine new water", func=self._refine)
    self._call(msg="Filter (q & B)",   func=self._filter_q_b)
    self._call(msg="Filter (dist only)",   func=self._filter_dist)
    #self._call(msg="Correct drifted",  func=self._correct_drifted_waters)

  def _call(self, msg, func = None):
    timer = user_plus_sys_time()
    self.ma.add_and_show(msg)
    self._assert_same_model()
    if(func is not None): func()
    self._get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps()
    self._assert_same_model()
    t = timer.elapsed()
    self.total_time += t
    self._add_to_message(this_step_time=t)

  def _add_to_message(self, this_step_time):
    rs="r_work=%6.4f r_free=%6.4f"%(self.fmodel.r_work(), self.fmodel.r_free())
    nw="n_water=%3d"%(self.n_water)
    self.ma.add_and_show("  %s | %s | time (s): %s (total time: %s)"%(rs, nw,
      ("%8.3f"%this_step_time).strip(), ("%8.3f"%self.total_time).strip()))

  def _get_and_set_n_water_and_sync_fmodel_and_model_and_update_maps(self):
    n_water = self.model.solvent_selection().count(True)
    if n_water!=self.n_water:
      self.fmodel.update_xray_structure(
        xray_structure = self.model.get_xray_structure(),
        update_f_calc  = True,
        update_f_mask  = True)
      self._get_maps()
    self.n_water = n_water

  def _assert_same_model(self):
    mmtbx.utils.assert_xray_structures_equal( # XXX MAKE METHOD OF XRS
      x1 = self.model.get_xray_structure(),
      x2 = self.fmodel.xray_structure,
      eps=1.e-3)

  def _get_maps(self):
    p = self.params
    self._maps = maps(
      fmodel               = self.fmodel,
      model                = self.model,
      mask_atoms_selection = p.mask_atoms_selection,
      difference_map_type  = p.primary_map_type,
      model_map_type       = p.secondary_map_and_map_cc_filter.cc_map_1_type,
      data_map_type        = p.secondary_map_and_map_cc_filter.cc_map_2_type)

  def _filter_dist_fix_altlocs(self):
    if self.params.include_altlocs:
      self.model = fix_altlocs_and_filter(
        model    = self.model,
        dist_min = self.params.dist_min)
    self.model = filter_by_distance(
      model                          = self.model,
      fix_altlocs_and_filter_was_run = self.params.include_altlocs,
      dist_min                       = self.params.dist_min,
      dist_max                       = self.params.dist_max)

  def _filter_dist(self):
    self.model = filter_by_distance(
      model                          = self.model,
      fix_altlocs_and_filter_was_run = self.params.include_altlocs,
      dist_min                       = self.params.dist_min,
      dist_max                       = self.params.dist_max)

  def _filter_q_b(self):
    self._filter(filter_occ=True, filter_adp=True)

  def _filter_map(self):
    self._filter(filter_map=True)

  def _filter(self,
              filter_map=False,
              filter_occ=False,
              filter_adp=False):
    mfp = self.params.secondary_map_and_map_cc_filter
    selection = flex.bool(self.model.size(), True)
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
            has_oxygen = False # catch only H/D with no O water
            for atom in r.atoms():
              if(atom.element.strip().upper()=="O"): has_oxygen = True
              if(atom.element_is_hydrogen()): continue
              i_seq = atom.i_seq
              # Occupancy
              if filter_occ:
                if(atom.occ > self.params.occupancy_max or
                   atom.occ < self.params.occupancy_min): keep = False
                assert approx_equal(atom.occ, occ[i_seq], 1.e-3)
              # ADP
              if filter_adp:
                if(anisotropy[i_seq] < self.params.anisotropy_min): keep = False
                if(b_isos[i_seq] < self.params.b_iso_min or
                   b_isos[i_seq] > self.params.b_iso_max): keep = False
              #
              if filter_map:
                good_map = self._maps.score_atom(
                  atom      = atom,
                  min_cc    = mfp.poor_cc_threshold,
                  min_value = mfp.poor_map_value_threshold)
                if(not good_map.result):
                  keep = False
            if(not has_oxygen): keep=False
            if(not keep):
              selection = selection.set_selected(i_seqs, False)
    self.model = self.model.select(selection)

  def _refine(self):
    if(self.params.mode == "filter_only"): return
    if(self.model.size() == self.model_size_init or self.n_water==0): return
    for i in range(self.params.n_cycles):
      self.refine_oat()

  def _find_peaks(self):
    if(self.params.mode == "filter_only"): return
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    self.find_peaks_params.max_number_of_peaks=self.model.get_number_of_atoms()
    assert self.params.primary_map_type is not None
    self._peaks = find_peaks.manager(
      xray_structure = self.fmodel.xray_structure,
      map_data       = self._maps.difference_map, # diff-map
      map_cutoff     = self.params.primary_map_cutoff,
      params         = self.find_peaks_params,
      log            = null_out()).peaks_mapped()

  def _write_pdb_file(self, file_name="tmp.pdb", sites_frac=None):
    if sites_frac is not None:
      fmt = "ATOM  %5d  O   HOH S%4d    %8.3f%8.3f%8.3f  1.00  0.00           O"
      uc = self.fmodel.xray_structure.crystal_symmetry().unit_cell()
      with open(file_name, "w") as fo:
        for i, sf in enumerate(sites_frac):
          sc = uc.orthogonalize(sf)
          print(fmt%(i,i,sc[0],sc[1],sc[2]), file=fo)
    else:
      with open(file_name, "w") as fo:
        fo.write(self.model.model_as_pdb())

  def _add_new_solvent(self):
    if(self._peaks is None): return
    sites, heights = self._peaks.sites, self._peaks.heights
    if(sites.size()==0): return
    if(self.params.mode == "filter_only"): return

    self.new_solvent_selection = flex.bool(self.model.size(), False)
    self.new_solvent_selection.extend(flex.bool(sites.size(), True))

    if self.params.refine_oat:
      b_solv = 0
      occ = 0.
    else:
      b_solv = 20
      if self.params.refine_occupancies:
        occ = 0.004
      else:
        occ = 1.

    if(self.params.new_solvent == "isotropic"):
      new_scatterers = flex.xray_scatterer(
        sites.size(),
        xray.scatterer(occupancy       = occ,
                       b               = b_solv,
                       scattering_type = self.params.scattering_type))
    elif(self.params.new_solvent == "anisotropic"):
      u_star = adptbx.u_iso_as_u_star(
        self.model.crystal_symmetry().unit_cell(), adptbx.b_as_u(b_solv))
      new_scatterers = flex.xray_scatterer(
        sites.size(),
        xray.scatterer(
          occupancy       = occ,
          u               = u_star,
          scattering_type = self.params.scattering_type))
    else: raise RuntimeError
    new_scatterers.set_sites(sites)
    solvent_xray_structure = xray.structure(
      special_position_settings = self.model.get_xray_structure(),
      scatterers                = new_scatterers)
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      conformer_indices      = None, #self._peaks.conformer_indices,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)
    ####
    # This is an ugly work-around to set altlocs and conformer_indices
    ####
    if self.params.include_altlocs:
      self.model = fix_altlocs_and_filter(model=self.model, fix_only=True)
      ss = self.model.solvent_selection()
      ms = self.model.select(ss)
      self.model = self.model.select(~ss)
      self.model.add_solvent(
        solvent_xray_structure = ms.get_xray_structure(),
        conformer_indices      = ms.get_hierarchy().get_conformer_indices(),
        residue_name           = self.params.output_residue_name,
        atom_name              = self.params.output_atom_name,
        chain_id               = self.params.output_chain_id,
        refine_occupancies     = self.params.refine_occupancies,
        refine_adp             = self.params.new_solvent)
    ###

  def refine_oat(self):
    if self.new_solvent_selection is None: return
    if(self.params.refine_oat and self.new_solvent_selection.count(True)>0):
      from phenix.programs import oat
      from cctbx import adptbx
      atoms = self.model.get_hierarchy().atoms()
      scatterers = self.fmodel.xray_structure.scatterers()
      for i, sel in enumerate(self.new_solvent_selection):
        if not sel: continue
        r_start = self.fmodel.r_work()
        scatterers[i].occupancy=0
        scatterers[i].u_iso=0
        self.fmodel.update_xray_structure(update_f_calc=True)
        r_omit = self.fmodel.r_work()
        fmodel_dc = self.fmodel.deep_copy()
        oo = oat.loop(
          fmodel  = fmodel_dc,
          site    = atoms[i].xyz,
          label   = "O",
          qs = flex.double([0.3, 0.6, 0.9]),
          bs = flex.double([10, 30, 60])
          )
        scatterers[i].occupancy = oo.o_best
        scatterers[i].u_iso     = adptbx.b_as_u(oo.b_best)
        self.fmodel.update_xray_structure(update_f_calc=True)
        r_final = self.fmodel.r_work()
        #print(atoms[i].quote(), "%8.6f %8.6f %8.6f %8.6f"%(
        #  r_start, r_omit, r_final, oo.rw_best), "|", oo.o_best, oo.b_best)
      self.model.adopt_xray_structure(
        xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), OAT, final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)
      #
    if([self.params.refine_adp , self.params.refine_occupancies].count(True)>0 and
       self.new_solvent_selection.count(True)>0):
      from mmtbx.refinement import wrappers
      o = wrappers.unrestrained_qbr_fsr(
        fmodel     = self.fmodel,
        model      = self.model,
        refine_xyz = False,
        refine_q   = self.params.refine_occupancies,
        refine_b   = self.params.refine_adp,
        selection  = self.new_solvent_selection,
        q_min      = 0.004,
        b_max      = 60,
        b_min      = self.params.b_iso_min,
        log        = self.log)
      self.model.adopt_xray_structure(
          xray_structure = self.fmodel.xray_structure)
      print("  ADP+occupancy (water only), MIN, final r_work=%6.4f r_free=%6.4f"%(
          self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def _correct_drifted_waters(self):
    if(self.params.mode != "filter_only"): return
    if(not self.params.correct_drifted_waters): return
    sol_sel = self.model.solvent_selection()
    hd_sel  = self.model.get_hd_selection()
    hd_sol = sol_sel & hd_sel
    if(hd_sol.count(True)>0): return
    map_cutoff = self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold/2
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    peaks = find_peaks.manager(
      xray_structure = self.fmodel.xray_structure,
      map_data       = self._maps.data_map,
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
