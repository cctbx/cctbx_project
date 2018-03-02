from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import iotbx.ccp4_map
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import sys
from cctbx import maptbx
from cctbx import crystal
import boost.python
cctbx_maptbx_ext = boost.python.import_ext("cctbx_maptbx_ext")

negate_map_table = {
  #"ala": False,
  "asn": 5,
  "asp": 5,
  "cys": 5,
  "gln": 6,
  "glu": 6,
  #"gly": False,
  "his": 6,
  "ile": 5,
  "leu": 5,
  "met": 6,
  "mse": 6,
  "phe": 6.5,
  #"pro": 5,#False,
  "ser": 5,
  "thr": 5,
  "trp": 7.4,
  "tyr": 7.7,
  "val": 5,
  "arg": 10,
  "lys": 8
}

class run(object):
  def __init__(self,
               pdb_hierarchy,
               map_data,
               crystal_symmetry,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               do_all=False,
               backbone_sample=True,
               diff_map_data=None,
               massage_map=True,
               tune_up_only=False,
               log = None):
    adopt_init_args(self, locals())
    self.number_of_outliers = None
    if(self.log is None): self.log = sys.stdout
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.special_position_settings = crystal.special_position_settings(
      crystal_symmetry = self.crystal_symmetry)
    self.special_position_indices = self.get_special_position_indices()
    # Even better would be to pass it here. Ideally just use model
    self.atom_selection_cache = self.pdb_hierarchy.atom_selection_cache()
    self.selection_water_as_set = set(self.atom_selection_cache.\
        selection(string = "water").iselection())
    if(self.massage_map):
      self.target_map = self.prepare_target_map()
    else:
      self.target_map = map_data
    print >> self.log, \
      "outliers start: %d (percent: %6.2f)"%self.count_outliers()
    #
    if(not self.tune_up_only):
      self.loop(function = self.one_residue_iteration)
      self.pdb_hierarchy.atoms().set_xyz(self.sites_cart)
      print >> self.log, \
        "outliers after map fit: %d (percent: %6.2f)"%self.count_outliers()
    print >> self.log, "tune up"
    self.loop(function = self.one_residue_tune_up)
    print >> self.log, \
      "outliers final: %d (percent: %6.2f)"%self.count_outliers()

  def get_special_position_indices(self):
    site_symmetry_table = \
      self.special_position_settings.site_symmetry_table(
        sites_cart = self.sites_cart,
        unconditional_general_position_flags=(
          self.pdb_hierarchy.atoms().extract_occ() != 1))
    return site_symmetry_table.special_position_indices()

  def compute_negate_mask(self, residue, radius):
    residue_i_selection = residue.atoms().extract_i_seq()
    residue_b_selection = flex.bool(self.sites_cart.size(), residue_i_selection)
    selection_around_residue = self.special_position_settings.pair_generator(
      sites_cart      = self.sites_cart,
      distance_cutoff = radius
        ).neighbors_of(primary_selection = residue_b_selection).iselection()
    selection_around_residue_minus_residue = flex.size_t(
      list(set(selection_around_residue).difference(
        set(residue_i_selection)).difference(self.selection_water_as_set)))
    sites_cart = self.sites_cart.select(selection_around_residue_minus_residue)
    sites_frac_p1 = self.crystal_symmetry.unit_cell().fractionalize(
      self.crystal_symmetry.expand_to_p1(sites_cart =
        self.sites_cart.select(selection_around_residue_minus_residue)))
    return cctbx_maptbx_ext.mask(
      sites_frac = sites_frac_p1,
      unit_cell  = self.crystal_symmetry.cell_equivalent_p1().unit_cell(),
      n_real     = self.target_map.all(),
      mask_value_inside_molecule = -10,
      mask_value_outside_molecule = 0,
      radii = flex.double(sites_frac_p1.size(), 2.))

  def on_special_position(self, residue):
    if(self.special_position_indices.size()==0): return False
    for i_seq in residue.atoms().extract_i_seq():
      if(i_seq in self.special_position_indices): return True
    return False

  def prepare_target_map(self): # XXX This may need to go external
    if self.map_data is None:
      # This just makes dummy map to allow functionality working without it.
      # Would prefer to just create all-zero map quickly, but couldn't find
      # how to do it.
      from cctbx import miller
      xrs = self.pdb_hierarchy.extract_xray_structure(crystal_symmetry=self.crystal_symmetry)
      xrs = xrs.select(self.atom_selection_cache.\
          selection("name C or name CA or name N or name O"))
      crystal_gridding = maptbx.crystal_gridding(
          unit_cell        = xrs.unit_cell(),
          space_group_info = xrs.space_group_info(),
          symmetry_flags   = maptbx.use_space_group_symmetry,
          d_min             = 10)
      fc = xrs.structure_factors(d_min = 10, algorithm = "fft").f_calc()
      fft_map = miller.fft_map(
          crystal_gridding=crystal_gridding,
          fourier_coefficients=fc)
      fft_map.apply_sigma_scaling()
      self.map_data = fft_map.real_map_unpadded(in_place=False)

    map_data = self.map_data.deep_copy()
    # truncate map
    selection = self.atom_selection_cache.selection(
      string = "element C or element O or element N")
    mean_atom = flex.double()
    for i_a, a in enumerate(list(self.pdb_hierarchy.atoms())):
      if(selection[i_a]):
        site_frac = self.crystal_symmetry.unit_cell().fractionalize(a.xyz)
        mean_atom.append(map_data.eight_point_interpolation(site_frac))
    mean_atom = flex.mean(mean_atom)
    map_data = map_data.set_selected(map_data>mean_atom, mean_atom)
    # Blend maps if applicable
    if(self.diff_map_data is not None):
      diff_map_data = self.diff_map_data.deep_copy()
      sel = diff_map_data < 2.
      diff_map_data = diff_map_data.set_selected(sel, 0)
      sel = diff_map_data > 3.
      diff_map_data = diff_map_data.set_selected(sel, 3)
      diff_map_data = diff_map_data/3.
      maptbx.combine_1(map_data=map_data, diff_map=diff_map_data)
    return map_data

  def one_residue_iteration(self, residue):
    if(residue.resname.strip().upper() in ["ALA","GLY","PRO"]): return
    if(self.on_special_position(residue=residue)): return
    if(not self.do_all):
      need_fix = mmtbx.refinement.real_space.need_sidechain_fit(
        residue           = residue,
        rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
        mon_lib_srv       = self.mon_lib_srv,
        unit_cell         = self.crystal_symmetry.unit_cell(),
        f_map             = self.map_data,
        fdiff_map         = self.diff_map_data)
      if(not need_fix): return
    negate_rad = negate_map_table[residue.resname.strip().lower()]
    if(not negate_rad): return
    mask = self.compute_negate_mask(residue=residue, radius=negate_rad)
    target_map_ = self.target_map + mask
    if 0:
      print r.residue.resseq, i_res
      iotbx.ccp4_map.write_ccp4_map(
        file_name      = "AMap.map",
        unit_cell      = self.crystal_symmetry.unit_cell(),
        space_group    = self.crystal_symmetry.space_group(),
        gridding_first = (0,0,0),
        gridding_last  = tuple(self.target_map.n_real()),
        map_data       = target_map_,
        labels         = flex.std_string(["DEBUG"]))
    mmtbx.refinement.real_space.fit_residue.run(
      residue           = residue,
      backbone_sample   = self.backbone_sample,
      unit_cell         = self.crystal_symmetry.unit_cell(),
      target_map        = target_map_,
      target_map_for_cb = self.target_map,
      mon_lib_srv       = self.mon_lib_srv,
      rotamer_manager   = self.rotamer_manager,
      sin_cos_table     = self.sin_cos_table)
    self.sites_cart = self.sites_cart.set_selected(
      residue.atoms().extract_i_seq(),
      residue.atoms().extract_xyz())

  def loop(self, function):
    get_class = iotbx.pdb.common_residue_names_get_class
    cntr = 0
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          if(len(conformers)>1): continue
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            if(get_class(residue.resname) != "common_amino_acid"): continue
            cntr += 1
            function(residue = residue)
    return cntr

  def count_outliers(self):
    self.number_of_outliers = 0
    def is_outlier(residue):
      re = self.rotamer_manager.rotamer_evaluator
      if(re.evaluate_residue(residue)=="OUTLIER"):
        self.number_of_outliers += 1
    total = self.loop(function = is_outlier)
    result = self.number_of_outliers
    self.number_of_outliers = None
    if(total==0): return 0, 0
    return result, result*100./total

  def one_residue_tune_up(self, residue):
    re = self.rotamer_manager.rotamer_evaluator
    if(re.evaluate_residue(residue)=="OUTLIER"):
      mmtbx.refinement.real_space.fit_residue.tune_up(
        residue         = residue,
        unit_cell       = self.crystal_symmetry.unit_cell(),
        target_map      = self.target_map,
        mon_lib_srv     = self.mon_lib_srv,
        rotamer_manager = self.rotamer_manager.rotamer_evaluator)

class fix_outliers(object):
  def __init__(self,
               pdb_hierarchy,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               f_map=None,
               fdiff_map=None,
               unit_cell=None,
               accept_only_if_max_shift_is_smaller_than=None,
               log = None):
    adopt_init_args(self, locals())
    assert [f_map, fdiff_map, unit_cell].count(None) in [0,3]
    ac = accept_only_if_max_shift_is_smaller_than
    if(self.log is None): self.log = sys.stdout
    get_class = iotbx.pdb.common_residue_names_get_class
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          if(len(conformers)>1): continue
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            id_str="%s%s%s"%(chain.id,residue.resname,residue.resseq.strip())
            if(get_class(residue.resname) == "common_amino_acid"):
              # Idealize rotamer: move to nearest rotameric state
              re = self.rotamer_manager.rotamer_evaluator
              if(re.evaluate_residue(residue)=="OUTLIER"):
                go=True
                if([f_map, fdiff_map, unit_cell].count(None)==0):
                  go = mmtbx.refinement.real_space.need_sidechain_fit(
                    residue           = residue,
                    rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
                    mon_lib_srv       = self.mon_lib_srv,
                    unit_cell         = self.unit_cell,
                    f_map             = f_map,
                    fdiff_map         = fdiff_map)
                if(go):
                  mmtbx.refinement.real_space.fit_residue.run(
                    residue         = residue,
                    backbone_sample = False,
                    mon_lib_srv     = self.mon_lib_srv,
                    rotamer_manager = self.rotamer_manager,
                    sin_cos_table   = self.sin_cos_table,
                    accept_only_if_max_shift_is_smaller_than = ac)
