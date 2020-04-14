from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import sys
from cctbx import maptbx
from cctbx import crystal
import boost.python
cctbx_maptbx_ext = boost.python.import_ext("cctbx_maptbx_ext")
from libtbx.test_utils import approx_equal
fit_ext = boost.python.import_ext("mmtbx_rotamer_fit_ext")

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

def get_mean_side_chain_density_value(hierarchy, map_data, unit_cell):
  results = flex.double()
  get_class = iotbx.pdb.common_residue_names_get_class
  main_chain = ["N","CA","CB","C","O"]
  for model in hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        if(len(conformers)>1): continue
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          if(get_class(residue.resname) != "common_amino_acid"): continue
          sites_cart = flex.vec3_double()
          atoms = residue.atoms()
          for a in atoms:
            if(a.name.strip() in main_chain): continue
            sites_cart.append(a.xyz)
          if(sites_cart.size()==0): continue
          mv = maptbx.real_space_target_simple(
            unit_cell   = unit_cell,
            density_map = map_data,
            sites_cart  = sites_cart)
          results.append(mv)
  return flex.mean_default(results,0)

def get_mean_side_chain_density_value_residue(residue, map_data, unit_cell):
  results = flex.double()
  get_class = iotbx.pdb.common_residue_names_get_class
  main_chain = ["N","CA","CB","C","O"]
  if(get_class(residue.resname) != "common_amino_acid"): return None
  sites_cart = flex.vec3_double()
  atoms = residue.atoms()
  for a in atoms:
    if(a.name.strip() in main_chain): continue
    sites_cart.append(a.xyz)
  if(sites_cart.size()==0): return None
  mv = maptbx.real_space_target_simple(
    unit_cell   = unit_cell,
    density_map = map_data,
    sites_cart  = sites_cart)
  results.append(mv)
  return flex.mean_default(results,0)

class run(object):
  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               rotatable_hd=None,
               outliers_only=False,
               bselection=None,
               map_data=None,
               vdw_radii=None,
               do_all=False,
               backbone_sample=False, # XXX
               diff_map_data=None,
               tune_up_only=False,
               log = None):
    adopt_init_args(self, locals())

#    ### DEBUG start
#    from libtbx import easy_pickle as ep
#    tmp = [
#               pdb_hierarchy,
#               crystal_symmetry,
#               rotamer_manager,
#               sin_cos_table,
#               mon_lib_srv,
#               outliers_only,
#               bselection,
#               map_data,
#               vdw_radii,
#               do_all,
#               backbone_sample,
#               diff_map_data,
#               tune_up_only,
#               log
#           ]
#    ep.dump(file_name="all.pkl", obj=tmp)
#    STOP()
#    ### DEBUG end
    self.did_it_for = 0
    if(self.do_all): assert not outliers_only
    if(self.outliers_only): assert not do_all
    self.number_of_outliers = None
    self.mean_side_chain_density = None
    if(self.map_data is not None):
      self.mean_side_chain_density = get_mean_side_chain_density_value(
        hierarchy = self.pdb_hierarchy,
        map_data  = self.map_data,
        unit_cell = self.crystal_symmetry.unit_cell())
    if(self.log is None): self.log = sys.stdout
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.atom_names = flex.std_string(
      [i.strip() for i in self.pdb_hierarchy.atoms().extract_name()])
    self.special_position_settings = None
    self.special_position_indices = None
    if(self.crystal_symmetry is not None):
      self.special_position_settings = crystal.special_position_settings(
        crystal_symmetry = self.crystal_symmetry)
      self.special_position_indices = self.get_special_position_indices()
    # Even better would be to pass it here. Ideally just use model
    self.atom_selection_cache = self.pdb_hierarchy.atom_selection_cache()
    self.selection_water_as_set = set(self.atom_selection_cache.\
        selection(string = "water").iselection())
    self.target_map = map_data
    print("outliers start: %d"%self.count_outliers(), file=self.log)
    #
    if(not self.tune_up_only):
      self.loop(function = self.one_residue_iteration)
      assert approx_equal(self.sites_cart,
        self.pdb_hierarchy.atoms().extract_xyz())
      print("outliers after map fit: %d"%self.count_outliers(), file=self.log)
    #
    if(not outliers_only):
      print("tune up", file=self.log)
      assert approx_equal(self.sites_cart,
        self.pdb_hierarchy.atoms().extract_xyz())
      self.loop(function = self.one_residue_tune_up)
    #
    print("outliers final: %d"%self.count_outliers(), file=self.log)
    assert approx_equal(self.sites_cart,
      self.pdb_hierarchy.atoms().extract_xyz())
    #
    print("Did it for:", self.did_it_for, file=self.log)
    #
#    ### DEBUG start
#    pdb_hierarchy.write_pdb_file(file_name="end1.pdb")
#    STOP()
#    ### DEBUG end

  def get_special_position_indices(self):
    site_symmetry_table = \
      self.special_position_settings.site_symmetry_table(
        sites_cart = self.sites_cart,
        unconditional_general_position_flags=(
          self.pdb_hierarchy.atoms().extract_occ() != 1))
    return site_symmetry_table.special_position_indices()

  def get_nonbonded_bumpers(self, residue, radius):
    if(self.special_position_settings is None): return None
    if(self.vdw_radii is None): return None
    residue_i_selection = flex.size_t()
    for a in residue.atoms():
      if(not a.name.strip() in ["N", "C", "O"]):
        residue_i_selection.append(a.i_seq)
    residue_b_selection = flex.bool(self.sites_cart.size(), residue_i_selection)
    selection_around_residue = self.special_position_settings.pair_generator(
      sites_cart      = self.sites_cart,
      distance_cutoff = radius
        ).neighbors_of(primary_selection = residue_b_selection).iselection()
    selection_around_residue_minus_residue = flex.size_t(
      list(set(selection_around_residue).difference(
        set(residue_i_selection)).difference(self.selection_water_as_set)))
    # exclude rotatable H
    selection_around_residue_minus_residue_minus_rotatableH = flex.size_t()
    for s in selection_around_residue_minus_residue:
      if(not self.rotatable_hd[s]):
        selection_around_residue_minus_residue_minus_rotatableH.append(s)
    sites_cart = self.sites_cart.select(
      selection_around_residue_minus_residue_minus_rotatableH)
    atom_names = self.atom_names.select(
      selection_around_residue_minus_residue_minus_rotatableH)
    #
    radii = flex.double()
    for an in atom_names:
      try: radii.append(self.vdw_radii[an]-0.25)
      except KeyError: radii.append(1.5) # XXX U, Uranium, OXT are problems!
    #
    return fit_ext.fixed(sites_cart = sites_cart, radii = radii)

  def on_special_position(self, residue):
    if(self.special_position_indices is None): return False
    if(self.special_position_indices.size()==0): return False
    for i_seq in residue.atoms().extract_i_seq():
      if(i_seq in self.special_position_indices): return True
    return False

  def one_residue_iteration(self, residue):
    if(residue.resname.strip().upper() in ["ALA","GLY","PRO"]): return
    if(self.on_special_position(residue=residue)): return
    if(not self.do_all and self.map_data is not None):
      need_fix = mmtbx.refinement.real_space.need_sidechain_fit(
        residue           = residue,
        rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
        mon_lib_srv       = self.mon_lib_srv,
        unit_cell         = self.crystal_symmetry.unit_cell(),
        outliers_only     = self.outliers_only,
        f_map             = self.map_data,
        fdiff_map         = self.diff_map_data)
      # XXX Totally ad hoc, to rationalize later!
      if(not self.outliers_only):
        mv_r = get_mean_side_chain_density_value_residue(
          residue   = residue,
          map_data  = self.map_data,
          unit_cell = self.crystal_symmetry.unit_cell())
        if(mv_r is not None and mv_r<self.mean_side_chain_density/2):
          need_fix = True
      #
      if(residue.resname in ["MSE", "MET"]): need_fix=True
      #
      if(not need_fix): return
    negate_rad = negate_map_table[residue.resname.strip().lower()]
    if(not negate_rad): return
    xyzrad_bumpers = self.get_nonbonded_bumpers(
      residue=residue, radius=negate_rad)
    unit_cell = None
    if(self.crystal_symmetry is not None):
      unit_cell = self.crystal_symmetry.unit_cell()
    self.did_it_for +=1
    mmtbx.refinement.real_space.fit_residue.run(
      residue           = residue,
      vdw_radii         = self.vdw_radii,
      xyzrad_bumpers    = xyzrad_bumpers,
      backbone_sample   = self.backbone_sample,
      unit_cell         = unit_cell,
      target_map        = self.target_map,
      target_map_for_cb = self.target_map,
      mon_lib_srv       = self.mon_lib_srv,
      rotamer_manager   = self.rotamer_manager,
      rotatable_hd      = self.rotatable_hd,
      sin_cos_table     = self.sin_cos_table,
      log               = self.log)
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()

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
            if self.bselection is not None:
              if not self.bselection[residue.atoms()[0].i_seq]: continue
            cntr += 1
            # Handy for debugging
            #print(chain.id, residue.resname, residue.resseq,"-"*10)
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
    if(total==0): return 0
    return result

  def one_residue_tune_up(self, residue):
    re = self.rotamer_manager.rotamer_evaluator
    if(re.evaluate_residue(residue)=="OUTLIER"):
      sites_cart_start = residue.atoms().extract_xyz()
      mv1 =  maptbx.real_space_target_simple(
        unit_cell   = self.crystal_symmetry.unit_cell(),
        density_map = self.target_map,
        sites_cart  = sites_cart_start)
      mmtbx.refinement.real_space.fit_residue.tune_up(
        residue         = residue,
        unit_cell       = self.crystal_symmetry.unit_cell(),
        target_map      = self.target_map,
        mon_lib_srv     = self.mon_lib_srv,
        rotamer_manager = self.rotamer_manager.rotamer_evaluator)
      mv2 =  maptbx.real_space_target_simple(
        unit_cell   = self.crystal_symmetry.unit_cell(),
        density_map = self.target_map,
        sites_cart  = residue.atoms().extract_xyz())
      den = ((abs(mv1)+abs(mv2))/2)*100.
      if(abs(den)>1.e-6):
        pshift = abs(abs(mv1)-abs(mv2)) / den
      else:
        pshift = None
      if(pshift is not None and pshift<3):
        self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      else:
        residue.atoms().set_xyz(sites_cart_start)

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
                    accept_only_if_max_shift_is_smaller_than = ac,
                    log = self.log)
