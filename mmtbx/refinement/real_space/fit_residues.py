from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import sys, time
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

class run(object):
  def __init__(self,
               pdb_hierarchy,
               crystal_symmetry,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               map_data_scale = 2.5,
               diff_map_data_threshold = -2.5,
               rotatable_hd=None,
               bselection=None,
               map_data=None,
               vdw_radii=None,
               backbone_sample=False, # XXX
               diff_map_data=None,
               log = None):
    adopt_init_args(self, locals())
    self.processed = 0
    self.total_time_residue_loop = 0
    t0 = time.time()
    self.mes = []
    if(self.bselection is None):
      o = mmtbx.refinement.real_space.side_chain_fit_evaluator(
        pdb_hierarchy           = self.pdb_hierarchy,
        crystal_symmetry        = self.crystal_symmetry,
        rotamer_evaluator       = self.rotamer_manager.rotamer_evaluator,
        map_data                = self.map_data,
        diff_map_data           = self.diff_map_data,
        map_data_scale          = self.map_data_scale,
        diff_map_data_threshold = self.diff_map_data_threshold)
      self.mes.extend(o.mes)
      self.bselection = o.sel_all() # or all_possible() ?
    if(self.log is None): self.log = sys.stdout
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.atom_names = flex.std_string(
      [i.strip() for i in self.pdb_hierarchy.atoms().extract_name()])
    self.special_position_settings = crystal.special_position_settings(
      crystal_symmetry = self.crystal_symmetry)
    # Even better would be to pass it here. Ideally just use model
    asc = self.pdb_hierarchy.atom_selection_cache()
    self.selection_water_as_set = \
      set(asc.selection(string = "water").iselection())
    self.target_map = map_data
    self.mes.append("outliers start: %d"%self.count_outliers())
    #
    self.loop(function = self.one_residue_iteration)
    assert approx_equal(self.sites_cart,
      self.pdb_hierarchy.atoms().extract_xyz())
    self.mes.append("outliers after fit: %d"%self.count_outliers())
    #
    self.mes.append("outliers final: %d"%self.count_outliers())
    assert approx_equal(self.sites_cart,
      self.pdb_hierarchy.atoms().extract_xyz())
    #
    self.mes.append("residues processed: %d"%self.processed)
    if(self.processed > 0):
      self.mes.append("average time/residue: %6.4f"%(
        self.total_time_residue_loop/self.processed))
    self.mes.append("time to fit residues: %6.4f"%(time.time()-t0))

  def show(self, prefix=""):
    for m in self.mes:
      print("%s%s"%(prefix,m), file=self.log)
    self.log.flush()

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

  def one_residue_iteration(self, residue):
    t0 = time.time()
    negate_rad = negate_map_table[residue.resname.strip().lower()]
    if(not negate_rad): return
    xyzrad_bumpers = self.get_nonbonded_bumpers(
      residue=residue, radius=negate_rad)
    self.processed +=1
    mmtbx.refinement.real_space.fit_residue.run(
      residue           = residue,
      vdw_radii         = self.vdw_radii,
      xyzrad_bumpers    = xyzrad_bumpers,
      backbone_sample   = self.backbone_sample,
      unit_cell         = self.crystal_symmetry.unit_cell(),
      target_map        = self.target_map,
      target_map_for_cb = self.target_map,
      mon_lib_srv       = self.mon_lib_srv,
      rotamer_manager   = self.rotamer_manager,
      rotatable_hd      = self.rotatable_hd,
      sin_cos_table     = self.sin_cos_table,
      log               = self.log)
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.total_time_residue_loop += (time.time()-t0)

  def loop(self, function):
    for model in self.pdb_hierarchy.models():
      for chain in model.chains():
        for residue_group in chain.residue_groups():
          conformers = residue_group.conformers()
          if(len(conformers)>1): continue
          for conformer in residue_group.conformers():
            residue = conformer.only_residue()
            if(not self.bselection[residue.atoms()[0].i_seq]): continue
            function(residue = residue)

  def count_outliers(self):
    o = mmtbx.refinement.real_space.side_chain_fit_evaluator(
      pdb_hierarchy     = self.pdb_hierarchy,
      crystal_symmetry  = self.crystal_symmetry,
      rotamer_evaluator = self.rotamer_manager.rotamer_evaluator)
    return o.cntr_outliers

# XXX Where is this used? Looks like severe duplication!
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
