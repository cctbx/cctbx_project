from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import iotbx.ccp4_map
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import time, sys
from libtbx.test_utils import approx_equal
from cctbx import maptbx

negate_map_table = {
  "ala": False,
  "asn": 5,
  "asp": 5,
  "cys": False,
  "gln": 6,
  "glu": 6,
  "gly": False,
  "his": 6,
  "ile": 5,
  "leu": 5,
  "met": 6,
  "mse": 6,
  "phe": 6.5,
  "pro": False,
  "ser": False,
  "thr": False,
  "trp": 7.4,
  "tyr": 7.7,
  "val": False,
  "arg": 10,
  "lys": 8
}

class manager(object):
  def __init__(self,
               structure_monitor,
               rotamer_manager,
               sin_cos_table,
               mon_lib_srv,
               backbone_sample=True,
               log = None):
    adopt_init_args(self, locals())
    self.unit_cell = self.structure_monitor.xray_structure.unit_cell()
    if(self.log is None): self.log = sys.stdout
    assert approx_equal(
      self.structure_monitor.xray_structure.sites_cart(),
      self.structure_monitor.pdb_hierarchy.atoms().extract_xyz())
    self.special_position_indices = \
      self.structure_monitor.xray_structure.special_position_indices()
    self.atom_radius_to_negate_map_within = None
    if(self.structure_monitor.target_map_object.miller_array.d_min()>3.5):
      self.atom_radius_to_negate_map_within = 2.0
    else:
      self.atom_radius_to_negate_map_within = 1.5
    #
    self.selection_good = self.structure_monitor.map_cc_per_atom > 0.8
    self.selection_water = self.structure_monitor.pdb_hierarchy.\
      atom_selection_cache().selection(string = "water")
    self.iselection_backbone=flex.size_t()
    for r in self.structure_monitor.residue_monitors:
      if(r.selection_backbone is not None):
        self.iselection_backbone.extend(r.selection_backbone)
    self.loop_over_residues()

  def on_special_position(self, sm_residue):
    if(self.special_position_indices.size()==0): return False
    for i_seq in sm_residue.selection_all:
      if(i_seq in self.special_position_indices): return True
    return False

  def prepare_target_map(self): # XXX This may need to go external
    sm = self.structure_monitor
    target_map = sm.target_map_object.data.deep_copy()
    # truncate map
    selection = sm.pdb_hierarchy.atom_selection_cache().selection(
        string = "element C or element O or element N")
    mean_atom = flex.double()
    for i_sc, sc in enumerate(sm.xray_structure.scatterers()):
      if(selection[i_sc]):
        mean_atom.append(target_map.eight_point_interpolation(sc.site))
    mean_atom = flex.mean(mean_atom)
    target_map = target_map.set_selected(target_map>mean_atom, mean_atom)
    # Blend maps if applicable
    if(sm.target_map_object.f_map_diff is not None):
      diff_map = sm.target_map_object.f_map_diff.deep_copy()
      sel = diff_map < 2.
      diff_map = diff_map.set_selected(sel, 0)
      sel = diff_map > 3.
      diff_map = diff_map.set_selected(sel, 3)
      diff_map = diff_map/3.
      maptbx.combine_1(map_data=target_map, diff_map=diff_map)
    return target_map

  def loop_over_residues(self):
    sm = self.structure_monitor
    xrs = sm.xray_structure.deep_copy_scatterers()
    sites_cart = sm.xray_structure.sites_cart()
    get_class = iotbx.pdb.common_residue_names_get_class
    target_map = self.prepare_target_map()
    for i_res, r in enumerate(sm.residue_monitors):
      if(get_class(r.residue.resname) != "common_amino_acid"): continue
      go = mmtbx.refinement.real_space.need_sidechain_fit(
        residue           = r.residue,
        rotamer_evaluator = self.rotamer_manager.rotamer_evaluator,
        mon_lib_srv       = self.mon_lib_srv,
        unit_cell         = self.unit_cell,
        f_map             = sm.target_map_object.data,
        fdiff_map         = sm.target_map_object.f_map_diff)
      if(go):
        if(self.on_special_position(sm_residue=r)): continue
        #t0=time.time()
        iselection_n_external=None
        iselection_c_external=None
        if(i_res!=0 and i_res!=len(sm.residue_monitors)-1):
          iselection_c_external = sm.residue_monitors[i_res-1].selection_c
          iselection_n_external = sm.residue_monitors[i_res+1].selection_n
        elif(i_res==0 and len(sm.residue_monitors)>1):
          iselection_n_external = sm.residue_monitors[i_res+1].selection_n
        elif(i_res==len(sm.residue_monitors)-1 and len(sm.residue_monitors)>1):
          iselection_c_external = sm.residue_monitors[i_res-1].selection_c
        #print r.residue.resname, r.residue.resseq
        target_map_ = target_map
        negate_rad = negate_map_table[r.residue.resname.strip().lower()]
        if(negate_rad):
          negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
            xray_structure          = xrs,
            selection_within_radius = 5,#XXX
            iselection              = r.residue.atoms().extract_i_seq(),
            selection_good          = self.selection_good,
            iselection_backbone     = self.iselection_backbone,
            iselection_n_external   = iselection_n_external,
            iselection_c_external   = iselection_c_external)
          for i_w, w in enumerate(self.selection_water):
            if(w): negate_selection[i_w]=False
          #print "sel: %6.4f"%(time.time()-t0)
          target_map_ = mmtbx.refinement.real_space.\
            negate_map_around_selected_atoms_except_selected_atoms(
              xray_structure          = xrs,
              map_data                = target_map,
              negate_selection        = negate_selection,
              atom_radius             = self.atom_radius_to_negate_map_within)
          #print "neg: %6.4f"%(time.time()-t0)
          if 0:
            iotbx.ccp4_map.write_ccp4_map(
              file_name      = "AMap.map",
              unit_cell      = sm.xray_structure.unit_cell(),
              space_group    = sm.xray_structure.space_group(),
              gridding_first = (0,0,0),
              gridding_last  = tuple(sm.target_map_object.crystal_gridding.n_real()),
              map_data       = target_map_,
              labels         = flex.std_string(["DEBUG"]))
        mmtbx.refinement.real_space.fit_residue.run(
          residue         = r.residue,
          backbone_sample = self.backbone_sample,
          unit_cell       = xrs.unit_cell(),
          target_map      = target_map_,
          mon_lib_srv     = self.mon_lib_srv,
          rotamer_manager = self.rotamer_manager,
          sin_cos_table   = self.sin_cos_table)
        sites_cart = sites_cart.set_selected(r.residue.atoms().extract_i_seq(),
          r.residue.atoms().extract_xyz())
        xrs.set_sites_cart(sites_cart)
        #print "ref: %6.4f"%(time.time()-t0), r.residue.resname
    sm.update(xray_structure = xrs, accept_as_is=True)
    assert approx_equal(
      sm.xray_structure.sites_cart(),
      sm.pdb_hierarchy.atoms().extract_xyz())
    #
    result = tune_up(
      pdb_hierarchy   = sm.pdb_hierarchy,
      rotamer_manager = self.rotamer_manager.rotamer_evaluator,
      target_map      = target_map,
      unit_cell       = self.unit_cell,
      mon_lib_srv     = self.mon_lib_srv,
      log             = self.log)
    xrs.set_sites_cart(result.pdb_hierarchy.atoms().extract_xyz())
    sm.update(xray_structure = xrs, accept_as_is=True)

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
                    target_map      = None,
                    mon_lib_srv     = self.mon_lib_srv,
                    rotamer_manager = self.rotamer_manager,
                    sin_cos_table   = self.sin_cos_table,
                    accept_only_if_max_shift_is_smaller_than = ac)

class tune_up(object):
  def __init__(self,
               pdb_hierarchy,
               rotamer_manager,
               target_map,
               unit_cell,
               mon_lib_srv,
               log = None):
    adopt_init_args(self, locals())
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
              re = self.rotamer_manager
              if(re.evaluate_residue(residue)=="OUTLIER"):
                mmtbx.refinement.real_space.fit_residue.tune_up(
                  residue         = residue,
                  unit_cell       = self.unit_cell,
                  target_map      = target_map,
                  mon_lib_srv     = self.mon_lib_srv,
                  rotamer_manager = self.rotamer_manager)
