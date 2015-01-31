from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import iotbx.ccp4_map
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
from mmtbx.rotamer.rotamer_eval import RotamerEval
import time, sys

class manager(object):
  def __init__(self,
               structure_monitor,
               mon_lib_srv,
               map_cc_all_threshold = 0.8,
               map_cc_sidechain_threshold=0.8,
               use_slope            = True,
               use_torsion_search   = True,
               use_rotamer_iterator = True,
               torsion_search_all_start = 0,
               torsion_search_all_stop  = 360,
               torsion_search_all_step  = 1,
               torsion_search_local_start = -50,
               torsion_search_local_stop  = 50,
               torsion_search_local_step  = 5,
               backbone_sample            = True,
               log                        = None):
    adopt_init_args(self, locals())
    if(self.log is None): self.log = sys.stdout
    self.special_position_indices = \
      self.structure_monitor.xray_structure.special_position_indices()
    self.rotamer_manager = RotamerEval()
    self.atom_radius_to_negate_map_within = None
    if(self.structure_monitor.target_map_object.miller_array.d_min()>3.5):
      self.atom_radius_to_negate_map_within = 2.0
    else:
      self.atom_radius_to_negate_map_within = 1.5
    #
    self.selection_good = self.structure_monitor.map_cc_per_atom > 0.8
    self.iselection_backbone=flex.size_t()
    for r in self.structure_monitor.residue_monitors:
      if(r.selection_backbone is not None):
        self.iselection_backbone.extend(r.selection_backbone)
    self.loop_over_residues()
    clash_list = self.structure_monitor.find_sidechain_clashes()
    if(clash_list.size()>0):
      self.loop_over_residues(iselection = clash_list, use_clash_filter=True)
      clash_list = self.structure_monitor.find_sidechain_clashes()
    if(clash_list.size()>0):
      self.loop_over_residues(iselection = clash_list, use_clash_filter=True,
        use_torsion_search=True, use_rotamer_iterator=False)
    ####
    # TODO: attemp to search valid rotamer outliers
    #for r in self.structure_monitor.residue_monitors:
    #  print r.id_str, r.map_cc_all, r.map_cc_backbone, r.map_cc_sidechain
    #  if(r.map_cc_backbone>0.9 and r.map_cc_sidechain<0.5):
    #    print r.id_str
    ####

  def on_special_position(self, sm_residue):
    if(self.special_position_indices.size()==0): return False
    for i_seq in sm_residue.selection_all:
      if(i_seq in self.special_position_indices): return True
    return False

  def loop_over_residues(self,
                         iselection           = None,
                         use_clash_filter     = False,
                         debug                = False,
                         use_slope            = None,
                         use_torsion_search   = None,
                         use_rotamer_iterator = None):
    if(use_slope is None): use_slope = self.use_slope
    if(use_torsion_search is None): use_torsion_search = self.use_torsion_search
    if(use_rotamer_iterator is None): use_rotamer_iterator = self.use_rotamer_iterator
    sm = self.structure_monitor
    xrs = sm.xray_structure.deep_copy_scatterers()
    sites_cart = sm.xray_structure.sites_cart()
    get_class = iotbx.pdb.common_residue_names_get_class
    for i_res, r in enumerate(sm.residue_monitors):
      if(iselection is not None and not i_res in iselection): continue
      if(get_class(r.residue.resname) != "common_amino_acid"): continue
      go = r.rotamer_status=="OUTLIER" or \
           r.map_cc_all < self.map_cc_all_threshold or \
           r.map_cc_sidechain < self.map_cc_sidechain_threshold or \
           use_clash_filter
      if(go):
        if(self.on_special_position(sm_residue=r)): continue
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
        #t0=time.time()
        negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
          xray_structure          = xrs,
          selection_within_radius = 5, # XXX make residue dependent !!!!
          iselection              = r.residue.atoms().extract_i_seq(),
          selection_good          = self.selection_good,
          iselection_backbone     = self.iselection_backbone,
          iselection_n_external   = iselection_n_external,
          iselection_c_external   = iselection_c_external)
        #print "sel: %6.4f"%(time.time()-t0)
        target_map_ = mmtbx.refinement.real_space.\
          negate_map_around_selected_atoms_except_selected_atoms(
            xray_structure          = xrs,
            map_data                = sm.target_map_object.data.deep_copy(),
            negate_selection        = negate_selection,
            atom_radius             = self.atom_radius_to_negate_map_within)
        #print "neg: %6.4f"%(time.time()-t0)
        if(debug): # DEBUG
          iotbx.ccp4_map.write_ccp4_map(
            file_name      = "AMap.map",
            unit_cell      = sm.xray_structure.unit_cell(),
            space_group    = sm.xray_structure.space_group(),
            gridding_first = (0,0,0),
            gridding_last  = tuple(sm.target_map_object.crystal_gridding.n_real()),
            map_data       = target_map_,
            labels         = flex.std_string(["DEBUG"]))
        mmtbx.refinement.real_space.fit_residue.manager(
          target_map           = target_map_,
          mon_lib_srv          = self.mon_lib_srv,
          special_position_settings = xrs.special_position_settings(),
          residue              = r.residue,
          sites_cart_all       = sites_cart,
          rotamer_manager      = self.rotamer_manager,
          use_clash_filter     = use_clash_filter,
          debug                = debug,
          use_slope            = use_slope,
          use_torsion_search   = use_torsion_search,
          use_rotamer_iterator = use_rotamer_iterator,
          torsion_search_all_start = self.torsion_search_all_start,
          torsion_search_all_stop  = self.torsion_search_all_stop ,
          torsion_search_all_step  = self.torsion_search_all_step ,
          torsion_search_local_start = self.torsion_search_local_start,
          torsion_search_local_stop  = self.torsion_search_local_stop,
          torsion_search_local_step  = self.torsion_search_local_step,
          backbone_sample            = self.backbone_sample,
          log                        = self.log)
        sites_cart = sites_cart.set_selected(r.residue.atoms().extract_i_seq(),
          r.residue.atoms().extract_xyz())
        xrs.set_sites_cart(sites_cart)
        #print "ref: %6.4f"%(time.time()-t0)
    sm.update(xray_structure = xrs, accept_as_is=True)
