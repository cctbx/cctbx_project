from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import iotbx.pdb
import iotbx.ccp4_map
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue

class manager(object):
  def __init__(self,
               structure_monitor,
               mon_lib_srv,
               refine_if_map_cc_all_less_than=1.0):
    adopt_init_args(self, locals())
    self.atom_radius_to_negate_map_within = None
    if(self.structure_monitor.target_map_object.miller_array.d_min()>3.5):
      self.atom_radius_to_negate_map_within = 2.0
    else:
      self.atom_radius_to_negate_map_within = 1.5
    #
    self.selection_good = self.structure_monitor.map_cc_per_atom > 0.8
    self.iselection_backbone=flex.size_t()
    for r in self.structure_monitor.residue_monitors:
      self.iselection_backbone.extend(r.selection_backbone)
    self.loop_over_residues()
    clash_list = self.structure_monitor.find_sidechain_clashes()
    self.loop_over_residues(iselection = clash_list, use_clash_filter=True)

  def loop_over_residues(self, iselection=None, use_clash_filter=False,
                         debug=False):
    sm = self.structure_monitor
    sites_cart = sm.xray_structure.sites_cart()
    get_class = iotbx.pdb.common_residue_names_get_class
    for i_res, r in enumerate(sm.residue_monitors):
      if(iselection is not None and not i_res in iselection): continue
      if(get_class(r.residue.resname) == "common_amino_acid"):
        #and str(r.residue.resseq).strip() in ["49"]):
        #print
        #print r.residue.resname, r.residue.resseq, "<"*40
        iselection_n_external=None
        iselection_c_external=None
        if(i_res!=0 and i_res!=len(sm.residue_monitors)-1):
          iselection_c_external = sm.residue_monitors[i_res-1].selection_c
          iselection_n_external = sm.residue_monitors[i_res+1].selection_n
        elif(i_res==0 and len(sm.residue_monitors)>1):
          iselection_n_external = sm.residue_monitors[i_res+1].selection_n
        elif(i_res==len(sm.residue_monitors)-1 and len(sm.residue_monitors)>1):
          iselection_c_external = sm.residue_monitors[i_res-1].selection_c
        negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
          xray_structure          = sm.xray_structure,
          selection_within_radius = 5, # XXX make residue dependent !!!!
          iselection              = r.residue.atoms().extract_i_seq(),
          selection_good          = self.selection_good,
          iselection_backbone     = self.iselection_backbone,
          iselection_n_external   = iselection_n_external,
          iselection_c_external   = iselection_c_external)
        target_map_ = mmtbx.refinement.real_space.\
          negate_map_around_selected_atoms_except_selected_atoms(
            xray_structure          = sm.xray_structure,
            map_data                = sm.target_map_object.data.deep_copy(),
            negate_selection        = negate_selection,
            atom_radius             = self.atom_radius_to_negate_map_within)
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
          special_position_settings = sm.xray_structure.special_position_settings(),
          residue              = r.residue,
          sites_cart_all       = sites_cart,
          use_clash_filter     = use_clash_filter,
          debug                = debug,# XXX
          use_slope            = True,
          use_torsion_search   = True,
          use_rotamer_iterator = True)
        sites_cart.set_selected(r.residue.atoms().extract_i_seq(),
          r.residue.atoms().extract_xyz())
        sm.xray_structure = sm.xray_structure.replace_sites_cart(sites_cart)
    sm.pdb_hierarchy.adopt_xray_structure(sm.xray_structure)
    sm.update(xray_structure = sm.xray_structure)
