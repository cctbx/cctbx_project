from __future__ import division
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
from scitbx import matrix
from mmtbx.utils import rotatable_bonds
from libtbx import adopt_init_args
import iotbx.pdb
import iotbx.ccp4_map
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import time

class manager(object):
  def __init__(self,
               structure_monitor,
               mon_lib_srv,
               refine_if_map_cc_all_less_than=1.0):
    adopt_init_args(self, locals())
    self.loop_over_residues()

  def loop_over_residues(self):
    sm = self.structure_monitor
    select_good = self.structure_monitor.map_cc_per_atom > 0.8
    sites_cart = sm.xray_structure.sites_cart()
    get_class = iotbx.pdb.common_residue_names_get_class
    for r in sm.residue_monitors:
      if(get_class(r.residue.resname) == "common_amino_acid"):
        #and str(r.residue.resseq).strip() == "39"):
        #print r.residue.resname, r.residue.resseq, "<"*10
        target_map_ = mmtbx.refinement.real_space.\
          negate_map_around_selected_atoms_except_selected_atoms(
            xray_structure          = sm.xray_structure,
            map_data                = sm.target_map_object.data,
            iselection              = r.residue.atoms().extract_i_seq(),
            selection_within_radius = 5,
            select_good             = select_good,
            atom_radius             = 1.5)
        if(0): # DEBUG
          iotbx.ccp4_map.write_ccp4_map(
            file_name      = "AMap.map",
            unit_cell      = sm.xray_structure.unit_cell(),
            space_group    = sm.xray_structure.space_group(),
            gridding_first = (0,0,0),
            gridding_last  = tuple(sm.target_map_object.crystal_gridding.n_real()),
            map_data       = target_map_,
            labels         = flex.std_string(["???"]))
        mmtbx.refinement.real_space.fit_residue.manager(
          target_map           = target_map_,
          mon_lib_srv          = self.mon_lib_srv,
          unit_cell            = sm.xray_structure.unit_cell(),
          residue              = r.residue,
          use_slope            = True,
          use_torsion_search   = True,
          use_rotamer_iterator = True)
        sites_cart.set_selected(r.residue.atoms().extract_i_seq(),
          r.residue.atoms().extract_xyz())
    sm.xray_structure = sm.xray_structure.replace_sites_cart(sites_cart)
    sm.pdb_hierarchy.adopt_xray_structure(sm.xray_structure)
    sm.update(xray_structure = sm.xray_structure)
