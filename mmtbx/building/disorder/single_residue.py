
# TODO avoid re-processing structure each time
# TODO tests

"""
Prototype for building alternate conformations into difference density.
The actual method is a variant of one that Pavel suggested, in combination
with the procedure in mmtbx.building.extend_sidechains: first, the
backbone atoms for a residue and its neighbors are refined into the
target map using minimization and/or annealing, then the sidechain is
replaced (using an idealized copy) and its placement optimized by a grid
search that also allows for backbone flexibility.
"""

from __future__ import division
from libtbx import adopt_init_args

class rebuild_residue (object) :
  """
  Callable wrapper class for rebuilding a single residue at a time.  This is
  not necessarily limited to modeling disorder, but it has been specifically
  designed to fit to a difference map in the presence of backbone and sidechain
  shifts.  Unlike some of the other tools, this method completely removes all
  sidechains within a sliding window before doing anything else.

  Only the target residue is returned; splitting of adjacent residues will be
  essential in most cases but is not handled here.
  """
  def __init__ (self,
      target_map,
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      d_min) :
    adopt_init_args(self, locals())
    from mmtbx.monomer_library import idealized_aa
    from mmtbx.rotamer import rotamer_eval
    import mmtbx.monomer_library.server
    self.ideal_dict = idealized_aa.residue_dict()
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.rotamer_manager = rotamer_eval.RotamerEval()

  def __call__ (self, atom_group, log, window_size=2, anneal=False) :
    from mmtbx.building import extend_sidechains
    from mmtbx import building
    from scitbx.array_family import flex
    assert (atom_group is not None)
    pdb_hierarchy = self.pdb_hierarchy.deep_copy()
    xray_structure = self.xray_structure.deep_copy_scatterers()
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    isel = building.extract_iselection([atom_group])
    atom_group = pdb_atoms[isel[0]].parent()
    residue_group = atom_group.parent()
    assert (len(residue_group.atom_groups()) == 1)
    chain = residue_group.parent()
    all_residues = chain.residue_groups()
    sel_residues = []
    for i_res, other_rg in enumerate(all_residues) :
      if (other_rg == residue_group) :
        for j_res in range(-window_size, window_size+1) :
          k_res = i_res + j_res
          if (k_res >= 0) :
            sel_residues.append(all_residues[k_res])
    assert len(sel_residues) > 0
    print >> log, "  %d residues extracted" % len(sel_residues)
    print >> log, "  removing sidechain atoms..."
    building.remove_sidechain_atoms(sel_residues)
    pdb_atoms = pdb_hierarchy.atoms()
    all_mc_sel = pdb_atoms.extract_i_seq()
    xrs_mc = xray_structure.select(all_mc_sel)
    pdb_atoms.reset_i_seq()
    window_mc_sel = building.extract_iselection(sel_residues)
    selection = flex.bool(pdb_atoms.size(), False).set_selected(window_mc_sel,
      True)
    restraints_manager = self.geometry_restraints_manager.select(all_mc_sel)
    box = building.box_build_refine_base(
      xray_structure=xrs_mc,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      geometry_restraints_manager=restraints_manager.geometry,
      d_min=self.d_min,
      out=log,
      debug=True)
    box.restrain_atoms(
      selection=box.others_in_box,
      reference_sigma=0.1)
    box.real_space_refine(selection=box.selection_in_box)
    if (anneal) : # TODO
      raise NotImplementedError()
    sites_new = box.update_original_coordinates()
    pdb_atoms.set_xyz(sites_new)
    # extend and replace existing residue
    new_atom_group = extend_sidechains.extend_residue(
      residue=atom_group,
      ideal_dict=self.ideal_dict,
      hydrogens=False,
      mon_lib_srv=self.mon_lib_srv,
      match_conformation=True)
    rg = atom_group.parent()
    rg.remove_atom_group(atom_group)
    rg.append_atom_group(new_atom_group)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    # get new box around this residue
    residue_sel = building.extract_iselection([ new_atom_group ])
    selection = flex.bool(pdb_atoms.size(), False).set_selected(residue_sel,
      True)
    xray_structure = pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=self.xray_structure)
    # FIXME this is horrendously inefficient
    box = building.box_build_refine_base(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      d_min=self.d_min,
      out=log,
      debug=True)
    box.fit_residue_in_box()
    sites_new = box.update_original_coordinates()
    pdb_hierarchy.atoms().set_xyz(sites_new)
    return new_atom_group.detached_copy()
