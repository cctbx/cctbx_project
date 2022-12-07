from __future__ import absolute_import, division, print_function
import iotbx.pdb
import sys
import scitbx.graph.tardy_tree
import mmtbx.monomer_library.rotamer_utils

def tardy_model_one_residue(residue, mon_lib_srv, log = None):
  try:
    resname = residue.resname
  except AttributeError:
    resname = residue.unique_resnames()
    assert resname.size()==1
    resname = resname[0]
  # FIXME atom_group objects have no id_str
  def format_atom_group_id_str(residue):
    assert (type(residue).__name__ == "atom_group")
    return '"%3s %s%4s%s"' % (residue.resname, residue.parent().parent().id,
      residue.parent().resseq, residue.parent().icode)
  if(log is None): log = sys.stdout
  comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  # create external clusters to address planes with missing atoms
  external_clusters = []
  for plane in comp_comp_id.get_planes():
    plane_atom_names = []
    for plane_atom in plane.plane_atoms:
      plane_atom_names.append(plane_atom.atom_id)
    plane_atom_i_seqs = []
    if 0: print(plane_atom_names)
    for i_atom, atom in enumerate(residue.atoms()):
      if(atom.name.strip() in plane_atom_names):
        plane_atom_i_seqs.append(i_atom)
    if(len(plane_atom_i_seqs)>0):
      external_clusters.append(plane_atom_i_seqs)
  #
  if 0: print("external_clusters1:",external_clusters)
  residue_atoms = residue.atoms()
  atom_names = residue_atoms.extract_name()
  matched_atom_names = iotbx.pdb.atom_name_interpretation.interpreters[
    resname].match_atom_names(atom_names=atom_names)
  mon_lib_atom_names = matched_atom_names.mon_lib_names()
  if(mon_lib_atom_names.count(None)>1):
    msg="Canot create TARDY: bad atom names in '%s %s':"%(
      resname, residue.resid())
    print(msg, file=log)
    print(list(residue.atoms().extract_name()), file=log)
    return None
  #
  rotamer_info = comp_comp_id.rotamer_info()
  bonds_to_omit = mmtbx.monomer_library.rotamer_utils.extract_bonds_to_omit(
    rotamer_info = rotamer_info)
  #
  external_edge_list = []
  special_cases = {
    "ASP": [(["OD1"],["HD1","1HD","DD1","1DD"]),
            (["OD2"],["HD2","2HD","DD2","2DD"])],
    "GLU": [(["OE1"],["HE1","1HE","DE1","1DE"]),
            (["OE2"],["HE2","2HE","DE2","2DE"])]
  }
  if resname.strip().upper() in special_cases:
    edge = []
    bonded_atom_names = special_cases[residue.resname.strip().upper()]
    for ban in bonded_atom_names:
      for ia1,atom1 in enumerate(residue_atoms):
        if(atom1.name.strip().upper() in ban[0]):
          for ia2,atom2 in enumerate(residue_atoms):
            if(atom2.name.strip().upper() in ban[1]):
              assert ia1 != ia2
              edge = [ia1, ia2]
      if(len(edge) > 0): external_edge_list.append(edge)
  #
  tardy_model = mmtbx.monomer_library.rotamer_utils.tardy_model(
    comp_comp_id       = comp_comp_id,
    input_atom_names   = atom_names,
    mon_lib_atom_names = mon_lib_atom_names,
    sites_cart         = residue_atoms.extract_xyz(),
    bonds_to_omit      = bonds_to_omit,
    constrain_dihedrals_with_sigma_less_than_or_equal_to = None,
    external_edge_list = external_edge_list,
    external_clusters  = external_clusters,
    return_none_if_unexpected_degrees_of_freedom = True)
  if(tardy_model is None):
    mes = "TARDY: cannot create tardy model for: %s (corrupted residue). Skipping it."
    if (hasattr(residue, "id_str")):
      print(mes%residue.id_str(suppress_segid=1)[-12:], file=log)
    else :
      print(mes%format_atom_group_id_str(residue), file=log)
    return None
  joint_dofs = tardy_model.degrees_of_freedom_each_joint()
  if(joint_dofs[0] != 0 or not joint_dofs[1:].all_eq(1)):
    mes = "TARDY error: unexpected degrees of freedom for %s. Skipping it."
    if (hasattr(residue, "id_str")):
      print(mes%residue.id_str(suppress_segid=1)[-12:], file=log)
    else :
      print(mes%format_atom_group_id_str(residue), file=log)
    return None
  return tardy_model

def axes_and_atoms_aa_specific(
      residue, mon_lib_srv,
      remove_clusters_with_all_h=False,
      include_labels=False,
      tardy_model=None,
      log=None):
  get_class = iotbx.pdb.common_residue_names_get_class
  try: # so it can be residue or residue_group
    resname = residue.resname
  except AttributeError:
    resname = residue.unique_resnames()
    assert resname.size()==1
    resname = resname[0]
  if(tardy_model is None):
    if(not (get_class(resname) == "common_amino_acid")): return None
    tardy_model = tardy_model_one_residue(residue = residue,
      mon_lib_srv = mon_lib_srv, log = log)
  if(tardy_model is None):
    if include_labels:
      return None, None
    else:
      return None
  clusters = tardy_model.tardy_tree.cluster_manager.clusters[1:]
  axes = tardy_model.tardy_tree.cluster_manager.hinge_edges[1:]
  assert len(clusters) == len(axes)
  if(len(axes)==0):
    if include_labels:
      return None, None
    else:
      return None
  if 0:
    print("clusters:", clusters)
    print("    axes:", axes)
    print()
  #
  ic = 0
  axes_and_atoms_to_rotate = []
  while ic < len(axes):
    axis = axes[ic]
    cluster = clusters[ic]
    n_branches = 0
    for i_axis_, axis_ in enumerate(axes):
      if(i_axis_ > ic and axis[0] == axis_[0]): n_branches += 1
    if(n_branches == 0):
      atoms_to_rotate = []
      for ci in clusters[ic:]:
        atoms_to_rotate.extend(ci)
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
    else:
      if(not n_branches <= 2):
        print(residue.id_str(), file=log)
        print("  Cannot extract axes_and_atoms_to_rotate:", file=log)
        print("    n_branches <= 2, n_branches:", n_branches, file=log)
        return []
      # branch 1
      icb = ic+1
      atoms_to_rotate = cluster
      while icb < len(axes):
        if(axis[1]==axes[icb][0]): atoms_to_rotate.extend(clusters[icb])
        icb += 1
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
      # branch 2
      axis = axes[ic]
      cluster = clusters[ic]
      icb = ic+1
      atoms_to_rotate = cluster
      while icb < len(axes):
        if(axis[0]==axes[icb][0]): atoms_to_rotate.extend(clusters[icb])
        icb += 1
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
  #
  if 0:
    print(axes_and_atoms_to_rotate)
    print()
  if(remove_clusters_with_all_h):
    tmp = []
    residue_atoms = residue.atoms()
    for axis, cluster in axes_and_atoms_to_rotate:
      count_h = 0
      for c_i in cluster:
        if(residue_atoms[c_i].element.strip().upper() in ["H","D"]): count_h +=1
      if(count_h != len(cluster)):
        tmp.append([axis, cluster])
    axes_and_atoms_to_rotate = tmp
  if include_labels:
    return axes_and_atoms_to_rotate, tardy_model.labels
  else:
    return axes_and_atoms_to_rotate
