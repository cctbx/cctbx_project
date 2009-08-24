import scitbx.rigid_body
import scitbx.graph.tardy_tree
from scitbx.array_family import flex
from libtbx.utils import sequence_index_dict

rotamer_info_master_phil_str = """\
tor_ids = None
  .type = strings
tor_atom_ids = None
  .type = strings
  .multiple = True
atom_ids_not_handled = None
  .type = strings
tree_generation_without_bond = None
  .type = strings
  .multiple = True
constrain_dihedrals_with_sigma_less_than_or_equal_to = 10
  .type=float
rotamer
  .multiple = True
{
 id = None
   .type = str
 frequency = None
   .type = float
 frequency_annotation = None
   .type = str
 angles = None
   .type = floats(allow_none_elements=True)
}
"""

def extract_bonds_to_omit(rotamer_info):
  result = set()
  if (rotamer_info is not None):
    for bond in rotamer_info.tree_generation_without_bond:
      assert len(bond) == 2
      bond = tuple(bond)
      if (bond in result):
        raise RuntimeError(
          "Duplicate tree_generation_without_bond definition: %s" % str(bond))
      result.add(bond)
  return result

def tardy_model(
      comp_comp_id,
      input_atom_names,
      mon_lib_atom_names,
      sites_cart,
      bonds_to_omit,
      constrain_dihedrals_with_sigma_less_than_or_equal_to):
  assert len(mon_lib_atom_names) == len(input_atom_names)
  assert len(sites_cart) == len(input_atom_names)
  atom_indices = sequence_index_dict(seq=mon_lib_atom_names)
  tree_root_atom_names = set(["N", "CA", "C", "O"])
  fixed_vertices = []
  for i,atom_name in enumerate(mon_lib_atom_names):
    if (atom_name in tree_root_atom_names):
      fixed_vertices.append(i)
  assert len(fixed_vertices) == len(tree_root_atom_names)
  bonds_omitted = set()
  edge_list = []
  for bond in comp_comp_id.bond_list:
    bond_atom_ids = bond.atom_ids()
    if (bond_atom_ids in bonds_to_omit):
      bonds_omitted.add(bond_atom_ids)
    else:
      ai = [atom_indices.get(atom_id) for atom_id in bond_atom_ids]
      if (ai.count(None) == 0):
        edge_list.append(tuple(sorted(ai)))
  unused = bonds_to_omit.difference(bonds_omitted)
  if (len(unused) != 0):
    raise RuntimeError(
      "tree_generation_without_bond does not match any bonds: %s"
        % str(unused))
  external_clusters = []
  if (constrain_dihedrals_with_sigma_less_than_or_equal_to is not None):
    for tor in comp_comp_id.tor_list:
      if (   tor.value_angle_esd
          <= constrain_dihedrals_with_sigma_less_than_or_equal_to):
        ai = [atom_indices.get(atom_id) for atom_id in tor.atom_ids()]
        if (ai.count(None) == 0):
          external_clusters.append(sorted(ai))
  for plane in comp_comp_id.get_planes():
    ai = []
    for atom_id in plane.plane_atoms:
      i = atom_indices.get(atom_id)
      if (i is not None):
        ai.append(i)
    external_clusters.append(sorted(ai))
  tardy_tree = scitbx.graph.tardy_tree.construct(
    n_vertices=len(mon_lib_atom_names),
    edge_list=edge_list,
    external_clusters=external_clusters,
    fixed_vertex_lists=[fixed_vertices]).build_tree()
  assert len(tardy_tree.cluster_manager.loop_edges) == 0
  tardy_model = scitbx.rigid_body.tardy_model(
    labels=input_atom_names,
    sites=sites_cart,
    masses=flex.double(len(input_atom_names), 1),
    tardy_tree=tardy_tree,
    potential_obj=None)
  joint_dofs = tardy_model.degrees_of_freedom_each_joint()
  assert joint_dofs[0] == 0
  assert joint_dofs[1:].all_eq(1)
  return tardy_model
