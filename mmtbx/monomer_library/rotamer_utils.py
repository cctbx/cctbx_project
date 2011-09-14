import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.math
from scitbx.array_family import flex
from libtbx.str_utils import show_string
from libtbx.utils import sequence_index_dict
import math

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
fine_sampling = False
  .type=bool
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

__rotamer_info_master_phil = None
def rotamer_info_master_phil():
  global __rotamer_info_master_phil
  if (__rotamer_info_master_phil is None):
    import libtbx.phil
    __rotamer_info_master_phil = libtbx.phil.parse(
      input_string=rotamer_info_master_phil_str)
  return __rotamer_info_master_phil

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
      constrain_dihedrals_with_sigma_less_than_or_equal_to,
      external_edge_list=None,
      external_clusters=None,
      return_none_if_unexpected_degrees_of_freedom=False,
      tree_root_atom_names=set(["N", "CA", "C", "O"]),
      terminal_backbone_atom_names=set(["OXT", "HXT", "H1", "H2", "H3"])):
  assert len(mon_lib_atom_names) == len(input_atom_names)
  assert len(sites_cart) == len(input_atom_names)
  atom_indices = sequence_index_dict(seq=mon_lib_atom_names)
  fixed_vertices = []
  for i,atom_name in enumerate(mon_lib_atom_names):
    if (atom_name in tree_root_atom_names):
      fixed_vertices.append(i)
  if(not (len(fixed_vertices) == len(tree_root_atom_names))): return None
  for i,atom_name in enumerate(mon_lib_atom_names):
    if (atom_name in terminal_backbone_atom_names):
      fixed_vertices.append(i)
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
  if (external_edge_list is not None):
    for eed in external_edge_list:
      edge_list.append(tuple(sorted(eed)))
  unused = bonds_to_omit.difference(bonds_omitted)
  if (len(unused) != 0):
    raise RuntimeError(
      "tree_generation_without_bond does not match any bonds: %s"
        % str(unused))
  if (constrain_dihedrals_with_sigma_less_than_or_equal_to is not None):
    if (external_clusters is None):
      external_clusters = []
    else:
      external_clusters = list(external_clusters)
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
  if (joint_dofs[0] != 0 or not joint_dofs[1:].all_eq(1)):
    if (return_none_if_unexpected_degrees_of_freedom):
      return None
    msg = ["Unexpected degrees of freedom:"]
    for dof,cluster in zip(joint_dofs,tardy_tree.cluster_manager.clusters):
      msg.append("  %s: %s" % (dof, [input_atom_names[i] for i in cluster]))
    raise RuntimeError("\n".join(msg))
  return tardy_model

def build_rotamer_tor_atom_ids_by_tor_id(comp_comp_id, rotamer_info):
  comp_tor_by_id = {}
  for tor in comp_comp_id.tor_list:
    assert tor.id not in comp_tor_by_id
    comp_tor_by_id[tor.id] = tor
  rotmer_info_tor_ids = set(rotamer_info.tor_ids)
  rotamer_tor_by_id = {}
  for tor_atom_ids in rotamer_info.tor_atom_ids:
    assert len(tor_atom_ids) == 5
    tor_id = tor_atom_ids[0]
    assert tor_id in rotmer_info_tor_ids
    assert tor_id not in rotamer_tor_by_id
    rotamer_tor_by_id[tor_id] = tuple(tor_atom_ids[1:])
  result = {}
  for tor_id in rotamer_info.tor_ids:
    atom_ids = rotamer_tor_by_id.get(tor_id)
    if (atom_ids is not None):
      result[tor_id] = atom_ids
    else:
      comp_tor = comp_tor_by_id.get(tor_id)
      if (comp_tor is not None):
        result[tor_id] = comp_tor.atom_ids()
      else:
        raise RuntimeError(
          "rotamer_info.tor_id %s is unknown." % show_string(tor_id))
  return result

def build_i_q_packed_by_tor_id(rotamer_tor_atom_ids_by_tor_id, tardy_model):
  tor_id_by_rotatable_bond_atom_names = {}
  for tor_id,atom_ids in rotamer_tor_atom_ids_by_tor_id.items():
    atom_names = tuple(sorted(atom_ids[1:3]))
    assert atom_names not in tor_id_by_rotatable_bond_atom_names
    tor_id_by_rotatable_bond_atom_names[atom_names] = tor_id
  result = {}
  number_of_trees = 0
  for i_body,he in enumerate(
                     tardy_model.tardy_tree.cluster_manager.hinge_edges):
    if (he[0] == -1):
      number_of_trees += 1
      continue
    hinge_atom_names = [tardy_model.labels[i].strip() for i in he]
    atom_names = tuple(sorted(hinge_atom_names))
    tor_id = tor_id_by_rotatable_bond_atom_names.get(atom_names)
    if (tor_id is None):
      raise RuntimeError(
        "rotatable bond atoms %s - %s (as defined by tardy_tree):"
        " no match in rotamer_info.tor_ids" % tuple(hinge_atom_names))
    result[tor_id] = i_body - 1
  assert number_of_trees == 1
  return result

def build_angle_start_by_tor_id(
      mon_lib_atom_names,
      sites_cart,
      rotamer_tor_atom_ids_by_tor_id,
      i_q_packed_by_tor_id):
  result = {}
  atom_indices = sequence_index_dict(seq=mon_lib_atom_names)
  for tor_id in i_q_packed_by_tor_id.keys():
    d_sites = []
    for atom_id in rotamer_tor_atom_ids_by_tor_id[tor_id]:
      i = atom_indices.get(atom_id)
      if (i is None):
        return (atom_id, tor_id)
      d_sites.append(sites_cart[i])
    dihe = scitbx.math.dihedral_angle(sites=d_sites, deg=True)
    assert dihe is not None
    assert tor_id not in result
    result[tor_id] = dihe
  return result

class rotamer_iterator(object):

  def __init__(O, comp_comp_id, atom_names, sites_cart, fine_sampling=False):
    assert sites_cart.size() == len(atom_names)
    O.problem_message = None
    O.rotamer_info = comp_comp_id.rotamer_info()
    if (O.rotamer_info is None):
      return
    if fine_sampling == True:
      O.rotamer_info.fine_sampling = True
    resname = comp_comp_id.chem_comp.id
    import iotbx.pdb.atom_name_interpretation
    matched_atom_names = iotbx.pdb.atom_name_interpretation.interpreters[
      resname].match_atom_names(atom_names=atom_names)
    names = matched_atom_names.unexpected
    if (len(names) != 0):
      O.problem_message = "resname=%s: unexpected atoms: %s" % (
        resname, " ".join(sorted(names)))
      return
    names = matched_atom_names.missing_atom_names(ignore_hydrogen=True)
    if (len(names) != 0):
      O.problem_message = "resname=%s: missing atoms: %s" % (
        resname, " ".join(sorted(names)))
      return
    O.mon_lib_atom_names = matched_atom_names.mon_lib_names()
    if (O.rotamer_info.atom_ids_not_handled is not None):
      atom_ids_not_handled = set(O.rotamer_info.atom_ids_not_handled)
      not_handled = []
      for atom_name, mon_lib_atom_name in zip(atom_names, O.mon_lib_atom_names):
        if (mon_lib_atom_name in atom_ids_not_handled):
          not_handled.append(atom_name.strip())
      if (len(not_handled) != 0):
        O.problem_message = \
          "%s: rotamer_info does not handle these atoms: %s" % (
            resname, " ".join(not_handled))
        return
    O.bonds_to_omit = extract_bonds_to_omit(rotamer_info=O.rotamer_info)
    O.tardy_model = tardy_model(
      comp_comp_id=comp_comp_id,
      input_atom_names=atom_names,
      mon_lib_atom_names=O.mon_lib_atom_names,
      sites_cart=sites_cart,
      bonds_to_omit=O.bonds_to_omit,
      constrain_dihedrals_with_sigma_less_than_or_equal_to
        =O.rotamer_info.constrain_dihedrals_with_sigma_less_than_or_equal_to)
    O.rotamer_tor_atom_ids_by_tor_id = build_rotamer_tor_atom_ids_by_tor_id(
      comp_comp_id=comp_comp_id,
      rotamer_info=O.rotamer_info)
    O.i_q_packed_by_tor_id = build_i_q_packed_by_tor_id(
      rotamer_tor_atom_ids_by_tor_id=O.rotamer_tor_atom_ids_by_tor_id,
      tardy_model=O.tardy_model)
    build_result = build_angle_start_by_tor_id(
      mon_lib_atom_names=O.mon_lib_atom_names,
      sites_cart=sites_cart,
      rotamer_tor_atom_ids_by_tor_id=O.rotamer_tor_atom_ids_by_tor_id,
      i_q_packed_by_tor_id=O.i_q_packed_by_tor_id)
    if (isinstance(build_result, dict)):
      O.angle_start_by_tor_id = build_result
    else:
      O.problem_message = 'resname=%s: missing atom "%s" for tor_id "%s"' % (
        (resname,) + build_result)
      return
    O.reset()

  def reset(O):
    O.__iterates = iter(O.rotamer_info.rotamer)

  def __iter__(O):
    return O

  def next(O):
    rotamer = O.__iterates.next()
    if O.rotamer_info.fine_sampling == False:
      while(rotamer.frequency_annotation == "for more uniform sampling"):
        rotamer = O.__iterates.next()
    q_packed_work = flex.double(O.tardy_model.q_packed_size, 0)
    for tor_id,angle in zip(O.rotamer_info.tor_ids, rotamer.angles):
      i_q_packed = O.i_q_packed_by_tor_id.get(tor_id)
      if (i_q_packed is not None and angle is not None):
        q_packed_work[i_q_packed] = math.radians(
          angle - O.angle_start_by_tor_id[tor_id])
    O.tardy_model.unpack_q(q_packed=q_packed_work)
    return rotamer, O.tardy_model.sites_moved()
