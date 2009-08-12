import mmtbx.monomer_library.server
import iotbx.pdb.amino_acid_codes
import scitbx.graph.tardy_tree
import sys

def process(mon_lib_srv, resname):
  print "resname:", resname
  comp = mon_lib_srv.get_comp_comp_id_direct(comp_id=resname)
  backbone_atom_names = set(["N", "CA", "C", "O"])
  fixed_vertices = []
  atom_indices = {}
  for i,atom in enumerate(comp.atom_list):
    atom_id = atom.atom_id
    assert atom_id not in atom_indices
    atom_indices[atom_id] = i
    if (atom_id in backbone_atom_names):
      fixed_vertices.append(i)
  assert len(fixed_vertices) == len(backbone_atom_names)
  edge_list = []
  for bond in comp.bond_list:
    edge_list.append(tuple(sorted([atom_indices[atom_id]
      for atom_id in bond.atom_id_1, bond.atom_id_2])))
  tor_dict = {}
  for tor in comp.tor_list:
    atom_names = tuple(sorted([tor.atom_id_2, tor.atom_id_3]))
    tor_dict.setdefault(atom_names, []).append(tor)
    print tor.id, ", ".join([
      tor.atom_id_1, tor.atom_id_2, tor.atom_id_3, tor.atom_id_4]), \
      tor.value_angle_esd
  for tors in tor_dict.values():
    if (len(tors) != 1):
      print "Warning: multiple tors:", ", ".join([tor.id for tor in tors])
  tardy_tree = scitbx.graph.tardy_tree.construct(
    n_vertices=len(comp.atom_list),
    edge_list=edge_list,
    fixed_vertex_lists=[fixed_vertices]).build_tree()
  assert len(tardy_tree.cluster_manager.loop_edges) == 0
  tor_hinge_matches = set()
  number_of_trees = 0
  for he in tardy_tree.cluster_manager.hinge_edges:
    if (he[0] == -1):
      number_of_trees += 1
      continue
    hinge_atom_names = [comp.atom_list[i].atom_id for i in he]
    atom_names = tuple(sorted(hinge_atom_names))
    tors = tor_dict.get(atom_names)
    if (tors is None):
      s = "Warning: no tor"
    else:
      for tor in tors:
        tor_hinge_matches.add(tor.id)
      s = ", ".join([tor.id for tor in tors])
      if (len(tors) != 1):
        s = "Warning: multiple tors: " + s
    print "hinge edge:", ", ".join(hinge_atom_names), s
  assert number_of_trees == 1
  non_const_tor_ids = set()
  for tor in comp.tor_list:
    if (tor.value_angle_esd == 0):
      assert tor.id.startswith("CONST_")
    else:
      non_const_tor_ids.add(tor.id)
  tors_not_hinge = non_const_tor_ids.difference(tor_hinge_matches)
  if (len(tors_not_hinge) != 0):
    print "tors_not_hinge:", ", ".join(sorted(tors_not_hinge))
  print

def run(args):
  assert len(args) == 0
  mon_lib_srv = mmtbx.monomer_library.server.server()
  amino_acid_resnames = sorted(
    iotbx.pdb.amino_acid_codes.three_letter_given_one_letter.values())
  for resname in amino_acid_resnames:
    process(mon_lib_srv=mon_lib_srv, resname=resname)
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
