from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import server
from libtbx import dict_with_default_0
import string
import os

def detect_unknown_type_energy(comp_id, comp_comp_id):
  n_unknown_type_energy = 0
  for atom in comp_comp_id.atom_list:
    if (atom.type_energy is None or len(atom.type_energy) == 0):
      n_unknown_type_energy += 1
  status = "ok"
  if (n_unknown_type_energy != 0):
    print("unknown type_energy:", comp_id, n_unknown_type_energy, \
      len(comp_comp_id.atom_list), end=' ')
    if (n_unknown_type_energy == len(comp_comp_id.atom_list)):
      status = "all"
    else:
      status = "some"
      print("only some unknown", end=' ')
    print()
  return status

def detect_missing_angle_definitions(comp_id, comp_comp_id):
  if (len(comp_comp_id.atom_list) <= 2): return "ok"
  if (len(comp_comp_id.angle_list) > 0): return "ok"
  print("missing bond angle definitions:", comp_id)
  return "all"

def detect_missing_bond_values(comp_id, comp_comp_id):
  n_missing_ideal = 0
  n_missing_esd = 0
  for bond in comp_comp_id.bond_list:
    if (bond.value_dist is None): n_missing_ideal += 1
    if (bond.value_dist_esd is None): n_missing_esd += 1
  if (n_missing_ideal != 0 or n_missing_esd != 0):
    print("missing bond values:", n_missing_ideal, n_missing_esd)
    return "bad"
  return "ok"

def detect_missing_angle_values(comp_id, comp_comp_id):
  n_missing_ideal = 0
  n_missing_esd = 0
  for bond in comp_comp_id.angle_list:
    if (bond.value_angle is None): n_missing_ideal += 1
    if (bond.value_angle_esd is None): n_missing_esd += 1
  if (n_missing_ideal != 0 or n_missing_esd != 0):
    print("missing angle values:", n_missing_ideal, n_missing_esd)
    return "bad"
  return "ok"

def exercise():
  list_cif = server.mon_lib_list_cif()
  srv = server.server(list_cif=list_cif)
  print("srv.root_path:", srv.root_path)
  table_of_contents = []
  n_get_comp_comp_id_successes = 0
  unknown_type_energy_counts = dict_with_default_0()
  missing_angle_definitions_counts = dict_with_default_0()
  missing_bond_values_counts = dict_with_default_0()
  missing_angle_values_counts = dict_with_default_0()
  for first_char in string.ascii_lowercase+string.digits:
    sub_dir = os.path.join(srv.root_path, first_char)
    if (not os.path.isdir(sub_dir)): continue
    for node in os.listdir(sub_dir):
      if (not node.lower().endswith(".cif")): continue
      comp_id = node[:-4]
      if (comp_id.endswith("_EL")): continue
      if (comp_id in ["CON_CON", "PRN_PRN"]):
        comp_id = comp_id[:3]
      if (comp_id.upper() != comp_id):
        print("Mixed case:", os.path.join(first_char, node))
      comp_comp_id = srv.get_comp_comp_id_direct(comp_id=comp_id)
      if (comp_comp_id is None):
        print("Error instantiating comp_comp_id %s (%s)" % (
          comp_id, os.path.join(sub_dir, node)))
      else:
        n_get_comp_comp_id_successes += 1
        table_of_contents.append(
          " ".join([comp_id.upper(), os.path.join(first_char, node)]))
        status = detect_unknown_type_energy(
          comp_id=comp_id, comp_comp_id=comp_comp_id)
        unknown_type_energy_counts[status] += 1
        status = detect_missing_angle_definitions(
          comp_id=comp_id, comp_comp_id=comp_comp_id)
        missing_angle_definitions_counts[status] += 1
        status = detect_missing_bond_values(
          comp_id=comp_id, comp_comp_id=comp_comp_id)
        missing_bond_values_counts[status] += 1
        status = detect_missing_angle_values(
          comp_id=comp_id, comp_comp_id=comp_comp_id)
        missing_angle_values_counts[status] += 1
        if (1 and status != "ok"):
          print('svn rm "%s"' % os.path.join(first_char, node))
  print("number of cif files read successfully:", n_get_comp_comp_id_successes)
  print("unknown type_energy counts:", unknown_type_energy_counts)
  print("missing bond angle definitions counts:", \
    missing_angle_definitions_counts)
  print("missing bond values counts:", missing_bond_values_counts)
  print("missing angle values counts:", missing_angle_values_counts)
  print("writing file table_of_contents")
  with open("table_of_contents", "w") as f:
    f.write("\n".join(table_of_contents)+"\n")

if (__name__ == "__main__"):
  exercise()
