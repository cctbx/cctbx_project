from __future__ import absolute_import, division, print_function
import iotbx.pdb
from libtbx.str_utils import split_keeping_spaces
import sys
import six

trans_dict = {}
for k,v in six.iteritems(iotbx.pdb.rna_dna_atom_names_reference_to_mon_lib_translation_dict):
  trans_dict[k.strip()] = v
trans_dict["H2'"] = "H2*"

def trans_field(flds, i):
  v3 = flds[i]
  v2 = trans_dict[v3]
  flds[i] = v2
  if (i+1 < len(flds)):
    l = len(flds[i+1])
    d = len(v2) - len(v3)
    assert l > d
    flds[i+1] = " " * (l-d)

def iter_until_loop(lines):
  for line in lines:
    if (   line.startswith("#")
        or line == "loop_"):
      print(line)
      return
    yield line

def rename_generic(lines, len_flds, i_list):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == len_flds
    for i in i_list:
      trans_field(flds, i)
    print("".join(flds))

def rename_atom(lines):
  rename_generic(lines, 10, [3])

def rename_tree(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 10
    for i in [3, 5, 7, 9]:
      if (flds[i] not in ["n/a", "START", "ADD", "END", "."]):
        trans_field(flds, i)
    print("".join(flds))

def rename_bond(lines):
  rename_generic(lines, 12, [3, 5])

def rename_angle(lines):
  rename_generic(lines, 12, [3, 5, 7])

def rename_tor(lines):
  rename_generic(lines, 18, [5, 7, 9, 11])

def rename_chir(lines):
  rename_generic(lines, 14, [5, 7, 9, 11])

def rename_plan(lines):
  rename_generic(lines, 8, [5])

def rename_link_bond(lines):
  rename_generic(lines, 16, [5, 9])

def rename_link_angle(lines):
  rename_generic(lines, 18, [5, 9, 13])

def rename_link_tor(lines):
  rename_generic(lines, 26, [7, 11, 15, 19])

def run(args):
  assert len(args) == 1
  lines = iter(open(args[0]).read().splitlines())
  for line in lines:
    print(line)
    if (line == "_chem_comp_atom.partial_charge"):
      rename_atom(lines)
    elif (line == "_chem_comp_tree.connect_type"):
      rename_tree(lines)
    elif (line == "_chem_comp_bond.value_dist_esd"):
      rename_bond(lines)
    elif (line == "_chem_comp_angle.value_angle_esd"):
      rename_angle(lines)
    elif (line == "_chem_comp_tor.period"):
      rename_tor(lines)
    elif (line == "_chem_comp_chir.volume_sign"):
      rename_chir(lines)
    elif (line == "_chem_comp_plane_atom.dist_esd"):
      rename_plan(lines)
    #
    elif (line == "_chem_link_bond.value_dist_esd"):
      rename_link_bond(lines)
    elif (line == "_chem_link_angle.value_angle_esd"):
      rename_link_angle(lines)
    elif (line == "_chem_link_tor.period"):
      rename_link_tor(lines)
    elif (line == "_chem_link_chir.volume_sign"):
      raise RuntimeError("Not implemented.")
    elif (line == "_chem_link_plane.dist_esd"):
      raise RuntimeError("Not implemented.")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
