import iotbx.pdb
from libtbx.str_utils import split_keeping_spaces
import sys

trans_dict = {}
for k,v in iotbx.pdb.rna_dna_atom_names_reference_to_mon_lib_translation_dict \
             .items():
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
    if (line == "loop_"):
      print line
      return
    yield line

def rename_atom(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 10
    trans_field(flds, 3)
    print "".join(flds)

def rename_tree(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 10
    for i in [3, 5, 7, 9]:
      if (flds[i] not in ["n/a", "START", "ADD", "END", "."]):
        trans_field(flds, i)
    print "".join(flds)

def rename_bond(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 12
    for i in [3, 5]:
      trans_field(flds, i)
    print "".join(flds)

def rename_angle(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 12
    for i in [3, 5, 7]:
      trans_field(flds, i)
    print "".join(flds)

def rename_tor(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 18
    for i in [5, 7, 9, 11]:
      trans_field(flds, i)
    print "".join(flds)

def rename_chir(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 14
    for i in [5, 7, 9, 11]:
      trans_field(flds, i)
    print "".join(flds)

def rename_plan(lines):
  for line in iter_until_loop(lines):
    flds = split_keeping_spaces(line)
    assert len(flds) == 8
    trans_field(flds, 5)
    print "".join(flds)

def run(args):
  assert len(args) == 1
  lines = iter(open(args[0]).read().splitlines())
  for line in lines:
    print line
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

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
