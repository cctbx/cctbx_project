"""
Symmetry copies of atoms in pdb input file.
Uses two-character chain ids in output pdb file.
Therefore it only works for space groups with less than
ten symmetry operations.
"""
from __future__ import absolute_import, division, print_function

def run(args):
  assert len(args) == 1, "pdb_file_name"
  import iotbx.pdb
  pdb_inp = iotbx.pdb.input(file_name=args[0])
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  cs = pdb_inp.crystal_symmetry()
  frac_mx = cs.unit_cell().fractionalization_matrix()
  orth_mx = cs.unit_cell().orthogonalization_matrix()
  sites_frac = frac_mx * pdb_inp.atoms().extract_xyz()
  pdb_hierarchies = [pdb_hierarchy]
  for i_sym_op,sym_op in enumerate(cs.space_group()):
    if (sym_op.is_unit_mx()): continue
    print(str(sym_op))
    pdb_hierarchy_copy = pdb_hierarchy.deep_copy()
    pdb_hierarchy_copy.atoms().set_xyz(
      new_xyz=orth_mx * (sym_op.as_rational().as_float() * sites_frac))
    pdb_hierarchies.append(pdb_hierarchy_copy)
  combined_pdb_hierarchies = iotbx.pdb.hierarchy.join_roots(
    roots=pdb_hierarchies)
  output_file_name = "symmetry_copies.pdb"
  print("Writing file:", output_file_name)
  print("""\
REMARK input file name: %s
REMARK original space group: %s
REMARK using two-character chain ids""" % (args[0], str(cs.space_group_info())), file=open(output_file_name, "w"))
  from cctbx import sgtbx
  combined_pdb_hierarchies.write_pdb_or_mmcif_file(
    target_filename=output_file_name,
    crystal_symmetry=cs.customized_copy(
      space_group_info=sgtbx.space_group_info(symbol="P1")),
    append_end=True)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
