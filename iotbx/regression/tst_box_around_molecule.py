from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_file

pdb_str = """\
ATOM      1  N   SER 1  71      10.595  29.158 114.534  1.00300.00           N
ATOM      2  CA  SER 1  71      10.156  30.582 114.477  1.00300.00           C
ATOM      3  C   SER 1  71      11.003  31.474 115.392  1.00300.00           C
ATOM      4  O   SER 1  71      10.892  31.400 116.622  1.00300.00           O
"""

# ------------------------------------------------------------------------------

cif_str = """\
data_default
loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1 N . PRO H 14 ? 5.39400 6.71500 9.04100 1.000 20.43000 N ? A ? 1 1
   ATOM 2 CA . PRO H 14 ? 6.20600 7.23200 10.14700 1.000 20.01000 C ? A ? 1 1
   ATOM 3 C . PRO H 14 ? 7.24800 8.27800 9.73800 1.000 20.76000 C ? A ? 1 1
   ATOM 4 O . PRO H 14 ? 7.79100 8.23500 8.63500 1.000 16.34000 O ? A ? 1 1
"""
# ------------------------------------------------------------------------------

def exercise(prefix="tst_pdb_box_around_mol"):
  # test pdb format
  with open("%s.pdb" % prefix, 'w') as f:
    f.write(pdb_str)
  assert not easy_run.call("iotbx.pdb.box_around_molecule %s.pdb" % (
      prefix))
  assert_lines_in_file(file_name="%s_box_000.pdb" % prefix,
    lines="""\
CRYST1   10.847   12.316   12.145  90.00  90.00  90.00 P 1
SCALE1      0.092191  0.000000  0.000000        0.00000
SCALE2      0.000000  0.081195  0.000000        0.00000
SCALE3      0.000000  0.000000  0.082338        0.00000
ATOM      1  N   SER 1  71       5.439   5.000   5.057  1.00300.00           N
ATOM      2  CA  SER 1  71       5.000   6.424   5.000  1.00300.00           C
ATOM      3  C   SER 1  71       5.847   7.316   5.915  1.00300.00           C
ATOM      4  O   SER 1  71       5.736   7.242   7.145  1.00300.00           O
TER""")
  # test mmcif format
  if libtbx.env.find_in_repositories(relative_path='chem_data') is not None:
    with open("%s.cif" % prefix, 'w') as f:
      f.write(cif_str)
    assert not easy_run.call("iotbx.pdb.box_around_molecule %s.cif" % (
        prefix))
    # XXX temp fix for python3 failure
    # XXX changes in python 3 because _cell.length_a and other order differs

    # Python 3 behavior looks better.
    # Reason for this discrepancy is that starting from Python 3.7
    # "the insertion-order preservation nature of dict objects has been
    # declared to be an official part of the Python language spec."
    # https://docs.python.org/3/whatsnew/3.7.html
    # Basically, dict.keys() in Py3 preserves the order in which keys were added
    # to the dict.
    #
    # Sorting of mmCIF keys here
    # cctbx_project/iotbx/cif/__init__.py, def category_sort_function()
    # should be improved to address such issues universally (sorting within
    # category), note the diference in _space_group and _symmetry categories
    # as well.

    import sys
    if sys.version_info.major == 3:
      assert_lines_in_file(file_name="%s_box_000.cif" % prefix,
      lines="""\
data_default
_cell.length_a                    12.397
_cell.length_b                    11.563
_cell.length_c                    11.512
_cell.angle_alpha                 90.000
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.volume                      1650.205
_space_group.crystal_system       triclinic
_space_group.IT_number            1
_space_group.name_H-M_alt         'P 1'
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.space_group_name_Hall   ' P 1'
_symmetry.Int_Tables_number       1
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z""")
    else: # version 2.7
      assert_lines_in_file(file_name="%s_box_000.cif" % prefix,
      lines="""\
data_default
_cell.angle_beta                  90.000
_cell.angle_gamma                 90.000
_cell.length_b                    11.563
_cell.length_c                    11.512
_cell.angle_alpha                 90.000
_cell.volume                      1650.205
_cell.length_a                    12.397
_space_group.crystal_system       triclinic
_space_group.name_H-M_alt         'P 1'
_space_group.IT_number            1
_space_group.name_Hall            ' P 1'
_symmetry.space_group_name_H-M    'P 1'
_symmetry.Int_Tables_number       1
_symmetry.space_group_name_Hall   ' P 1'
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z""")

  else:
    print('chem_data is not available, skipping')

  print("OK")

# ------------------------------------------------------------------------------

if __name__=="__main__":
  exercise()
