from __future__ import absolute_import, division, print_function
from iotbx.pdb.mmcif import cif_input
import iotbx.cif
from iotbx.cif.builders import crystal_symmetry_builder, CifBuilderError
from cctbx import sgtbx, uctbx

def tst_bad_symmetry_1():
  """ NMR structure with question mark in unit cell parameters.
  Expected outcome - as if absent completely.
  """
  cif_txt_1 = """\
data_1JZC
_cell.entry_id           1JZC
_cell.length_a           ?
_cell.length_b           ?
_cell.length_c           ?
_cell.angle_alpha        ?
_cell.angle_beta         ?
_cell.angle_gamma        ?
_cell.Z_PDB              1
_cell.pdbx_unique_axis   ?
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1   O "O5'"  . G A 1 1  ? -4.337  -11.083 -6.294 1.00 2.13 ? 1  G A "O5'"  1
ATOM 2   C "C5'"  . G A 1 1  ? -5.629  -11.603 -6.617 1.00 2.29 ? 1  G A "C5'"  1
ATOM 3   C "C4'"  . G A 1 1  ? -6.489  -11.784 -5.370 1.00 2.23 ? 1  G A "C4'"  1
"""
  c_inp = cif_input(source_info=None, lines=cif_txt_1)
  cs = c_inp.crystal_symmetry()
  assert cs.unit_cell() == None
  assert cs.space_group() == None
  assert cs.is_empty()

def tst_bad_symmetry_2():
  """ Same as #1, but no parameters at all.
  Expected outcome - as if absent completely.
  """
  cif_txt_1 = """\
data_1JZC
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1   O "O5'"  . G A 1 1  ? -4.337  -11.083 -6.294 1.00 2.13 ? 1  G A "O5'"  1
ATOM 2   C "C5'"  . G A 1 1  ? -5.629  -11.603 -6.617 1.00 2.29 ? 1  G A "C5'"  1
ATOM 3   C "C4'"  . G A 1 1  ? -6.489  -11.784 -5.370 1.00 2.23 ? 1  G A "C4'"  1
"""
  c_inp = cif_input(source_info=None, lines=cif_txt_1)
  cs = c_inp.crystal_symmetry()
  assert cs.unit_cell() == None
  assert cs.space_group() == None
  assert cs.is_empty()


def build(cif_txt, strict=False):
  """Crystal symmetry from the first data block of a CIF string."""
  block = list(iotbx.cif.reader(input_string=cif_txt).model().values())[0]
  return crystal_symmetry_builder(block, strict=strict).crystal_symmetry

monoclinic_cell = """\
data_x
_cell.length_a 52.7
_cell.length_b 26.18
_cell.length_c 51.71
_cell.angle_alpha 90
_cell.angle_beta 93
_cell.angle_gamma 90
"""

c_unique_cell = monoclinic_cell.replace(
  "_cell.angle_beta 93", "_cell.angle_beta 90").replace(
  "_cell.angle_gamma 90", "_cell.angle_gamma 93")

def tst_pdb_alternative_settings():
  """ H-M symbols the wwPDB uses for alternative settings, mapped through
  the table of the CRYST1 reader; the Hall symbols are the table entries.
  """
  template = """\
data_x
_cell.length_a %s
_cell.length_b %s
_cell.length_c %s
_cell.angle_alpha %s
_cell.angle_beta %s
_cell.angle_gamma %s
_symmetry.space_group_name_H-M '%s'
_symmetry.Int_Tables_number %d
"""
  cases = [
    ("C 1 21 1", 5, (52.7, 26.18, 51.71, 90, 93, 90),
     "C 2y (x+1/4,y+1/4,z)"),
    ("I 1 21 1", 5, (124.49, 42.28, 92.27, 90, 92.62, 90),
     "C 2y (x+1/4,y+1/4,-x+z-1/4)"),
    ("F 4 2 2", 97, (146.8, 146.8, 155.4, 90, 90, 90),
     "I 4 2 (1/2*x+1/2*y,-1/2*x+1/2*y,z)"),
    ("C 4 21 2", 90, (134.5, 134.5, 85.9, 90, 90, 90),
     "P 4 2 (1/2*x+1/2*y-1/4,-1/2*x+1/2*y-1/4,z)"),
    ("P 21 21 2 A", 18, (104.013, 74.078, 69.61, 90, 90, 90),
     "P 2 2ab (x+1/4,y+1/4,z)"),
  ]
  for symbol, number, cell, hall in cases:
    cs = build(template % (cell + (symbol, number)))
    expected = sgtbx.space_group_info("Hall: " + hall)
    assert cs.space_group() == expected.group(), (symbol, cs.space_group_info())
    assert cs.unit_cell().is_similar_to(uctbx.unit_cell(cell))

def tst_unrecognized_symbol_is_ignored():
  """ An H-M symbol neither sgtbx nor the table knows is skipped, as
  before the table lookup existed.
  """
  cs = build(
    monoclinic_cell + "_symmetry.space_group_name_H-M 'P 1 2 1 (No. 3)'\n")
  assert cs.space_group() is None

def tst_alternative_setting_symbol_yields_to_operators():
  """ C 1 21 1 next to a loop of standard C 1 2 1 operators: the operators
  define the group, as before the table lookup existed.
  """
  cif_txt = monoclinic_cell + """\
_symmetry.space_group_name_H-M 'C 1 21 1'
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 x,y,z
2 -x,y,-z
3 x+1/2,y+1/2,z
4 -x+1/2,y+1/2,-z
"""
  cs = build(cif_txt)
  assert cs.space_group() == sgtbx.space_group_info("C 1 2 1").group()

def tst_alternative_setting_symbol_yields_to_hall():
  """ C 1 21 1 next to a Hall symbol that parses: the Hall symbol defines
  the group, as before the table lookup existed.
  """
  cif_txt = monoclinic_cell + """\
_symmetry.space_group_name_Hall ' C 2y'
_symmetry.space_group_name_H-M 'C 1 21 1'
"""
  cs = build(cif_txt)
  assert cs.space_group() == sgtbx.space_group_info("C 1 2 1").group()

def tst_alternative_setting_symbol_not_fitting_cell():
  """ The table gives a b-unique group for C 1 21 1; with a c-unique cell
  there is no space group, as before the table lookup existed, rather
  than an incompatibility error.
  """
  cs = build(c_unique_cell + "_symmetry.space_group_name_H-M 'C 1 21 1'\n")
  assert cs.space_group() is None
  assert cs.unit_cell() is not None


def tst_bad_symop_is_cif_builder_error():
  """ A symmetry operator sgtbx cannot parse is reported as a
  CifBuilderError quoting the operator as written, not as the raw sgtbx
  exception. Only one trailing comma is forgiven.
  """
  cif_txt = """\
data_x
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 x,y,z,,
"""
  try:
    build(cif_txt)
  except CifBuilderError as e:
    assert "Error interpreting symmetry operator" in str(e), str(e)
    assert "x,y,z,," in str(e), str(e)
  else:
    raise AssertionError("CifBuilderError expected")

def tst_symop_trailing_comma():
  """ PDB entry 1l6j writes every operator with one trailing comma
  ("+X,+Y,+Z,"). Tolerate it.
  """
  cif_txt = """\
data_x
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 +X,+Y,+Z,
2 -X,-Y,+Z,
"""
  cs = build(cif_txt)
  assert cs.space_group() == sgtbx.space_group_info("P 1 1 2").group()


def tst_space_group_from_it_number_only():
  """ _symmetry.Int_Tables_number alone (no symbol, no operators) gives
  the standard setting of that number.
  """
  cif_txt = """\
data_x
_cell.length_a 10
_cell.length_b 20
_cell.length_c 30
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.Int_Tables_number 19
"""
  cs = build(cif_txt)
  assert cs.space_group() == sgtbx.space_group_info(number=19).group()

def tst_it_number_yields_to_operators():
  """ The number only knows the standard setting. Operators in another
  setting (I 1 2 1 for number 5) define the group; the number is ignored
  rather than reported as inconsistent.
  """
  cif_txt = monoclinic_cell + """\
_symmetry.Int_Tables_number 5
loop_
_symmetry_equiv.id
_symmetry_equiv.pos_as_xyz
1 x,y,z
2 -x,y,-z
3 x+1/2,y+1/2,z+1/2
4 -x+1/2,y+1/2,-z+1/2
"""
  cs = build(cif_txt)
  assert cs.space_group() == sgtbx.space_group_info("I 1 2 1").group()

def tst_it_number_not_fitting_cell():
  """ Number 4 is P 1 21 1, but the cell is c-unique: no space group,
  as when the number was not read at all. In strict mode the misfit is
  reported, as it was for a core CIF number.
  """
  cif_txt = c_unique_cell + "_symmetry.Int_Tables_number 4\n"
  cs = build(cif_txt)
  assert cs.space_group() is None
  assert cs.unit_cell() is not None
  try:
    build(cif_txt, strict=True)
  except CifBuilderError as e:
    assert "incompatible with unit cell" in str(e), str(e)
  else:
    raise AssertionError("CifBuilderError expected")

def tst_it_number_with_partial_cell():
  """ A number cannot complete a partial cell: the usual error, not a
  raw exception from the cell inference.
  """
  cif_txt = """\
data_x
_cell.length_a 10
_cell.length_b 20
_cell.length_c 30
_symmetry.Int_Tables_number 4
"""
  try:
    build(cif_txt)
  except CifBuilderError as e:
    assert "Not all unit cell parameters" in str(e), str(e)
  else:
    raise AssertionError("CifBuilderError expected")


if __name__ == '__main__':
  tst_bad_symmetry_1()
  tst_bad_symmetry_2()
  tst_pdb_alternative_settings()
  tst_unrecognized_symbol_is_ignored()
  tst_alternative_setting_symbol_yields_to_operators()
  tst_alternative_setting_symbol_yields_to_hall()
  tst_alternative_setting_symbol_not_fitting_cell()
  tst_bad_symop_is_cif_builder_error()
  tst_symop_trailing_comma()
  tst_space_group_from_it_number_only()
  tst_it_number_yields_to_operators()
  tst_it_number_not_fitting_cell()
  tst_it_number_with_partial_cell()
