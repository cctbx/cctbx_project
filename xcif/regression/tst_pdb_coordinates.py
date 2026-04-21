from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_pdb_coordinates.py
#
# Step 13: Real-world validation — PDB coordinate file (1yjp.cif).
# Parses a real mmCIF coordinate file with both xcif and iotbx.cif.reader,
# then verifies that all extracted data matches.

import os
import xcif
import xcif_ext
import iotbx.cif
import libtbx.load_env
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal


def _get_cif_path():
  dist_dir = libtbx.env.dist_path("xcif")
  return os.path.join(dist_dir, "regression", "1yjp.cif")


def _parse_both(cif_path):
  """Parse with both engines, return (xcif_doc, iotbx_model)."""
  doc = xcif_ext.parse_file(cif_path)
  cif_model = iotbx.cif.reader(file_path=cif_path).model()
  return doc, cif_model


def exercise_block_discovery():
  """Verify block name and case-insensitive lookup."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  # Exactly one block
  assert len(doc) == 1, "Expected 1 block, got %d" % len(doc)
  assert len(cif_model) == 1, "Expected 1 block in iotbx, got %d" % len(cif_model)

  # Block name preserved
  assert doc[0].name == "1YJP", "Block name: '%s'" % doc[0].name
  assert list(cif_model.keys())[0] == "1YJP"

  # Case-insensitive lookup
  assert doc.find_block("1yjp") is not None
  assert doc.find_block("1YJP") is not None
  assert cif_model.get("1yjp") is not None


def exercise_cell_parameters():
  """Unit cell parameters match exactly between parsers."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xb = doc[0]
  ib = list(cif_model.values())[0]

  cell_tags = [
    ("_cell.length_a", "21.937"),
    ("_cell.length_b", "4.866"),
    ("_cell.length_c", "23.477"),
    ("_cell.angle_alpha", "90.00"),
    ("_cell.angle_beta", "107.08"),
    ("_cell.angle_gamma", "90.00"),
  ]
  for tag, expected in cell_tags:
    xval = xb.find_value(tag)
    ival = ib[tag]
    assert xval == ival, \
      "%s: xcif='%s' vs iotbx='%s'" % (tag, xval, ival)
    assert xval == expected, \
      "%s: got '%s', expected '%s'" % (tag, xval, expected)


def exercise_symmetry():
  """Space group and symmetry number match."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xb = doc[0]
  ib = list(cif_model.values())[0]

  xsg = xb.find_value("_symmetry.space_group_name_H-M")
  isg = ib["_symmetry.space_group_name_H-M"]
  assert xsg == isg, "space_group: xcif='%s' vs iotbx='%s'" % (xsg, isg)
  assert xsg == "P 1 21 1", "space_group: '%s'" % xsg

  xnum = xb.find_value("_symmetry.Int_Tables_number")
  inum = ib["_symmetry.Int_Tables_number"]
  assert xnum == inum, \
    "Int_Tables_number: xcif='%s' vs iotbx='%s'" % (xnum, inum)
  assert xnum == "4", "Int_Tables_number: '%s'" % xnum


def exercise_atom_site_loop():
  """_atom_site loop dimensions and all string columns match."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_atom_site")
  iblock = list(cif_model.values())[0]
  iloop = iblock.get_loop("_atom_site")

  assert xloop is not None, "xcif: _atom_site loop not found"
  assert iloop is not None, "iotbx: _atom_site loop not found"

  # Dimensions
  assert xloop.width == 21, "width: %d" % xloop.width
  assert xloop.length == 66, "length: %d" % xloop.length
  assert xloop.width == iloop.n_columns(), \
    "width: xcif=%d vs iotbx=%d" % (xloop.width, iloop.n_columns())
  assert xloop.length == iloop.n_rows(), \
    "length: xcif=%d vs iotbx=%d" % (xloop.length, iloop.n_rows())

  # Tag names match
  xtags = list(xloop.tags)
  itags = list(iloop.keys())
  assert xtags == itags, \
    "tags differ: xcif=%s vs iotbx=%s" % (xtags[:5], itags[:5])

  # All string columns match cell-by-cell
  for tag in xtags:
    xcol = xloop.column(tag)
    icol = list(iloop[tag])
    assert xcol == icol, \
      "Column %s row 0: xcif='%s' vs iotbx='%s'" % (tag, xcol[0], icol[0])


def exercise_atom_site_numeric():
  """Numeric extraction from _atom_site: coordinates, occupancy, B-factor, id."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_atom_site")
  iblock = list(cif_model.values())[0]

  # Cartesian coordinates as flex.double
  for tag, expected_first in [("_atom_site.Cartn_x", -9.009),
                               ("_atom_site.Cartn_y", 4.612),
                               ("_atom_site.Cartn_z", 6.102)]:
    xarr = xloop.column_as_flex_double(tag)
    assert xarr.size() == 66, "%s size: %d" % (tag, xarr.size())
    assert approx_equal(xarr[0], expected_first), \
      "%s[0]: got %f, expected %f" % (tag, xarr[0], expected_first)
    # Compare against iotbx string values
    iarr_str = iblock[tag]
    for j in range(xarr.size()):
      assert approx_equal(xarr[j], float(iarr_str[j])), \
        "%s[%d]: xcif=%f vs iotbx=%f" % (tag, j, xarr[j], float(iarr_str[j]))

  # Occupancy and B-factor
  occ = xloop.column_as_flex_double("_atom_site.occupancy")
  assert occ.size() == 66
  assert approx_equal(occ[0], 1.00)

  biso = xloop.column_as_flex_double("_atom_site.B_iso_or_equiv")
  assert biso.size() == 66
  assert approx_equal(biso[0], 16.77)

  # Atom ID as flex.int
  ids = xloop.column_as_flex_int("_atom_site.id")
  assert ids.size() == 66
  assert ids[0] == 1
  assert ids[65] == 66

  # Type symbol spot-check
  syms = xloop.column_as_flex_string("_atom_site.type_symbol")
  assert syms.size() == 66
  assert syms[0] == "N", "First type_symbol: '%s'" % syms[0]

  # group_PDB spot-check: first 59 are ATOM, last 7 are HETATM
  groups = xloop.column_as_flex_string("_atom_site.group_PDB")
  assert groups[0] == "ATOM"
  assert groups[58] == "ATOM"
  assert groups[59] == "HETATM"
  assert groups[65] == "HETATM"


def exercise_other_loops():
  """Verify loops and get_loop_or_row match between parsers."""
  cif_path = _get_cif_path()
  m_xcif = xcif.reader(file_path=cif_path).model()
  m_ucif = iotbx.cif.reader(file_path=cif_path).model()

  bx = m_xcif["1YJP"]
  bu = m_ucif["1YJP"]

  # Real loops: get_loop() returns matching dimensions and data
  real_loop_tags = [
    ("_database_2", "_database_2.database_id", 3),
    ("_audit_author", "_audit_author.name", 7),
    ("_entity", "_entity.id", 2),
    ("_chem_comp", "_chem_comp.id", 5),
  ]
  for cat, tag, expected_rows in real_loop_tags:
    lx = bx.get_loop(cat)
    lu = bu.get_loop(cat)
    assert lx is not None, "xcif: get_loop(%s) is None" % cat
    assert lu is not None, "iotbx: get_loop(%s) is None" % cat
    assert lx.n_rows() == lu.n_rows(), \
      "%s n_rows: xcif=%d vs iotbx=%d" % (cat, lx.n_rows(), lu.n_rows())
    assert lx.n_columns() == lu.n_columns(), \
      "%s n_cols: xcif=%d vs iotbx=%d" % (cat, lx.n_columns(), lu.n_columns())
    col_x = list(lx[tag])
    col_u = list(lu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs iotbx=%s" % (tag, col_x[:3], col_u[:3])

  # Scalar categories: get_loop_or_row() synthesizes single-row loops
  scalar_cats = ["_cell", "_symmetry", "_entry"]
  for cat in scalar_cats:
    lx = bx.get_loop_or_row(cat)
    lu = bu.get_loop_or_row(cat)
    assert lx is not None, "xcif: get_loop_or_row(%s) is None" % cat
    assert lu is not None, "iotbx: get_loop_or_row(%s) is None" % cat
    assert lx.n_rows() == lu.n_rows(), \
      "%s n_rows: xcif=%d vs iotbx=%d" % (cat, lx.n_rows(), lu.n_rows())
    assert lx.n_columns() == lu.n_columns(), \
      "%s n_cols: xcif=%d vs iotbx=%d" % (cat, lx.n_columns(), lu.n_columns())
    # Compare all column values
    for tag in lx.keys():
      val_x = list(lx[tag])
      val_u = list(lu[tag])
      assert val_x == val_u, \
        "%s: xcif=%s vs iotbx=%s" % (tag, val_x, val_u)


def exercise_model_level():
  """Model-level API: xcif.reader vs iotbx.cif.reader."""
  cif_path = _get_cif_path()
  m_xcif = xcif.reader(file_path=cif_path).model()
  m_ucif = iotbx.cif.reader(file_path=cif_path).model()

  # Block names
  assert list(m_xcif.keys()) == list(m_ucif.keys())

  bx = m_xcif["1YJP"]
  bu = m_ucif["1YJP"]

  # Scalar tags
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_cell.angle_alpha", "_cell.angle_beta", "_cell.angle_gamma",
              "_symmetry.space_group_name_H-M",
              "_symmetry.Int_Tables_number",
              "_entry.id"]:
    assert bx[tag] == bu[tag], \
      "Tag %s: xcif='%s' vs ucif='%s'" % (tag, bx[tag], bu[tag])

  # Loop columns via model dict access -> flex.std_string
  for tag in ["_atom_site.Cartn_x", "_atom_site.Cartn_y",
              "_atom_site.Cartn_z", "_atom_site.type_symbol",
              "_atom_site.label_comp_id", "_atom_site.id"]:
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs ucif=%s" % (tag, col_x[:3], col_u[:3])

  # get_loop_or_row for scalar categories returns matching data
  for cat in ["_cell", "_symmetry"]:
    lx = bx.get_loop_or_row(cat)
    lu = bu.get_loop_or_row(cat)
    assert lx is not None, "xcif: get_loop_or_row(%s) None" % cat
    assert lu is not None, "iotbx: get_loop_or_row(%s) None" % cat
    assert sorted(lx.keys()) == sorted(lu.keys()), \
      "%s keys differ: xcif=%s vs ucif=%s" % (cat, sorted(lx.keys()), sorted(lu.keys()))


if __name__ == "__main__":
  exercise_block_discovery()
  exercise_cell_parameters()
  exercise_symmetry()
  exercise_atom_site_loop()
  exercise_atom_site_numeric()
  exercise_other_loops()
  exercise_model_level()
  print(format_cpu_times())
  print("OK")
