from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_pdb_sf.py
#
# Step 13: Real-world validation — PDB structure factor file (1yjp-sf.cif).
# Parses a real mmCIF structure factor file with both xcif and iotbx.cif.reader,
# then verifies that all extracted data matches.

import os
import math
import xcif
import xcif_ext
import iotbx.cif
import libtbx.load_env
from scitbx.array_family import flex
from libtbx.utils import format_cpu_times
from libtbx.test_utils import approx_equal


def _get_cif_path():
  dist_dir = libtbx.env.dist_path("xcif")
  return os.path.join(dist_dir, "regression", "1yjp-sf.cif")


def _parse_both(cif_path):
  """Parse with both engines, return (xcif_doc, iotbx_model)."""
  doc = xcif_ext.parse_file(cif_path)
  cif_model = iotbx.cif.reader(file_path=cif_path).model()
  return doc, cif_model


def exercise_block_discovery():
  """Verify block name and case-insensitive lookup."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  assert len(doc) == 1, "Expected 1 block, got %d" % len(doc)
  assert len(cif_model) == 1

  assert doc[0].name == "r1yjpsf", "Block name: '%s'" % doc[0].name
  assert list(cif_model.keys())[0] == "r1yjpsf"

  # Case-insensitive
  assert doc.find_block("R1YJPSF") is not None
  assert cif_model.get("R1YJPSF") is not None


def exercise_cell_symmetry():
  """Cell and symmetry parameters match coordinate file values."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xb = doc[0]
  ib = list(cif_model.values())[0]

  # Cell parameters
  cell_tags = [
    ("_cell.length_a", "21.937"),
    ("_cell.length_b", "4.866"),
    ("_cell.length_c", "23.477"),
    ("_cell.angle_alpha", "90.000"),
    ("_cell.angle_beta", "107.080"),
    ("_cell.angle_gamma", "90.000"),
  ]
  for tag, expected in cell_tags:
    xval = xb.find_value(tag)
    ival = ib[tag]
    assert xval == ival, \
      "%s: xcif='%s' vs iotbx='%s'" % (tag, xval, ival)
    assert xval == expected, \
      "%s: got '%s', expected '%s'" % (tag, xval, expected)

  # Space group
  xsg = xb.find_value("_symmetry.space_group_name_H-M")
  isg = ib["_symmetry.space_group_name_H-M"]
  assert xsg == isg, "space_group: xcif='%s' vs iotbx='%s'" % (xsg, isg)
  assert xsg == "P 1 21 1", "space_group: '%s'" % xsg


def exercise_refln_loop():
  """_refln loop dimensions and all string columns match."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_refln")
  iblock = list(cif_model.values())[0]
  iloop = iblock.get_loop("_refln")

  assert xloop is not None, "xcif: _refln loop not found"
  assert iloop is not None, "iotbx: _refln loop not found"

  # Dimensions
  assert xloop.width == 9, "width: %d" % xloop.width
  assert xloop.length == 495, "length: %d" % xloop.length
  assert xloop.width == iloop.n_columns(), \
    "width: xcif=%d vs iotbx=%d" % (xloop.width, iloop.n_columns())
  assert xloop.length == iloop.n_rows(), \
    "length: xcif=%d vs iotbx=%d" % (xloop.length, iloop.n_rows())

  # Tag names match
  xtags = list(xloop.tags)
  itags = list(iloop.keys())
  assert xtags == itags, \
    "tags differ: xcif=%s vs iotbx=%s" % (xtags, itags)

  # All columns match cell-by-cell
  for tag in xtags:
    xcol = xloop.column(tag)
    icol = list(iloop[tag])
    assert xcol == icol, \
      "Column %s row 0: xcif='%s' vs iotbx='%s'" % (tag, xcol[0], icol[0])


def exercise_miller_indices():
  """Miller indices extracted as flex.int."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_refln")
  iblock = list(cif_model.values())[0]

  for tag, expected_first in [("_refln.index_h", -12),
                               ("_refln.index_k", 0),
                               ("_refln.index_l", 2)]:
    xarr = xloop.column_as_flex_int(tag)
    assert xarr.size() == 495, "%s size: %d" % (tag, xarr.size())
    assert xarr[0] == expected_first, \
      "%s[0]: got %d, expected %d" % (tag, xarr[0], expected_first)
    # Compare against iotbx
    iarr_str = iblock[tag]
    for j in range(xarr.size()):
      assert xarr[j] == int(iarr_str[j]), \
        "%s[%d]: xcif=%d vs iotbx=%d" % (tag, j, xarr[j], int(iarr_str[j]))


def exercise_amplitudes():
  """Structure factor amplitudes and sigmas as flex.double."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_refln")
  iblock = list(cif_model.values())[0]

  for tag, expected_first in [("_refln.F_meas_au", 3.49),
                               ("_refln.F_meas_sigma_au", 2.34)]:
    xarr = xloop.column_as_flex_double(tag)
    assert xarr.size() == 495, "%s size: %d" % (tag, xarr.size())
    assert approx_equal(xarr[0], expected_first), \
      "%s[0]: got %f, expected %f" % (tag, xarr[0], expected_first)
    # No unexpected NaN
    for j in range(xarr.size()):
      assert not math.isnan(xarr[j]), "%s[%d] is NaN" % (tag, j)
    # Compare against iotbx
    iarr_str = iblock[tag]
    for j in range(xarr.size()):
      assert approx_equal(xarr[j], float(iarr_str[j])), \
        "%s[%d]: xcif=%f vs iotbx=%f" % (tag, j, xarr[j], float(iarr_str[j]))


def exercise_status_column():
  """Status column as flex.std_string — mix of 'o' and 'f'."""
  cif_path = _get_cif_path()
  doc, cif_model = _parse_both(cif_path)

  xloop = doc[0].find_loop("_refln")
  iblock = list(cif_model.values())[0]

  xstatus = xloop.column_as_flex_string("_refln.status")
  istatus = iblock["_refln.status"]

  assert xstatus.size() == 495
  assert list(xstatus) == list(istatus), "status columns differ"

  # Verify mix of 'o' and 'f'
  n_obs = list(xstatus).count("o")
  n_free = list(xstatus).count("f")
  assert n_obs > 0, "No observed reflections"
  assert n_free > 0, "No free-set reflections"
  assert n_obs + n_free == 495, \
    "obs(%d) + free(%d) != 495" % (n_obs, n_free)


def exercise_model_level():
  """Model-level API: xcif.reader vs iotbx.cif.reader."""
  cif_path = _get_cif_path()
  m_xcif = xcif.reader(file_path=cif_path).model()
  m_ucif = iotbx.cif.reader(file_path=cif_path).model()

  assert list(m_xcif.keys()) == list(m_ucif.keys())

  bx = m_xcif["r1yjpsf"]
  bu = m_ucif["r1yjpsf"]

  # Scalar tags
  for tag in ["_cell.length_a", "_cell.length_b", "_cell.length_c",
              "_symmetry.space_group_name_H-M",
              "_entry.id"]:
    assert bx[tag] == bu[tag], \
      "Tag %s: xcif='%s' vs ucif='%s'" % (tag, bx[tag], bu[tag])

  # Reflection columns via model dict access
  for tag in ["_refln.index_h", "_refln.index_k", "_refln.index_l",
              "_refln.F_meas_au", "_refln.F_meas_sigma_au",
              "_refln.status"]:
    col_x = list(bx[tag])
    col_u = list(bu[tag])
    assert col_x == col_u, \
      "Column %s: xcif=%s vs ucif=%s" % (tag, col_x[:3], col_u[:3])


if __name__ == "__main__":
  exercise_block_discovery()
  exercise_cell_symmetry()
  exercise_refln_loop()
  exercise_miller_indices()
  exercise_amplitudes()
  exercise_status_column()
  exercise_model_level()
  print(format_cpu_times())
  print("OK")
