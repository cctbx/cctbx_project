from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_string_interning.py
#
# String interning tests for xcif Python bindings.
#
# Verifies that the C++/Python boundary interns repeated strings so that
# identical tag names and common values share the same Python str object.
# This reduces memory consumption when accessing large CIF files from
# Python (many atoms share the same element symbols, residue names, etc).

from libtbx.utils import format_cpu_times

def exercise_tag_name_interning():
  """Repeated calls to loop.tags must return the same Python str objects.

  When the binding interns tag names, id(tag_a) == id(tag_b) for the
  same tag string across separate .tags calls.
  """
  import xcif_ext
  cif_text = (
    "data_intern\n"
    "loop_\n"
    "_atom_site.id\n"
    "_atom_site.type_symbol\n"
    "_atom_site.Cartn_x\n"
    "1 C 1.0\n"
    "2 N 2.0\n"
  )
  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_atom_site")

  tags_a = loop.tags
  tags_b = loop.tags
  for i in range(len(tags_a)):
    assert tags_a[i] is tags_b[i], (
      "Tag '%s' not interned: id %d vs %d"
      % (tags_a[i], id(tags_a[i]), id(tags_b[i])))

def exercise_value_interning():
  """Repeated identical values in a column should share Python str objects.

  For example, element symbols like 'C' appearing in many rows should
  all be the same Python object after interning.
  """
  import xcif_ext
  # Build a loop with many repeated element symbols
  lines = ["data_intern_val\n", "loop_\n",
           "_atom_site.id\n", "_atom_site.type_symbol\n"]
  for i in range(100):
    elem = "C" if i % 3 == 0 else ("N" if i % 3 == 1 else "O")
    lines.append("%d %s\n" % (i, elem))
  cif_text = "".join(lines)

  doc = xcif_ext.parse(cif_text)
  loop = doc[0].find_loop("_atom_site")
  col = loop.column("_atom_site.type_symbol")

  # All "C" values should be the same object
  c_values = [v for v in col if v == "C"]
  assert len(c_values) >= 2, "Need at least 2 'C' values for interning test"
  first_c = c_values[0]
  for v in c_values[1:]:
    assert v is first_c, (
      "Value 'C' not interned: id %d vs %d" % (id(first_c), id(v)))

  # Same for "N"
  n_values = [v for v in col if v == "N"]
  assert len(n_values) >= 2
  first_n = n_values[0]
  for v in n_values[1:]:
    assert v is first_n, (
      "Value 'N' not interned: id %d vs %d" % (id(first_n), id(v)))

def exercise_find_value_interning():
  """Repeated find_value calls for the same tag should return the same str."""
  import xcif_ext
  cif_text = (
    "data_fv\n"
    "_cell.length_a 10.5\n"
  )
  doc = xcif_ext.parse(cif_text)
  block = doc[0]

  v1 = block.find_value("_cell.length_a")
  v2 = block.find_value("_cell.length_a")
  assert v1 is v2, (
    "find_value not interned: id %d vs %d" % (id(v1), id(v2)))

def exercise_block_name_interning():
  """Block name from repeated accesses should be the same str object."""
  import xcif_ext
  cif_text = "data_myblock\n_x 1\n"
  doc = xcif_ext.parse(cif_text)

  n1 = doc[0].name
  n2 = doc[0].name
  assert n1 is n2, (
    "Block name not interned: id %d vs %d" % (id(n1), id(n2)))

if __name__ == "__main__":
  exercise_tag_name_interning()
  exercise_value_interning()
  exercise_find_value_interning()
  exercise_block_name_interning()
  print(format_cpu_times())
  print("OK")
