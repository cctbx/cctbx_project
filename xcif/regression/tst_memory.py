from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_memory.py
#
# Memory safety and ownership tests for xcif Python bindings.
#
# Verifies that:
#   1. Block objects keep the Document alive (no dangling pointers).
#   2. Loop objects keep the Block/Document alive.
#   3. Repeated parsing does not leak memory.

import gc
import os
import libtbx.load_env
from libtbx.utils import format_cpu_times

CIF_TEXT = (
  "data_mem\n"
  "_cell.length_a 10.5\n"
  "loop_\n"
  "_atom_site.id\n"
  "_atom_site.type_symbol\n"
  "_atom_site.Cartn_x\n"
  "1 C 1.5\n"
  "2 N 3.5\n"
  "3 O 5.5\n"
)

def exercise_block_survives_document():
  """Block must remain valid after the Document Python wrapper is deleted.

  The binding must use return_internal_reference or shared_ptr so that
  the C++ Document (which owns the buffer) stays alive while any Block
  wrapper exists.
  """
  import xcif_ext
  doc = xcif_ext.parse(CIF_TEXT)
  block = doc[0]
  del doc
  gc.collect()
  # Block should still work - buffer must still be alive
  assert block.name == "mem", "Block.name after doc deletion: '%s'" % block.name
  assert block.has_tag("_cell.length_a")
  assert block.find_value("_cell.length_a") == "10.5"

def exercise_loop_survives_block():
  """Loop must remain valid after Block and Document wrappers are deleted."""
  import xcif_ext
  doc = xcif_ext.parse(CIF_TEXT)
  loop = doc[0].find_loop("_atom_site")
  del doc
  gc.collect()
  # Loop should still work
  assert loop.width == 3, "Loop.width after doc deletion: %d" % loop.width
  assert loop.length == 3
  assert loop.value(0, 0) == "1"
  col = loop.column("_atom_site.type_symbol")
  assert col == ["C", "N", "O"]

def exercise_find_block_survives():
  """Block obtained via find_block must also keep Document alive."""
  import xcif_ext
  doc = xcif_ext.parse(CIF_TEXT)
  block = doc.find_block("mem")
  del doc
  gc.collect()
  assert block is not None
  assert block.name == "mem"

def exercise_repeated_parse_no_leak():
  """Parsing the same data many times must not grow memory without bound.

  This is a smoke test, not a precise RSS measurement. It verifies that
  1000 parse calls complete without OOM or crash, which catches gross
  reference-counting bugs.
  """
  import xcif_ext
  for i in range(1000):
    doc = xcif_ext.parse(CIF_TEXT)
    block = doc[0]
    loop = block.find_loop("_atom_site")
    _ = loop.column("_atom_site.Cartn_x")
  gc.collect()
  # If we get here without OOM, the test passes.

def exercise_file_parse_no_leak():
  """File-based parsing must not leak mmap handles or buffers."""
  import xcif_ext
  dist_dir = libtbx.env.dist_path("xcif")
  cif_file = os.path.join(dist_dir, "regression", "example.cif")
  for i in range(200):
    doc = xcif_ext.parse_file(cif_file)
    _ = len(doc)
  gc.collect()

if __name__ == "__main__":
  exercise_block_survives_document()
  exercise_loop_survives_block()
  exercise_find_block_survives()
  exercise_repeated_parse_no_leak()
  exercise_file_parse_no_leak()
  print(format_cpu_times())
  print("OK")
