from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_example_parse.py
#
# Regression guard for the example_parse demo.
#
# Runs the compiled example_parse binary against example.cif and checks:
#   1. The process exits cleanly.
#   2. Structural facts (block count, loop dimensions, tag values).
#   3. Numeric conversion results (as_double, as_int).
#   4. Null value detection (unknown ?, inapplicable .).

import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.utils import format_cpu_times

def exercise():
  build_dir = libtbx.env.under_build("xcif")
  dist_dir  = libtbx.env.dist_path("xcif")

  binary = os.path.join(build_dir, "regression", "cpp", "example_parse")
  if not os.path.isfile(binary) and os.path.isfile(binary + ".exe"):
    binary += ".exe"
  assert os.path.isfile(binary), \
    "example_parse binary not found: %s\n(run libtbx.scons first)" % binary

  cif_file = os.path.join(dist_dir, "regression", "example.cif")
  assert os.path.isfile(cif_file), \
    "example.cif not found: %s" % cif_file

  result = easy_run.fully_buffered(
    command='"%s" "%s"' % (binary, cif_file))

  assert result.return_code == 0, (
    "example_parse exited with code %d\nstderr:\n%s" % (
      result.return_code, "\n".join(result.stderr_lines)))

  out = "\n".join(result.stdout_lines)

  # -----------------------------------------------------------------------
  # Structure
  # -----------------------------------------------------------------------
  assert "Blocks: 1" in out, "Expected 1 block"
  assert "Block: 1UBQ" in out, "Expected block name 1UBQ"
  assert "Loops: 5" in out, "Expected 5 loops"
  assert "13 tags x 10 rows" in out, "Expected atom loop dimensions"
  assert "6 tags x 1 rows" in out, "Expected reflns loop dimensions"

  # -----------------------------------------------------------------------
  # Tag-value lookups
  # -----------------------------------------------------------------------
  assert "_cell.length_a   = 50.840(10)" in out, \
    "Expected cell length_a value with SU"
  assert "P 21 21 21" in out, \
    "Expected space group name (quoted string)"

  # -----------------------------------------------------------------------
  # Numeric conversions
  # -----------------------------------------------------------------------
  assert "as_double: 50.840" in out, \
    "Expected as_double for cell length"
  assert "as_double: 0.2410" in out, \
    "Expected as_double for R_free"
  assert "as_int: 19" in out, \
    "Expected as_int for space group number"

  # -----------------------------------------------------------------------
  # Standard uncertainty (SU) extraction
  # -----------------------------------------------------------------------
  assert "as_double_with_su: value=50.840  su=0.010" in out, \
    "Expected as_double_with_su for cell length_a"
  assert "value=42.770  su=0.001" in out, \
    "Expected SU for cell length_b"
  assert "value=28.950  su=2.000" in out, \
    "Expected SU for cell length_c"

  # -----------------------------------------------------------------------
  # Null value handling
  # -----------------------------------------------------------------------
  assert "[inapplicable]" in out, \
    "Expected ? detected as inapplicable"
  assert "[null" in out and "NaN: yes" in out, \
    "Expected . detected as null → NaN"
  assert "has_tag(_cell.volume)  = true" in out, \
    "Expected has_tag true"
  assert "has_tag(_nonexistent)  = false" in out, \
    "Expected has_tag false"

  # -----------------------------------------------------------------------
  # Atom coordinate table
  # -----------------------------------------------------------------------
  assert "27.340" in out and "24.430" in out and "2.614" in out, \
    "Expected first atom coordinates"
  assert "column_index" in out, \
    "Expected column_index in output header"
  assert "5 more atoms" in out, \
    "Expected truncation message for atoms"

  # -----------------------------------------------------------------------
  # Column extraction
  # -----------------------------------------------------------------------
  assert "as_double = 1.80" in out, \
    "Expected column extraction result"

if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
