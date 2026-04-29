from __future__ import absolute_import, division, print_function
# cctbx_project/xcif/regression/tst_example_tokenize.py
#
# Regression guard for the example_tokenize demo and example.cif fixture.
#
# Runs the compiled example_tokenize binary against example.cif and checks:
#   1. The process exits cleanly.
#   2. The total token count and per-type breakdown are unchanged.
#   3. Specific token texts are present (block header, quoted string,
#      semicolon field, scientific notation, unknown/missing values).
#
# If the tokenizer logic changes behaviour on example.cif, or if example.cif
# is edited in a way that alters the token stream, this test will catch it.

import os
import libtbx.load_env
from libtbx import easy_run
from libtbx.utils import format_cpu_times

# ---------------------------------------------------------------------------
# Expected golden values — update here when example.cif is intentionally
# changed, and document the reason in the commit message.
# ---------------------------------------------------------------------------
EXPECTED_TOTAL          = 216
EXPECTED_TAG            = 45
EXPECTED_VALUE          = 165
EXPECTED_BLOCK_HEADER   = 1
EXPECTED_LOOP           = 5

def exercise():
  build_dir = libtbx.env.under_build("xcif")
  dist_dir  = libtbx.env.dist_path("xcif")

  # Locate the compiled binary (handle Windows .exe suffix).
  binary = os.path.join(build_dir, "regression", "cpp", "example_tokenize")
  if not os.path.isfile(binary) and os.path.isfile(binary + ".exe"):
    binary += ".exe"
  assert os.path.isfile(binary), \
    "example_tokenize binary not found: %s\n(run libtbx.scons first)" % binary

  # Locate the fixture CIF in the source tree.
  cif_file = os.path.join(dist_dir, "regression", "example.cif")
  assert os.path.isfile(cif_file), \
    "example.cif not found: %s" % cif_file

  # Run the demo.
  result = easy_run.fully_buffered(
    command='"%s" "%s"' % (binary, cif_file))

  assert result.return_code == 0, (
    "example_tokenize exited with code %d\nstderr:\n%s" % (
      result.return_code, "\n".join(result.stderr_lines)))

  out = "\n".join(result.stdout_lines)

  # -----------------------------------------------------------------------
  # Golden token counts — the primary regression guard.
  # Any change to the tokenizer or to example.cif that shifts these numbers
  # is caught here.
  # -----------------------------------------------------------------------
  assert ("%d tokens total" % EXPECTED_TOTAL) in out, \
    "Total token count changed (expected %d).\nOutput:\n%s" % (EXPECTED_TOTAL, out)

  def check_count(type_name, expected):
    # The summary line looks like:  "  TAG             : 42"
    needle = "%s : %d" % (type_name.ljust(15), expected)
    assert needle in out, \
      "%s count changed (expected %d).\nOutput:\n%s" % (type_name, expected, out)

  check_count("TAG",          EXPECTED_TAG)
  check_count("VALUE",        EXPECTED_VALUE)
  check_count("BLOCK_HEADER", EXPECTED_BLOCK_HEADER)
  check_count("LOOP",         EXPECTED_LOOP)

  # -----------------------------------------------------------------------
  # Spot-check specific token texts — verifies the tokenizer produces the
  # right text for each important token class.
  # -----------------------------------------------------------------------
  checks = [
    ("data_1UBQ",              "block header"),
    ("P 21 21 21",             "single-quoted string (space group)"),
    ("Structure of ubiquitin", "semicolon text field (first line)"),
    ("3.60e-02",               "scientific notation"),
    ("?",                      "unknown value"),
    ("_cell.length_a",         "tag name"),
  ]
  for text, description in checks:
    assert text in out, \
      "Expected %s token %r not found in output.\nOutput:\n%s" % (
        description, text, out)

if __name__ == "__main__":
  exercise()
  print(format_cpu_times())
  print("OK")
