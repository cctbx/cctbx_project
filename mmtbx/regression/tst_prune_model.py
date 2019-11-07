from __future__ import absolute_import, division, print_function

import libtbx.load_env
from libtbx import easy_run
import os

def exercise():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/prune.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/prune_data.mtz",
    test=os.path.isfile)
  if (None in [pdb_file, mtz_file]):
    print("Input files not found, skipping test")
    return
  cmd = "mmtbx.prune_model \"%s\" \"%s\" mainchain=True" % (pdb_file, mtz_file)
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  if (not "Removed 3 residues and 1 sidechains" in result.stdout_lines):
    raise RuntimeError(("Program output differs from expected - last 20 lines "+
      "shown below:\n\n%s") % "\n".join(result.stdout_lines[-20:]))

if (__name__ == "__main__"):
  exercise()
  print("OK")
