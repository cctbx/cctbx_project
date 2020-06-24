
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import os

def exercise():
  # XXX can't use libtbx.env.find_in_repositories for this in our current
  # nightly build system
  script_file = os.path.join(abs(libtbx.env.module_dist_paths['mmtbx']),
    "examples", "simple_command_line_cc.py")
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1yjp_h.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  if (None in [pdb_file, mtz_file]):
    print("phenix_regression not found, skipping")
    return
  assert (script_file is not None)
  args = [script_file, pdb_file, mtz_file]
  result = easy_run.fully_buffered("mmtbx.python \"%s\" \"%s\" \"%s\"" %
    (script_file, pdb_file, mtz_file)).raise_if_errors()
  for l in result.stdout_lines:
    print(l)
  assert ("CC(obs-calc): 0.950" in result.stdout_lines)
  print("OK")

if (__name__ == "__main__"):
  exercise()
