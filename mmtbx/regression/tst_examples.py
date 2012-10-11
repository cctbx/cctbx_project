
from __future__ import division
from libtbx import easy_run
import libtbx.load_env
import os

def exercise () :
  script_file = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/examples/simple_command_line_cc.py",
    test=os.path.isfile)
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1yjp_h.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not found, skipping"
    return
  args = [script_file, pdb_file, mtz_file]
  result = easy_run.fully_buffered(
    "mmtbx.python %s" % " ".join(args)).raise_if_errors()
  assert ("CC(obs-calc): 0.953" in result.stdout_lines)
  print "OK"

if (__name__ == "__main__") :
  exercise()
