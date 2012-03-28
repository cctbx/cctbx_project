
import libtbx.load_env
from libtbx import easy_run
import os

def exercise ():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/prune.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/prune_data.mtz",
    test=os.path.isfile)
  if (None in [pdb_file, mtz_file]) :
    print "Input files not found, skipping test"
    return
  result = easy_run.fully_buffered("mmtbx.prune_model %s %s" % (pdb_file,
    mtz_file))
  assert ("Removed 3 residues and 1 sidechains" in result.stdout_lines)

if (__name__ == "__main__") :
  exercise()
  print "OK"
