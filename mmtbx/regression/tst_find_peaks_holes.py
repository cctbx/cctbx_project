
import libtbx.load_env
from libtbx.test_utils import contains_substring
from libtbx.utils import null_out
from cStringIO import StringIO
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1yjp_h.pdb",
    test=os.path.isfile)
  mtz_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not available, skipping test."
    return
  from mmtbx.command_line import find_peaks_holes
  out = StringIO()
  peaks_holes = find_peaks_holes.run(
    args=[pdb_file, mtz_file],
    out=out)
  assert contains_substring(out.getvalue(), "  mFo-DFc >  3      :      2")
  assert contains_substring(out.getvalue(), "  mFo-DFc max       :   3.45")

  peaks_holes.save_pdb_file(file_name="%s.pdb" % os.getpid(), log=null_out())
  from iotbx.file_reader import any_file
  pdbh = any_file("%s.pdb" % os.getpid()).file_object.construct_hierarchy()
  assert (len(pdbh.atoms()) == 4)
  assert (pdbh.atoms()[0].b == 3.45)
  # filter by 2fo-fc
  out = StringIO()
  peaks_holes = find_peaks_holes.run(
    args=[pdb_file, mtz_file, "filter_peaks_by_2fofc=1.0"],
    out=out)
  assert contains_substring(out.getvalue(), "  mFo-DFc max       :   None")

if (__name__ == "__main__") :
  exercise()
  print "OK"
