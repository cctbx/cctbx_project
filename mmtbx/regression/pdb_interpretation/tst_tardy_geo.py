from __future__ import absolute_import, division, print_function

from mmtbx.regression import model_1yjp
from libtbx.test_utils import assert_lines_in_file
from libtbx import easy_run
import libtbx.load_env

def exercise_01(prefix="tst_tardy_geo"):
  fname = "%s.pdb" % prefix
  with open(fname, 'w') as f:
    f.write(model_1yjp)
  cmd = "mmtbx.pdb_interpretation %s write_geo_files=True write_tardy_geo_files=True" % fname
  assert not easy_run.call(cmd)
  assert_lines_in_file(file_name="%s.pdb.geo" % prefix, lines="Bond restraints: 59")
  with open("%s.pdb.tardy.geo" % prefix, 'r') as f:
    tardy_geo_cont = f.read()
    # print (tardy_geo_cont.find("Bond restraints"))
    # there should not be any bond restraints.
    assert tardy_geo_cont.find("Bond restraints") == -1

if(__name__ == "__main__"):
  if libtbx.env.find_in_repositories(relative_path="chem_data") is None:
    print("Skipping exercise_01(): chem_data directory not available")
  else:
    exercise_01()
    print('OK')
