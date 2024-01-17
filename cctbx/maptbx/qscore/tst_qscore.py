from __future__ import absolute_import, division, print_function
import os
from cctbx.array_family import flex
from iotbx.cli_parser import run_program
import libtbx
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
from iotbx.data_manager import DataManager
from cctbx.programs import qscore




expected_results = {

}
# make flex arrays
expected_results = {key:flex.double(val) for key,val in list(expected_results.items())}

def exercise(test_name):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.pdb",
    test=os.path.isfile)

  map_file = libtbx.env.find_in_repositories(
    relative_path=f"phenix_regression/real_space_refine/data/tst_{test_name}.ccp4",
    test=os.path.isfile)

  result = run_program(
    program_class=qscore.Program,
    args = [pdb_file,map_file,"probe_allocation_method=progressive"],
    logger=null_out(),
  )

  
  try:
    expected_result = expected_results[test_name]
    assert approx_equal(expected_result, result.qscore,eps=1.e-2)
  except:
    for val in result.qscore:
      print(str(val)+",")
    raise

if (__name__ == "__main__"):
  """
  Test random files to verify basic functionality remains unchanged
  Data from phenix_regression/real_space_refine/data
  """
  for test_name in [17,42,48]:
    exercise(test_name)
  print("OK")
