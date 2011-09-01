
import libtbx.load_env
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
import os

def exercise () :
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/1yjp.mtz",
    test=os.path.isfile)
  if (hkl_file is None) :
    print "Skipping"
    return
  from iotbx.command_line import reindex
  from iotbx import file_reader
  output_file = reindex.run(args=[hkl_file, "change_of_basis=c,b,a",
    "output_file=tmp666.mtz"], out=null_out())
  assert os.path.isfile(output_file)
  arrays_in = file_reader.any_file(hkl_file).file_server.miller_arrays
  arrays_out = file_reader.any_file(output_file).file_server.miller_arrays
  assert (approx_equal(arrays_in[0].unit_cell().parameters(),
    (21.937, 4.866, 23.477, 90.0, 107.08, 90.0), eps=0.001))
  assert (approx_equal(arrays_out[0].unit_cell().parameters(),
    (23.477, 4.866, 21.937, 90.0, 107.08, 90.0), eps=0.001))
  assert (arrays_out[0].info().labels == ["FOBS_X","SIGFOBS_X"])

if (__name__ == "__main__") :
  exercise()
  print "OK"
