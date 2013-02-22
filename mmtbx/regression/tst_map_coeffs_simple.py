
from __future__ import division
from libtbx.utils import null_out
import os

def exercise () :
  if (os.path.isfile("tst_map_coeffs_simple_anom_residual.mtz")) :
    os.remove("tst_map_coeffs_simple_anom_residual.mtz")
  from mmtbx.command_line import compute_map_coefficients
  from mmtbx.regression import make_fake_anomalous_data
  from iotbx import file_reader
  mtz_file, pdb_file = make_fake_anomalous_data.generate_cd_cl_inputs(
    file_base="tst_map_coeffs_simple")
  args = [ mtz_file, pdb_file, "map_type=anom_residual", "wavelength=1.116",
           "skip_twin_detection=True", ]
  compute_map_coefficients.run(args=args, out=null_out())
  assert os.path.isfile("tst_map_coeffs_simple_anom_residual.mtz")
  mtz_in = file_reader.any_file("tst_map_coeffs_simple_anom_residual.mtz")
  assert (mtz_in.file_server.miller_arrays[0].info().label_string() ==
          "ANOM_DIFF,PHANOM_DIFF")
  print "OK"

if (__name__ == "__main__") :
  exercise()
