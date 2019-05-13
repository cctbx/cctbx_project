
from __future__ import division
from __future__ import print_function
from cStringIO import StringIO
import os

def exercise():
  from mmtbx.command_line import compute_map_coefficients
  from mmtbx.regression.make_fake_anomalous_data import generate_cd_cl_inputs
  from iotbx import file_reader
  mtz_file, pdb_file = generate_cd_cl_inputs(
    file_base="tst_map_coeffs_simple")
  args = [mtz_file, pdb_file, "map_type=anom_residual",
          "skip_twin_detection=True"]
  out = StringIO()
  compute_map_coefficients.run(args=args, out=out)
  assert """Using wavelength = 1.116 from PDB header""" in out.getvalue()
  assert os.path.isfile("tst_map_coeffs_simple_anom_residual.mtz")
  mtz_in = file_reader.any_file("tst_map_coeffs_simple_anom_residual.mtz")
  assert (mtz_in.file_server.miller_arrays[0].info().label_string() ==
          "ANOM_DIFF,PHANOM_DIFF")
  print("OK")

if (__name__ == "__main__"):
  exercise()
