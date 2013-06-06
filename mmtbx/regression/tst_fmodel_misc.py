
from __future__ import division
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from libtbx import easy_run
import os

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  from mmtbx.command_line import fmodel
  from iotbx import file_reader
  from scitbx.array_family import flex
  if (os.path.isfile("tst_fmodel_anomalous.mtz")) :
    os.remove("tst_fmodel_anomalous.mtz")
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(
    file_base="tst_fmodel_anomalous")
  # phenix.fmodel (with wavelength)
  args = [
    pdb_file,
    "high_resolution=1.0",
    "wavelength=1.116",
    "label=F",
    "type=real",
    "output.file_name=tst_fmodel_anomalous.mtz",
    "r_free_flags_fraction=0.1",
  ]
  fmodel.run(args=args, log=null_out())
  assert os.path.isfile("tst_fmodel_anomalous.mtz")
  mtz_in = file_reader.any_file("tst_fmodel_anomalous.mtz")
  array = mtz_in.file_server.miller_arrays[0]
  assert (array.anomalous_flag())
  anom_diffs = array.anomalous_differences()
  assert approx_equal(flex.max(anom_diffs.data()), 5.72, eps=0.01)
  # mmtbx.fmodel_simple
  result = easy_run.fully_buffered(
    "mmtbx.fmodel_simple \"%s\" tst_fmodel_anomalous.mtz" % pdb_file
    ).raise_if_errors()
  assert (result.return_code == 0)
  print "OK"

if (__name__ == "__main__") :
  exercise()
