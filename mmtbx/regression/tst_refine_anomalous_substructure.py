
from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out

def exercise():
  from mmtbx.command_line import refine_anomalous_substructure
  from mmtbx.regression.make_fake_anomalous_data import generate_cd_cl_inputs
  pdb_file, mtz_file = generate_cd_cl_inputs(file_base="tst_anom_ref")
  groups = refine_anomalous_substructure.run(
    args=[pdb_file, mtz_file],
    out=null_out())
  assert len(groups) == 2

if (__name__ == "__main__"):
  exercise()
