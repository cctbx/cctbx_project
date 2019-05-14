
from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out

def exercise():
  from mmtbx.command_line import refine_anomalous_substructure
  from mmtbx.regression.make_fake_anomalous_data import generate_cd_cl_inputs
  mtz_file, pdb_file = generate_cd_cl_inputs(
    file_base = "tst_mmtbx_substructure")
  # XXX for some reason even though I'm using synthetic data, it ends up with
  # an extra water with a peak height of 3.1
  anom_groups = refine_anomalous_substructure.run(
    args=[mtz_file, pdb_file, "skip_twin_detection=True", "map_sigma_min=4"],
    out=null_out())
  assert (len(anom_groups) == 2)
  print("OK")

if (__name__ == "__main__"):
  exercise()
