
from __future__ import division
from libtbx import easy_run
import time
import os

def exercise () :
  from mmtbx.regression.make_fake_anomalous_data import generate_calcium_inputs
  from iotbx.file_reader import any_file
  anom_mtz_file, pdb_file = generate_calcium_inputs(
    file_base = "tst_isomorphous_difference_misc_anom", anonymize = False)
  hoh_file = generate_calcium_inputs(
    file_base = "tst_isomorphous_difference_misc_hoh", anonymize = True)[1]
  mtz_file = "tst_isomorphous_difference_misc.mtz"
  args = [
    "phenix.fmodel",
    hoh_file,
    "high_resolution=1.5",
    "r_free_flags_fraction=0.1",
    "type=real",
    "output.file_name=" + mtz_file
    ]
  easy_run.fully_buffered(args).raise_if_errors()
  time.sleep(2)
  assert os.path.isfile(mtz_file)
  minus_mtz_file = "tst_isomorphous_difference_misc_anon_minus.mtz"
  args = [
    "phenix.fobs_minus_fobs_map",
    "f_obs_1_file=" + anom_mtz_file,
    "f_obs_2_file=" + mtz_file,
    pdb_file,
    "omit_selection=\"element CA\"",
    "multiscale=True",
    "output_file=" + minus_mtz_file,
  ]
  result = easy_run.fully_buffered(args).raise_if_errors()
  assert os.path.isfile(minus_mtz_file)
  assert ("1 atoms selected for removal" in result.stdout_lines)
  f = any_file(minus_mtz_file)
  coeffs = f.file_server.miller_arrays[0]
  fft_map = coeffs.fft_map(resolution_factor=0.25)
  minus_map = fft_map.apply_sigma_scaling().real_map_unpadded()
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  for i_seq, atom in enumerate(hierarchy.atoms()) :
    if (atom.element == "CA") :
      site_frac = xrs.sites_frac()[i_seq]
      val = minus_map.eight_point_interpolation(site_frac)
      assert (val > 40)
  print "OK"

if (__name__ == "__main__") :
  exercise()
