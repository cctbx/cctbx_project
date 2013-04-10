
from __future__ import division
from libtbx import easy_run
import os

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  from iotbx.file_reader import any_file
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs(
    "tmp_fobs_minus_fobs")
  assert os.path.isfile("tmp_fobs_minus_fobs.pdb")
  args = [
    "phenix.fmodel",
    pdb_file,
    "high_resolution=1.5",
    "r_free_flags_fraction=0.1",
    "type=real",
    "output.file_name=tmp_fobs_minus_fobs_hoh.mtz",
  ]
  easy_run.fully_buffered(args).raise_if_errors()
  assert os.path.isfile("tmp_fobs_minus_fobs_hoh.mtz")
  args = [
    "phenix.fobs_minus_fobs_map",
    "f_obs_1_file=tmp_fobs_minus_fobs.mtz",
    "f_obs_2_file=tmp_fobs_minus_fobs_hoh.mtz",
    "tmp_fobs_minus_fobs.pdb",
    "omit_selection=\"element CA\"",
    "multiscale=True",
    "output_file=tmp_fobs_minus_fobs_map_coeffs.mtz",
  ]
  result = easy_run.fully_buffered(args).raise_if_errors()
  assert os.path.isfile("tmp_fobs_minus_fobs_map_coeffs.mtz")
  assert ("1 atoms selected for removal" in result.stdout_lines)
  f = any_file("tmp_fobs_minus_fobs_map_coeffs.mtz")
  coeffs = f.file_server.miller_arrays[0]
  fft_map = coeffs.fft_map(resolution_factor=0.25)
  map = fft_map.apply_sigma_scaling().real_map_unpadded()
  pdb_in = any_file("tmp_fobs_minus_fobs.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  for i_seq, atom in enumerate(hierarchy.atoms()) :
    if (atom.element == "CA") :
      site_frac = xrs.sites_frac()[i_seq]
      val = map.eight_point_interpolation(site_frac)
      assert (val > 40)
  print "OK"

if (__name__ == "__main__") :
  exercise()
