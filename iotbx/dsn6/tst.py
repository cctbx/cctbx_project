
from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx.development import debug_utils
from libtbx.math_utils import iceil
import os.path as op
import os
import sys

def exercise_writer(space_group_info, n_sites=100, d_min=1.5):
  if op.isfile("tst_iotbx_dsn6.omap"):
    os.remove("tst_iotbx_dsn6.omap")
  xrs = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=(("O","N","C")*(n_sites//3+1))[:n_sites],
    volume_per_atom=50,
    min_distance=1.5)
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  fft_map = fc.fft_map(resolution_factor=1/3).apply_sigma_scaling()
  n_real = fft_map.n_real()
  n_blocks = iceil(n_real[0]/8) * iceil(n_real[1]/8) * iceil(n_real[2]/8)
  fft_map.as_dsn6_map(file_name="tst_iotbx_dsn6.omap")
  assert op.isfile("tst_iotbx_dsn6.omap")
  size = op.getsize("tst_iotbx_dsn6.omap")
  assert (size == (n_blocks+1)*512)
  os.remove("tst_iotbx_dsn6.omap")
  fft_map.as_dsn6_map(
    file_name="tst_iotbx_dsn6.omap",
    gridding_first=[-5,-5,-5],
    gridding_last=n_real)
  n_blocks = iceil(1+n_real[0]/8) * iceil(1+n_real[1]/8) * iceil(1+n_real[2]/8)
  size = op.getsize("tst_iotbx_dsn6.omap")
  assert (size == (n_blocks+1)*512)

def run_call_back(flags, space_group_info):
  exercise_writer(space_group_info)

if (__name__ == "__main__"):
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
