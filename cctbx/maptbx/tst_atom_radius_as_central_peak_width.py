from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from libtbx.test_utils import approx_equal

def run():
  rad = maptbx.atom_radius_as_central_peak_width(
    element="C", b_iso=100, d_min=2.0, scattering_table="electron")
  assert approx_equal(rad, 3.69)
  rad = maptbx.atom_radius_as_central_peak_width(
    element="C", b_iso=100, d_min=2.5, scattering_table="electron")
  assert approx_equal(rad, 3.77)

if (__name__ == "__main__"):
  run()
  print("OK")
