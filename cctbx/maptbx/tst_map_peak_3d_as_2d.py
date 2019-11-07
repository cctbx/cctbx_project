from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import iotbx.pdb
import mmtbx.real_space
from cctbx import maptbx
from libtbx.test_utils import approx_equal

pdb_str = """
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM  120  CA  CA  A   1       5.000   5.000   5.000  1.00 10.00          CA
TER
END
"""

def run(b_iso=10):
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  xrs = xrs.set_b_iso(value=b_iso)
  xrs.scattering_type_registry(
    table = "n_gaussian",
    d_min = 0.,
    types_without_a_scattering_contribution=["?"])
  n_real = [100,100,100]
  pixel_volume = xrs.unit_cell().volume()/(n_real[0]*n_real[1]*n_real[2])
  map_data_3d = mmtbx.real_space.sampled_model_density(
    xray_structure = xrs,
    n_real         = n_real).data()*pixel_volume
  dist, map_data_2d = maptbx.map_peak_3d_as_2d(
    map_data    = map_data_3d,
    unit_cell   = xrs.unit_cell(),
    center_cart = xrs.sites_cart()[0],
    radius      = 3.0)
  #
  map_data_2d_exact = flex.double()
  ed = xrs._scattering_type_registry.gaussian("Ca")
  for r in dist:
    map_data_2d_exact.append(ed.electron_density(r, b_iso))
  map_data_2d_exact = map_data_2d_exact * pixel_volume
  #
  assert approx_equal(flex.sum(map_data_3d), 20, 0.1) # Page 556, Int.Tables.
  assert flex.mean(abs(map_data_2d-map_data_2d_exact)) < 1.e-4

if (__name__ == "__main__"):
  run()
