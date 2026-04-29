from __future__ import division
import iotbx.pdb
import mmtbx.model
from scitbx.matrix import col
from scitbx.array_family import flex

pdb_str = """
CRYST1   12.650   11.560   10.720  90.00  90.00  90.00 P 1
ATOM      1  N   GLY B   1       7.650   6.063   5.176  1.00 10.00           N
ATOM      2  CA  GLY B   1       6.683   5.000   5.000  1.00 10.00           C
ATOM      3  C   GLY B   1       5.326   5.383   5.551  1.00 10.00           C
ATOM      4  O   GLY B   1       5.000   6.560   5.720  1.00 10.00           O
TER
END
"""

def tst_01():
  # Load model
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  # Compute Fourier map (Fcalc)
  fc = xrs.structure_factors(d_min=1).f_calc()

  mm = fc.as_map_manager()

  result = mm.peak_search(max_peaks = 3)

  # expected sites
  sc = xrs.sites_cart()

  # Make sure we have sites and sites_cart ok
  sc_found = flex.vec3_double(result.sites_cart)
  sites_frac_found = flex.vec3_double(result.sites)
  sc_found_2= mm.crystal_symmetry().unit_cell().orthogonalize(sites_frac_found)
  assert (sc_found_2 - sc_found).rms_length() < 0.001

  sc_found_absolute = result.sites_cart_absolute
  assert (sc_found_absolute -
      (sc_found - col(mm.shift_cart()))).rms_length() < 0.001

  # List peaks with their heights and match to expected
  for site_cart, height in zip( result.sites_cart, result.heights):
    one_site_list = flex.vec3_double([site_cart])
    dist, ii,jj = sc.min_distance_between_any_pair_with_id(one_site_list)
    assert dist < 0.2

if(__name__ == "__main__"):
  tst_01()
  print ("OK")
