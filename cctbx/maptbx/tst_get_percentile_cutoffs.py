from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import time
from cctbx import maptbx

pdb_str="""
CRYST1    8.228   11.366   10.991  90.00  90.00  90.00 P 1
SCALE1      0.121536  0.000000  0.000000        0.00000
SCALE2      0.000000  0.087982  0.000000        0.00000
SCALE3      0.000000  0.000000  0.090984        0.00000
ATOM      1  N  AGLY B   1       3.105   3.876   3.794  1.00  1.00           N
ATOM      2  CA AGLY B   1       3.987   4.960   4.314  1.00  1.00           C
ATOM      3  C  AGLY B   1       5.456   4.601   4.162  1.00  1.00           C
ATOM      4  O  AGLY B   1       5.759   3.433   3.921  1.00  1.00           O
ATOM      5  H1 AGLY B   1       2.686   3.651   4.546  1.00  1.00           H
ATOM      6  H2 AGLY B   1       2.814   4.265   3.048  1.00  1.00           H
ATOM      7  H3 AGLY B   1       3.289   3.065   3.477  1.00  1.00           H
ATOM      8  HA2AGLY B   1       3.792   5.776   3.837  1.00  1.00           H
ATOM      9  HA3AGLY B   1       3.796   5.095   5.255  1.00  1.00           H
ATOM     10  N  BGLY B   1       4.736   2.369   3.428  1.00  1.00           N
ATOM     11  CA BGLY B   1       5.850   3.229   3.918  1.00  1.00           C
ATOM     12  C  BGLY B   1       5.409   4.663   4.139  1.00  1.00           C
ATOM     13  O  BGLY B   1       4.213   4.956   4.163  1.00  1.00           O
ATOM     14  H1 BGLY B   1       4.969   1.993   2.655  1.00  1.00           H
ATOM     15  H2 BGLY B   1       4.565   1.734   4.027  1.00  1.00           H
ATOM     16  H3 BGLY B   1       4.008   2.867   3.308  1.00  1.00           H
ATOM     17  HA2BGLY B   1       6.186   2.877   4.757  1.00  1.00           H
ATOM     18  HA3BGLY B   1       6.572   3.225   3.270  1.00  1.00           H
ATOM     19  N   CYS B   2       6.381   5.557   4.294  1.00  1.00           N
ATOM     20  CA  CYS B   2       6.108   6.984   4.532  1.00  1.00           C
ATOM     21  C   CYS B   2       5.166   7.276   5.709  1.00  1.00           C
ATOM     22  O   CYS B   2       5.514   7.018   6.862  1.00  1.00           O
ATOM     23  CB  CYS B   2       5.630   7.674   3.244  1.00  1.00           C
ATOM     24  SG  CYS B   2       5.160   9.406   3.466  1.00  1.00           S
ATOM     25  H   CYS B   2       7.222   5.385   4.257  1.00  1.00           H
ATOM     26  HA  CYS B   2       6.955   7.399   4.759  1.00  1.00           H
ATOM     27  HB2 CYS B   2       6.356   7.653   2.601  1.00  1.00           H
ATOM     28  HB3 CYS B   2       4.866   7.207   2.879  1.00  1.00           H
ATOM     29  N   CYS B   3       3.983   7.810   5.419  1.00  1.00           N
ATOM     30  CA  CYS B   3       3.014   8.128   6.460  1.00  1.00           C
ATOM     31  C   CYS B   3       2.428   6.855   7.061  1.00  1.00           C
ATOM     32  O   CYS B   3       1.975   6.848   8.205  1.00  1.00           O
ATOM     33  CB  CYS B   3       1.892   9.000   5.895  1.00  1.00           C
ATOM     34  SG  CYS B   3       2.443  10.598   5.253  1.00  1.00           S
ATOM     35  H   CYS B   3       3.708   7.997   4.627  1.00  1.00           H
ATOM     36  HA  CYS B   3       3.456   8.623   7.168  1.00  1.00           H
ATOM     37  HB2 CYS B   3       1.463   8.523   5.168  1.00  1.00           H
ATOM     38  HB3 CYS B   3       1.247   9.172   6.599  1.00  1.00           H
TER
ATOM     39  O   HOH C   1       1.194   0.871   7.026  1.00  7.85           O
ATOM     40  O   HOH C   2       2.387   2.996   6.915  1.00  6.79           O
ATOM     41  O   HOH C   3       6.169   4.601   8.563  0.50  5.59           O
ATOM     42  O   HOH C   4       3.645   3.402   9.119  1.00  6.04           O
ATOM     43  O   HOH C   5       4.716   5.971   0.442  1.00  5.95           O
ATOM     44  O   HOH C   6       2.271   7.570   2.088  1.00  5.17           O
ATOM     45  O   HOH C   7       4.683   8.148   9.569  1.00  4.74           O
ATOM     46  O   HOH C   8       7.084   7.769   9.388  1.00  4.51           O
ATOM     47  O   HOH C   9       6.636   9.772   7.633  1.00  5.99           O
TER
"""

def percentile_cutoffs_inefficient(
      map_data,
      vol_cutoff_plus_percent,
      vol_cutoff_minus_percent):
  s = flex.sort_permutation(map_data.as_1d())
  map_data_sorted = map_data.select(s)
  i = map_data.size()-1-int(map_data.size()*(vol_cutoff_plus_percent/100.))
  cutoffp = map_data_sorted[i]
  j = int(map_data.size()*(vol_cutoff_minus_percent/100.))
  cutoffm = map_data_sorted[j]
  return cutoffp, cutoffm

def exercise():
  from iotbx.pdb.utils import get_pdb_input
  pdb_inp = get_pdb_input(text = pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  fc = xrs.structure_factors(d_min=1.5).f_calc()
  fft_map = fc.fft_map(resolution_factor=0.1)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  for vcp in [1,50,99]:
    for vcm in [1,50,99]:
      # Doing it old inefficient way
      po, mo = percentile_cutoffs_inefficient(
        map_data                 = map_data,
        vol_cutoff_plus_percent  = vcp,
        vol_cutoff_minus_percent = vcm)
      # Doing it new memory efficient way
      hist = maptbx.histogram(map_data, n_bins=min(10000, map_data.size()))
      pn, mn = hist.get_percentile_cutoffs(
        map                      = map_data,
        vol_cutoff_plus_percent  = vcp,
        vol_cutoff_minus_percent = vcm)
      assert approx_equal(po,pn, 0.01)
      assert approx_equal(mo,mn, 0.01)
      rp = (map_data >= po).count(True)*100./map_data.size()
      rm = (map_data <= mo).count(True)*100./map_data.size()
      assert approx_equal(rp, vcp, 1.e-3)
      assert approx_equal(rm, vcm, 1.e-3)
      assert approx_equal(
        rp,
        (map_data >= pn).count(True)*100./map_data.size(), 0.1)
      assert approx_equal(
        rm,
        (map_data <= mn).count(True)*100./map_data.size(), 0.1)

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f" % (time.time()-t0))
