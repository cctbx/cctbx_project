from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import boost_adaptbx.boost.python as bp
asu_map_ext = bp.import_ext("cctbx_asymmetric_map_ext")
import iotbx.pdb
from mmtbx.maps import correlation
from libtbx.test_utils import approx_equal
from cctbx import maptbx

pdb_str="""
CRYST1   60.410   44.706   51.507  90.00  90.00  90.00 P 1
SCALE1      0.016554  0.000000  0.000000        0.00000
SCALE2      0.000000  0.022368  0.000000        0.00000
SCALE3      0.000000  0.000000  0.019415        0.00000
ATOM      1  N   GLY A   1      20.043  24.706  27.193  1.00 20.00           N
ATOM      2  CA  GLY A   1      20.000  24.301  25.742  1.00 20.00           C
ATOM      3  C   GLY A   1      21.037  23.234  25.510  1.00 20.00           C
ATOM      4  O   GLY A   1      21.529  22.615  26.472  1.00 20.00           O
ATOM      5  N   ASN A   2      21.396  23.017  24.246  1.00 20.00           N
ATOM      6  CA  ASN A   2      22.530  22.132  23.922  1.00 20.00           C
ATOM      7  C   ASN A   2      23.811  22.631  24.518  1.00 20.00           C
ATOM      8  O   ASN A   2      24.074  23.836  24.517  1.00 20.00           O
ATOM      9  CB  ASN A   2      22.706  21.975  22.432  1.00 20.00           C
ATOM     10  CG  ASN A   2      21.468  21.436  21.783  1.00 20.00           C
ATOM     11  OD1 ASN A   2      21.027  20.321  22.107  1.00 20.00           O
ATOM     12  ND2 ASN A   2      20.848  22.249  20.922  1.00 20.00           N
ATOM     13  N   ASN A   3      24.614  21.684  24.996  1.00 20.00           N
ATOM     14  CA  ASN A   3      25.859  21.998  25.680  1.00 20.00           C
ATOM     15  C   ASN A   3      27.097  21.426  24.986  1.00 20.00           C
ATOM     16  O   ASN A   3      27.180  20.213  24.739  1.00 20.00           O
ATOM     17  CB  ASN A   3      25.793  21.472  27.133  1.00 20.00           C
ATOM     18  CG  ASN A   3      27.046  21.833  27.952  1.00 20.00           C
ATOM     19  OD1 ASN A   3      27.350  23.019  28.163  1.00 20.00           O
ATOM     20  ND2 ASN A   3      27.781  20.809  28.397  1.00 20.00           N
ATOM     21  N   GLN A   4      28.047  22.322  24.689  1.00 20.00           N
ATOM     22  CA  GLN A   4      29.436  21.982  24.290  1.00 20.00           C
ATOM     23  C   GLN A   4      30.487  22.700  25.179  1.00 20.00           C
ATOM     24  O   GLN A   4      30.599  23.937  25.206  1.00 20.00           O
ATOM     25  CB  GLN A   4      29.708  22.242  22.802  1.00 20.00           C
ATOM     26  CG  GLN A   4      30.996  21.552  22.304  1.00 20.00           C
ATOM     27  CD  GLN A   4      31.556  22.138  21.002  1.00 20.00           C
ATOM     28  OE1 GLN A   4      31.796  23.362  20.901  1.00 20.00           O
ATOM     29  NE2 GLN A   4      31.802  21.255  20.000  1.00 20.00           N
ATOM     30  N   GLN A   5      31.206  21.915  25.962  1.00 20.00           N
ATOM     31  CA  GLN A   5      32.322  22.455  26.731  1.00 20.00           C
ATOM     32  C   GLN A   5      33.646  21.862  26.263  1.00 20.00           C
ATOM     33  O   GLN A   5      33.820  20.640  26.145  1.00 20.00           O
ATOM     34  CB  GLN A   5      32.108  22.277  28.238  1.00 20.00           C
ATOM     35  CG  GLN A   5      30.881  23.044  28.738  1.00 20.00           C
ATOM     36  CD  GLN A   5      30.396  22.508  30.045  1.00 20.00           C
ATOM     37  OE1 GLN A   5      29.826  21.419  30.093  1.00 20.00           O
ATOM     38  NE2 GLN A   5      30.601  23.281  31.130  1.00 20.00           N
ATOM     39  N   ASN A   6      34.566  22.758  25.947  1.00 20.00           N
ATOM     40  CA  ASN A   6      35.883  22.404  25.409  1.00 20.00           C
ATOM     41  C   ASN A   6      36.906  22.855  26.415  1.00 20.00           C
ATOM     42  O   ASN A   6      37.271  24.037  26.465  1.00 20.00           O
ATOM     43  CB  ASN A   6      36.117  23.110  24.084  1.00 20.00           C
ATOM     44  CG  ASN A   6      35.013  22.829  23.094  1.00 20.00           C
ATOM     45  OD1 ASN A   6      34.850  21.698  22.642  1.00 20.00           O
ATOM     46  ND2 ASN A   6      34.247  23.841  22.770  1.00 20.00           N
ATOM     47  N   TYR A   7      37.344  21.911  27.238  1.00 20.00           N
ATOM     48  CA  TYR A   7      38.211  22.238  28.390  1.00 20.00           C
ATOM     49  C   TYR A   7      39.655  22.425  27.976  1.00 20.00           C
ATOM     50  O   TYR A   7      40.093  21.905  26.946  1.00 20.00           O
ATOM     51  CB  TYR A   7      38.113  21.159  29.460  1.00 20.00           C
ATOM     52  CG  TYR A   7      36.717  21.023  29.993  1.00 20.00           C
ATOM     53  CD1 TYR A   7      35.823  20.115  29.418  1.00 20.00           C
ATOM     54  CD2 TYR A   7      36.262  21.850  31.011  1.00 20.00           C
ATOM     55  CE1 TYR A   7      34.532  20.000  29.887  1.00 20.00           C
ATOM     56  CE2 TYR A   7      34.956  21.743  31.507  1.00 20.00           C
ATOM     57  CZ  TYR A   7      34.099  20.823  30.922  1.00 20.00           C
ATOM     58  OH  TYR A   7      32.818  20.683  31.382  1.00 20.00           O
ATOM     59  OXT TYR A   7      40.410  23.093  28.703  1.00 20.00           O
END
"""

def exercise_d99():
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  ph = pdb_inp.construct_hierarchy()
  xrs = ph.extract_xray_structure(crystal_symmetry = pdb_inp.crystal_symmetry())
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell         = xrs.unit_cell(),
    space_group_info  = xrs.space_group_info(),
    resolution_factor = 0.25,
    d_min             = 2.,
    symmetry_flags    = maptbx.use_space_group_symmetry)
  fc = xrs.structure_factors(d_min=2.).f_calc()
  fft_map = fc.fft_map(crystal_gridding=crystal_gridding)
  map = fft_map.real_map_unpadded()
  #
  o = maptbx.d99(map=map, crystal_symmetry=xrs.crystal_symmetry())
  assert approx_equal(o.result.d99, 2.13, 0.03)

def exercise_five_cc():
  pdb_inp = iotbx.pdb.input(source_info=None, lines = pdb_str)
  ph = pdb_inp.construct_hierarchy()
  #ph.write_pdb_file(file_name="m1.pdb")
  xrs = ph.extract_xray_structure(crystal_symmetry = pdb_inp.crystal_symmetry())
  f_calc = xrs.structure_factors(d_min=2.0).f_calc()#.resolution_filter(d_max=6)
  fft_map = f_calc.fft_map(resolution_factor=0.25)
  m1 = fft_map.real_map_unpadded()
  #
  #xrs = xrs.set_b_iso(value=50)
  sc = xrs.sites_cart()
  #sc = sc + flex.vec3_double(sc.size(), [.5,.5,.5])
  xrs = xrs.replace_sites_cart(new_sites=sc)
  ph.adopt_xray_structure(xrs)
  #ph.write_pdb_file(file_name="m2.pdb")
  o = correlation.five_cc(map=m1, xray_structure=xrs, d_min=2.0,
    compute_cc_box=True, compute_cc_image=True, box=False).result
  assert approx_equal(o.cc_box, 1.0)
  assert approx_equal(o.cc_image  , 1.0)
  assert approx_equal(o.cc_mask   , 1.0)
  assert approx_equal(o.cc_peaks  , 1.0)
  assert approx_equal(o.cc_volume , 1.0)

if (__name__ == "__main__"):
  exercise_d99()
  exercise_five_cc()
