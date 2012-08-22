from __future__ import division

from libtbx.test_utils import approx_equal
from mmtbx.geometry_restraints import reference_coordinate
import iotbx.pdb
from cctbx.array_family import flex
from cctbx import adp_restraints # import dependency

def simple_pdb () :
  import iotbx.pdb
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
""")
  return pdb_in

def exercise_simple () :
  pdb_in = simple_pdb()
  pdb_hierarchy = pdb_in.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  sites_cart = pdb_hierarchy.atoms().extract_xyz()
  reference_coordinate_proxies = \
    reference_coordinate.build_proxies(
      sites_cart=sites_cart,
      pdb_hierarchy=pdb_hierarchy,
      c_alpha_only=True).reference_coordinate_proxies
  grads = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
  residual = reference_coordinate.target_and_gradients(
               proxies=reference_coordinate_proxies,
               sites_cart=sites_cart,
               gradient_array=grads)
  assert approx_equal(residual, 0.0)

if (__name__ == "__main__") :
  exercise_simple()
  print "OK"
