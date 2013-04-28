
from __future__ import division
from libtbx.utils import null_out

def exercise () :
  from mmtbx.building import extend_sidechains
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string="""
ATOM     65  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM     66  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM     67  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM     68  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM     69  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM     70  CG  LYS A   7       4.345   7.215   0.830  1.00 33.58           C
ATOM     71  CD  LYS A   7       3.213   7.570  -0.123  1.00 41.39           C
ATOM     72  CE  LYS A   7       2.976   6.471  -1.165  1.00 48.81           C
""")
  extend_sidechains.extend_protein_model(
    pdb_hierarchy=pdb_in.hierarchy,
    modify_segids=False,
    log=null_out())
  assert (pdb_in.hierarchy.as_pdb_string() == """\
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      5  O   LYS A   7       5.430   3.767   4.222  1.00 15.83           O
ATOM      6  CG  LYS A   7       4.373   7.174   0.814  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.240   7.520  -0.138  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.007   6.412  -1.152  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.899   6.737  -2.090  1.00 48.81           N
TER
""")
  print "OK"

if (__name__ == "__main__") :
  exercise()
