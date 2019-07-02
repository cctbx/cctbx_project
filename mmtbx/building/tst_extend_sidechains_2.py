from __future__ import absolute_import, division, print_function

from mmtbx.building import extend_sidechains
import mmtbx.monomer_library
from scitbx.array_family import flex
import iotbx.pdb
from libtbx.test_utils import approx_equal

mon_lib_srv = mmtbx.monomer_library.server.server()

###

pdb_str0_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
ATOM      0  HA  LYS A   7       5.646   6.153   2.918  1.00 18.23           H   new
remark ATOM      0  HB2 LYS A   7       3.514   5.375   1.233  1.00 23.56           H   new
remark ATOM      0  HB3 LYS A   7       3.315   6.378   2.411  1.00 23.56           H   new
remark ATOM      0  HG2 LYS A   7       4.599   7.982   1.386  1.00 33.58           H   new
remark ATOM      0  HG3 LYS A   7       5.161   6.961   0.350  1.00 33.58           H   new
remark ATOM      0  HD2 LYS A   7       2.426   7.717   0.393  1.00 41.39           H   new
remark ATOM      0  HD3 LYS A   7       3.456   8.385  -0.570  1.00 41.39           H   new
remark ATOM      0  HE2 LYS A   7       3.826   6.296  -1.622  1.00 48.81           H   new
remark ATOM      0  HE3 LYS A   7       2.797   5.628  -0.660  1.00 48.81           H   new
remark ATOM      0  HZ1 LYS A   7       1.802   6.117  -2.648  1.00 48.81           H   new
remark ATOM      0  HZ2 LYS A   7       1.154   6.908  -1.616  1.00 48.81           H   new
remark ATOM      0  HZ3 LYS A   7       2.113   7.530  -2.513  1.00 48.81           H   new
"""

pdb_str0_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
ATOM      0  HA  LYS A   7       5.646   6.153   2.918  1.00 18.23           H   new
ATOM      0  HB2 LYS A   7       3.514   5.375   1.233  1.00 23.56           H   new
ATOM      0  HB3 LYS A   7       3.315   6.378   2.411  1.00 23.56           H   new
ATOM      0  HG2 LYS A   7       4.599   7.982   1.386  1.00 33.58           H   new
ATOM      0  HG3 LYS A   7       5.161   6.961   0.350  1.00 33.58           H   new
ATOM      0  HD2 LYS A   7       2.426   7.717   0.393  1.00 41.39           H   new
ATOM      0  HD3 LYS A   7       3.456   8.385  -0.570  1.00 41.39           H   new
ATOM      0  HE2 LYS A   7       3.826   6.296  -1.622  1.00 48.81           H   new
ATOM      0  HE3 LYS A   7       2.797   5.628  -0.660  1.00 48.81           H   new
ATOM      0  HZ1 LYS A   7       1.802   6.117  -2.648  1.00 48.81           H   new
ATOM      0  HZ2 LYS A   7       1.154   6.908  -1.616  1.00 48.81           H   new
ATOM      0  HZ3 LYS A   7       2.113   7.530  -2.513  1.00 48.81           H   new
"""

###

pdb_str1_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
remark ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
ATOM      0  HA  LYS A   7       5.646   6.153   2.918  1.00 18.23           H   new
ATOM      0  HB2 LYS A   7       3.514   5.375   1.233  1.00 23.56           H   new
ATOM      0  HB3 LYS A   7       3.315   6.378   2.411  1.00 23.56           H   new
ATOM      0  HG2 LYS A   7       4.599   7.982   1.386  1.00 33.58           H   new
ATOM      0  HG3 LYS A   7       5.161   6.961   0.350  1.00 33.58           H   new
ATOM      0  HD2 LYS A   7       2.426   7.717   0.393  1.00 41.39           H   new
ATOM      0  HD3 LYS A   7       3.456   8.385  -0.570  1.00 41.39           H   new
ATOM      0  HE2 LYS A   7       3.826   6.296  -1.622  1.00 48.81           H   new
ATOM      0  HE3 LYS A   7       2.797   5.628  -0.660  1.00 48.81           H   new
ATOM      0  HZ1 LYS A   7       1.802   6.117  -2.648  1.00 48.81           H   new
ATOM      0  HZ2 LYS A   7       1.154   6.908  -1.616  1.00 48.81           H   new
ATOM      0  HZ3 LYS A   7       2.113   7.530  -2.513  1.00 48.81           H   new
"""

pdb_str1_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
ATOM      0  HA  LYS A   7       5.646   6.153   2.918  1.00 18.23           H   new
ATOM      0  HB2 LYS A   7       3.514   5.375   1.233  1.00 23.56           H   new
ATOM      0  HB3 LYS A   7       3.315   6.378   2.411  1.00 23.56           H   new
ATOM      0  HG2 LYS A   7       4.599   7.982   1.386  1.00 33.58           H   new
ATOM      0  HG3 LYS A   7       5.161   6.961   0.350  1.00 33.58           H   new
ATOM      0  HD2 LYS A   7       2.426   7.717   0.393  1.00 41.39           H   new
ATOM      0  HD3 LYS A   7       3.456   8.385  -0.570  1.00 41.39           H   new
ATOM      0  HE2 LYS A   7       3.826   6.296  -1.622  1.00 48.81           H   new
ATOM      0  HE3 LYS A   7       2.797   5.628  -0.660  1.00 48.81           H   new
ATOM      0  HZ1 LYS A   7       1.802   6.117  -2.648  1.00 48.81           H   new
ATOM      0  HZ2 LYS A   7       1.154   6.908  -1.616  1.00 48.81           H   new
ATOM      0  HZ3 LYS A   7       2.113   7.530  -2.513  1.00 48.81           H   new
"""

###

pdb_str2_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
remark ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

pdb_str2_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

###

pdb_str3_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
remark ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

pdb_str3_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

###

pdb_str4_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
remark ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

pdb_str4_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

###

pdb_str5_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
remark ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

pdb_str5_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

###

pdb_str6_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
remark ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
TER
"""

pdb_str6_good="""\
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
TER
"""

###

pdb_str7_good="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
ATOM      0  HA  LYS A   7       5.646   6.153   2.918  1.00 18.23           H   new
ATOM      0  HB2 LYS A   7       3.514   5.375   1.233  1.00 23.56           H   new
ATOM      0  HB3 LYS A   7       3.315   6.378   2.411  1.00 23.56           H   new
ATOM      0  HG2 LYS A   7       4.599   7.982   1.386  1.00 33.58           H   new
ATOM      0  HG3 LYS A   7       5.161   6.961   0.350  1.00 33.58           H   new
ATOM      0  HD2 LYS A   7       2.426   7.717   0.393  1.00 41.39           H   new
ATOM      0  HD3 LYS A   7       3.456   8.385  -0.570  1.00 41.39           H   new
ATOM      0  HE2 LYS A   7       3.826   6.296  -1.622  1.00 48.81           H   new
ATOM      0  HE3 LYS A   7       2.797   5.628  -0.660  1.00 48.81           H   new
ATOM      0  HZ1 LYS A   7       1.802   6.117  -2.648  1.00 48.81           H   new
ATOM      0  HZ2 LYS A   7       1.154   6.908  -1.616  1.00 48.81           H   new
ATOM      0  HZ3 LYS A   7       2.113   7.530  -2.513  1.00 48.81           H   new
"""

pdb_str7_bad="""
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.368   7.206   0.851  1.00 33.58           C
remark ATOM      7  CD  LYS A   7       3.242   7.559  -0.108  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.009   6.453  -1.124  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.909   6.785  -2.070  1.00 48.81           N
"""

pdb_str_do_no_move_if_complete="""
ATOM      1  N   LYS A   7       6.187   5.271   2.497  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.198   6.138   3.127  1.00 18.23           C
ATOM      3  C   LYS A   7       4.507   5.426   4.285  1.00 14.78           C
ATOM      4  O   LYS A   7       3.859   4.396   4.096  1.00 15.83           O
ATOM      5  CB  LYS A   7       4.161   6.608   2.104  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.553   7.868   1.345  1.00 33.58           C
ATOM      7  CD  LYS A   7       5.574   7.573   0.257  1.00 41.39           C
ATOM      8  CE  LYS A   7       4.958   6.765  -0.874  1.00 48.81           C
ATOM      9  NZ  LYS A   7       5.929   6.522  -1.976  1.00 48.81           N
TER
END
"""

###

def check(answer, result, bad):
  a = list(answer.atoms().extract_name())
  a.sort()
  #
  r = list(result.atoms().extract_name())
  r.sort()
  #
  b = list(bad.atoms().extract_name())
  b.sort()
  #
  diff = list(set(r).symmetric_difference(set(a)))
  assert len(diff) in [0,1], diff
  if(len(diff)==1):
    assert diff[0]==' H  ', diff
  assert b!=a
  #
  xyz_a = flex.vec3_double()
  xyz_r = flex.vec3_double()
  for aa in answer.atoms():
    for ar in result.atoms():
      if(aa.name == ar.name):
        xyz_a.append(aa.xyz)
        xyz_r.append(ar.xyz)
  print(flex.max(flex.sqrt((xyz_a - xyz_r).dot())))

def exercise_extend_sidechains(pdb_str_bad, pdb_str_good, i, Sorry_msg, add_h):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_bad)
  pdb_h = pdb_inp.construct_hierarchy()
  pdb_h.write_pdb_file(file_name="m_in_%d.pdb"%i)
  pdb_h_bad = pdb_h.deep_copy()
  e = None
  try:
    extend_sidechains.extend_protein_model(
      pdb_hierarchy = pdb_h,
      mon_lib_srv   = mon_lib_srv,
      add_hydrogens = add_h)
  except Exception as e:
    if(Sorry_msg is not None):
      assert str(e).find(Sorry_msg)>-1
      return
  #
  pdb_h.write_pdb_file(file_name="m_completed_%d.pdb"%i)
  pdb_h_result = pdb_h.deep_copy()
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str_good)
  pdb_h_answer = pdb_inp.construct_hierarchy()
  pdb_h_answer.write_pdb_file(file_name="m_good_%d.pdb"%i)
  #
  check(answer=pdb_h_answer, result=pdb_h_result, bad=pdb_h_bad)

def exercise_do_not_move_if_complete():
  pdb_inp = iotbx.pdb.input(source_info=None,
    lines=pdb_str_do_no_move_if_complete)
  pdb_h = pdb_inp.construct_hierarchy()
  xyz1 = pdb_h.atoms().extract_xyz()
  pdb_h_completed = pdb_h.deep_copy()
  extend_sidechains.extend_protein_model(
    pdb_hierarchy = pdb_h_completed,
    mon_lib_srv   = mon_lib_srv,
    add_hydrogens = False)
  xyz2 = pdb_h_completed.atoms().extract_xyz()
  assert approx_equal(flex.max(flex.sqrt((xyz1 - xyz2).dot())),0.0)

if(__name__ == "__main__"):
  exercise_do_not_move_if_complete()
  #
  for i, t in enumerate([
     (pdb_str0_bad, pdb_str0_good, None, None),
     (pdb_str1_bad, pdb_str1_good, None, None),
     (pdb_str2_bad, pdb_str2_good, None, None),
     (pdb_str3_bad, pdb_str3_good, None, None),
     (pdb_str4_bad, pdb_str4_good, "Main chain must be complete.", None),
     (pdb_str5_bad, pdb_str5_good, "Main chain must be complete.", None),
     (pdb_str6_bad, pdb_str6_good, None, None),
     (pdb_str7_bad, pdb_str7_good, None, True),
    ]):
    exercise_extend_sidechains(
      pdb_str_bad=t[0], pdb_str_good=t[1], i=i, Sorry_msg=t[2], add_h=t[3])
  print("OK")
