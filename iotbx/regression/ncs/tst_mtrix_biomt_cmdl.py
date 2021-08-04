from __future__ import absolute_import, division, print_function
import libtbx
from libtbx import easy_run
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff

pdb_str_0 = """\
REMARK 350   BIOMT1   1  1.000000 -0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1 -0.000000  1.000000 -0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000 -0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -0.809017 -0.500000  0.309017        0.00000
REMARK 350   BIOMT2   2 -0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   2  0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT1   3 -0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT2   3 -0.000000 -0.000000 -1.000000        0.00000
REMARK 350   BIOMT3   3 -1.000000  0.000000 -0.000000        0.00000
REMARK 350   BIOMT1   4  0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   4 -0.500000 -0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   4  0.309017  0.809017 -0.500000       -0.00000
REMARK 350   BIOMT1   5  0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT2   5 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT3   5 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT1   6 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT2   6 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT3   6  0.500000  0.309017 -0.809017       -0.00000
REMARK 350   BIOMT1   7 -0.809017  0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   7  0.500000  0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   7 -0.309017 -0.809017 -0.500000        0.00000
REMARK 350   BIOMT1   8 -0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT2   8  0.500000 -0.309017 -0.809017        0.00000
REMARK 350   BIOMT3   8  0.309017 -0.809017  0.500000        0.00000
REMARK 350   BIOMT1   9 -0.309017  0.809017 -0.500000        0.00000
REMARK 350   BIOMT2   9 -0.809017 -0.500000 -0.309017        0.00000
REMARK 350   BIOMT3   9 -0.500000  0.309017  0.809017        0.00000
HELIX    1   1 PRO A   28  LEU A   33  1                                   6
HELIX    7   7 ASP B   74  ALA B   79  1                                   6
HELIX   12  12 LYS C  151  TYR C  159  1                                   9
SHEET    1   A 2 PRO A  28  LEU A  33  0
SHEET    2   A 2 ASP B  74  ALA B  79 -1  O  ASP B  74   N  PRO A  28
CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1
ATOM     62  N   PRO A  28      33.884  22.417 110.984 14.01 14.01           N
ATOM     63  CA  PRO A  28      34.083  23.080 112.268 12.01 12.01           C
ATOM     64  C   PRO A  28      32.908  24.004 112.604 12.01 12.01           C
ATOM     65  O   PRO A  28      32.254  24.556 111.723 16.00 16.00           O
ATOM     66  CB  PRO A  28      35.381  23.873 112.096 12.01 12.01           C
ATOM     67  CG  PRO A  28      35.369  24.249 110.611 12.01 12.01           C
ATOM     68  CD  PRO A  28      34.692  23.045 109.949 12.01 12.01           C
ATOM     69  N   ASN A  29      32.694  24.235 113.899 14.01 14.01           N
ATOM     70  CA  ASN A  29      31.700  25.147 114.485 12.01 12.01           C
ATOM     71  C   ASN A  29      30.215  24.812 114.276 12.01 12.01           C
ATOM     72  O   ASN A  29      29.406  25.211 115.111 16.00 16.00           O
ATOM     73  CB  ASN A  29      32.013  26.616 114.122 12.01 12.01           C
ATOM     74  CG  ASN A  29      33.381  27.082 114.575 12.01 12.01           C
ATOM     75  OD1 ASN A  29      34.061  26.468 115.379 16.00 16.00           O
ATOM     76  ND2 ASN A  29      33.841  28.204 114.073 14.01 14.01           N
ATOM     77  N   GLU A  30      30.550  21.845 114.346 14.01 14.01           N
ATOM     78  CA  GLU A  30      29.386  21.544 115.199 12.01 12.01           C
ATOM     79  C   GLU A  30      29.802  20.720 116.428 12.01 12.01           C
ATOM     80  O   GLU A  30      29.237  19.675 116.751 16.00 16.00           O
ATOM     81  CB  GLU A  30      28.306  20.825 114.365 12.01 12.01           C
ATOM     82  CG  GLU A  30      27.795  21.680 113.200 12.01 12.01           C
ATOM     83  CD  GLU A  30      26.582  21.020 112.532 12.01 12.01           C
ATOM     84  OE1 GLU A  30      26.801  20.115 111.689 16.00 16.00           O
ATOM     85  OE2 GLU A  30      25.450  21.412 112.878 16.00 16.00           O
ATOM     86  N   LEU A  31      29.828  23.005 116.792 14.01 14.01           N
ATOM     87  CA  LEU A  31      29.299  22.073 115.816 12.01 12.01           C
ATOM     88  C   LEU A  31      28.149  21.279 116.426 12.01 12.01           C
ATOM     89  O   LEU A  31      28.305  20.133 116.853 16.00 16.00           O
ATOM     90  CB  LEU A  31      30.413  21.167 115.268 12.01 12.01           C
ATOM     91  CG  LEU A  31      31.519  21.894 114.483 12.01 12.01           C
ATOM     92  CD1 LEU A  31      32.517  20.866 113.958 12.01 12.01           C
ATOM     93  CD2 LEU A  31      30.995  22.696 113.286 12.01 12.01           C
ATOM     94  N   GLY A  32      26.947  21.837 116.303 14.01 14.01           N
ATOM     95  CA  GLY A  32      25.679  21.115 116.378 12.01 12.01           C
ATOM     96  C   GLY A  32      25.506  20.004 115.323 12.01 12.01           C
ATOM     97  O   GLY A  32      24.382  19.638 114.991 16.00 16.00           O
ATOM     98  N   LEU A  33      26.624  19.470 114.815 14.01 14.01           N
ATOM     99  CA  LEU A  33      26.791  18.378 113.866 12.01 12.01           C
ATOM    100  C   LEU A  33      27.494  17.166 114.509 12.01 12.01           C
ATOM    101  O   LEU A  33      27.549  16.110 113.890 16.00 16.00           O
ATOM    102  CB  LEU A  33      27.595  18.877 112.644 12.01 12.01           C
ATOM    103  CG  LEU A  33      27.192  20.239 112.037 12.01 12.01           C
ATOM    104  CD1 LEU A  33      28.098  20.578 110.853 12.01 12.01           C
ATOM    105  CD2 LEU A  33      25.745  20.281 111.549 12.01 12.01           C
ATOM   2090  N   ASP B  74      21.209  34.482 113.075 14.01 14.01           N
ATOM   2091  CA  ASP B  74      20.506  33.939 114.233 12.01 12.01           C
ATOM   2092  C   ASP B  74      21.375  33.841 115.474 12.01 12.01           C
ATOM   2093  O   ASP B  74      22.182  32.923 115.611 16.00 16.00           O
ATOM   2094  CB  ASP B  74      19.906  32.571 113.907 12.01 12.01           C
ATOM   2095  CG  ASP B  74      18.979  32.068 114.997 12.01 12.01           C
ATOM   2096  OD1 ASP B  74      19.110  32.533 116.149 16.00 16.00           O
ATOM   2097  OD2 ASP B  74      18.121  31.208 114.708 16.00 16.00           O
ATOM   2098  N   MET B  75      21.194  34.811 116.364 14.01 14.01           N
ATOM   2099  CA  MET B  75      21.934  34.884 117.610 12.01 12.01           C
ATOM   2100  C   MET B  75      22.089  33.512 118.239 12.01 12.01           C
ATOM   2101  O   MET B  75      23.123  33.213 118.837 16.00 16.00           O
ATOM   2102  CB  MET B  75      21.234  35.824 118.583 12.01 12.01           C
ATOM   2103  CG  MET B  75      21.525  37.288 118.322 12.01 12.01           C
ATOM   2104  SD  MET B  75      23.242  37.703 118.659 32.07 32.07           S
ATOM   2105  CE  MET B  75      23.058  38.659 120.162 12.01 12.01           C
ATOM   2106  N   ILE B  76      21.067  32.672 118.091 14.01 14.01           N
ATOM   2107  CA  ILE B  76      21.116  31.324 118.645 12.01 12.01           C
ATOM   2108  C   ILE B  76      22.422  30.683 118.212 12.01 12.01           C
ATOM   2109  O   ILE B  76      22.994  29.858 118.922 16.00 16.00           O
ATOM   2110  CB  ILE B  76      19.929  30.457 118.191 12.01 12.01           C
ATOM   2111  CG1 ILE B  76      20.068  30.070 116.720 12.01 12.01           C
ATOM   2112  CG2 ILE B  76      18.618  31.182 118.436 12.01 12.01           C
ATOM   2113  CD1 ILE B  76      18.979  29.144 116.230 12.01 12.01           C
ATOM   2114  N   LYS B  77      22.902  31.098 117.044 14.01 14.01           N
ATOM   2115  CA  LYS B  77      24.160  30.591 116.525 12.01 12.01           C
ATOM   2116  C   LYS B  77      25.209  30.734 117.613 12.01 12.01           C
ATOM   2117  O   LYS B  77      26.219  30.038 117.612 16.00 16.00           O
ATOM   2118  CB  LYS B  77      24.577  31.362 115.281 12.01 12.01           C
ATOM   2119  CG  LYS B  77      23.624  31.194 114.112 12.01 12.01           C
ATOM   2120  CD  LYS B  77      24.208  30.282 113.048 12.01 12.01           C
ATOM   2121  CE  LYS B  77      25.656  30.627 112.741 12.01 12.01           C
ATOM   2122  NZ  LYS B  77      25.986  32.047 113.045 14.01 14.01           N
ATOM   2123  N   ILE B  78      24.954  31.644 118.547 14.01 14.01           N
ATOM   2124  CA  ILE B  78      25.852  31.863 119.664 12.01 12.01           C
ATOM   2125  C   ILE B  78      25.432  30.891 120.752 12.01 12.01           C
ATOM   2126  O   ILE B  78      26.264  30.367 121.489 16.00 16.00           O
ATOM   2127  CB  ILE B  78      25.788  33.314 120.171 12.01 12.01           C
ATOM   2128  CG1 ILE B  78      26.541  34.227 119.215 12.01 12.01           C
ATOM   2129  CG2 ILE B  78      26.405  33.447 121.552 12.01 12.01           C
ATOM   2130  CD1 ILE B  78      27.879  33.670 118.785 12.01 12.01           C
ATOM   2131  N   ALA B  79      24.131  30.644 120.844 14.01 14.01           N
ATOM   2132  CA  ALA B  79      23.615  29.721 121.839 12.01 12.01           C
ATOM   2133  C   ALA B  79      24.401  28.426 121.758 12.01 12.01           C
ATOM   2134  O   ALA B  79      25.260  28.159 122.591 16.00 16.00           O
ATOM   2135  CB  ALA B  79      22.140  29.461 121.605 12.01 12.01           C
ATOM   4569  N   LYS C 151      65.128  35.635 115.156 14.01 14.01           N
ATOM   4570  CA  LYS C 151      65.530  35.389 113.777 12.01 12.01           C
ATOM   4571  C   LYS C 151      64.991  36.467 112.836 12.01 12.01           C
ATOM   4572  O   LYS C 151      65.493  37.600 112.882 16.00 16.00           O
ATOM   4573  CB  LYS C 151      65.144  33.970 113.358 12.01 12.01           C
ATOM   4574  CG  LYS C 151      65.912  32.877 114.088 12.01 12.01           C
ATOM   4575  CD  LYS C 151      67.128  32.399 113.304 12.01 12.01           C
ATOM   4576  CE  LYS C 151      67.638  31.055 113.817 12.01 12.01           C
ATOM   4577  NZ  LYS C 151      68.536  30.367 112.845 14.01 14.01           N
ATOM   4578  N   PRO C 152      63.957  36.119 111.972 14.01 14.01           N
ATOM   4579  CA  PRO C 152      63.508  37.232 111.112 12.01 12.01           C
ATOM   4580  C   PRO C 152      63.109  38.448 111.934 12.01 12.01           C
ATOM   4581  O   PRO C 152      63.440  39.574 111.564 16.00 16.00           O
ATOM   4582  CB  PRO C 152      62.297  36.661 110.371 12.01 12.01           C
ATOM   4583  CG  PRO C 152      62.611  35.226 110.237 12.01 12.01           C
ATOM   4584  CD  PRO C 152      63.111  34.901 111.606 12.01 12.01           C
ATOM   4585  N   ALA C 153      62.412  38.211 113.040 14.01 14.01           N
ATOM   4586  CA  ALA C 153      61.963  39.291 113.912 12.01 12.01           C
ATOM   4587  C   ALA C 153      63.098  40.266 114.150 12.01 12.01           C
ATOM   4588  O   ALA C 153      63.069  41.394 113.658 16.00 16.00           O
ATOM   4589  CB  ALA C 153      61.447  38.743 115.234 12.01 12.01           C
ATOM   4590  N   TYR C 154      64.115  39.823 114.883 14.01 14.01           N
ATOM   4591  CA  TYR C 154      65.267  40.665 115.165 12.01 12.01           C
ATOM   4592  C   TYR C 154      65.756  41.261 113.864 12.01 12.01           C
ATOM   4593  O   TYR C 154      66.400  42.302 113.847 16.00 16.00           O
ATOM   4594  CB  TYR C 154      66.376  39.866 115.831 12.01 12.01           C
ATOM   4595  CG  TYR C 154      67.016  40.606 116.973 12.01 12.01           C
ATOM   4596  CD2 TYR C 154      67.159  40.015 118.220 12.01 12.01           C
ATOM   4597  CD1 TYR C 154      67.468  41.902 116.804 12.01 12.01           C
ATOM   4598  CE2 TYR C 154      67.743  40.697 119.266 12.01 12.01           C
ATOM   4599  CE1 TYR C 154      68.055  42.593 117.842 12.01 12.01           C
ATOM   4600  CZ  TYR C 154      68.190  41.987 119.072 12.01 12.01           C
ATOM   4601  OH  TYR C 154      68.775  42.675 120.111 16.00 16.00           O
ATOM   4602  N   GLY C 155      65.432  40.583 112.771 14.01 14.01           N
ATOM   4603  CA  GLY C 155      65.788  41.039 111.449 12.01 12.01           C
ATOM   4604  C   GLY C 155      65.033  42.277 111.006 12.01 12.01           C
ATOM   4605  O   GLY C 155      64.881  43.241 111.751 16.00 16.00           O
ATOM   4606  N   GLN C 156      64.550  42.233 109.773 14.01 14.01           N
ATOM   4607  CA  GLN C 156      63.853  43.353 109.162 12.01 12.01           C
ATOM   4608  C   GLN C 156      62.850  43.968 110.123 12.01 12.01           C
ATOM   4609  O   GLN C 156      62.506  45.141 110.012 16.00 16.00           O
ATOM   4610  CB  GLN C 156      63.138  42.899 107.894 12.01 12.01           C
ATOM   4611  CG  GLN C 156      64.021  42.114 106.945 12.01 12.01           C
ATOM   4612  CD  GLN C 156      63.662  40.645 106.916 12.01 12.01           C
ATOM   4613  OE1 GLN C 156      62.516  40.271 107.160 16.00 16.00           O
ATOM   4614  NE2 GLN C 156      64.642  39.802 106.616 14.01 14.01           N
ATOM   4615  N   ALA C 157      62.387  43.173 111.075 14.01 14.01           N
ATOM   4616  CA  ALA C 157      61.412  43.650 112.046 12.01 12.01           C
ATOM   4617  C   ALA C 157      61.892  44.877 112.817 12.01 12.01           C
ATOM   4618  O   ALA C 157      61.110  45.516 113.514 16.00 16.00           O
ATOM   4619  CB  ALA C 157      61.039  42.537 113.014 12.01 12.01           C
ATOM   4620  N   LYS C 158      63.168  45.211 112.703 14.01 14.01           N
ATOM   4621  CA  LYS C 158      63.684  46.365 113.420 12.01 12.01           C
ATOM   4622  C   LYS C 158      63.183  47.677 112.829 12.01 12.01           C
ATOM   4623  O   LYS C 158      62.833  48.598 113.561 16.00 16.00           O
ATOM   4624  CB  LYS C 158      65.215  46.350 113.441 12.01 12.01           C
ATOM   4625  CG  LYS C 158      65.885  46.986 112.230 12.01 12.01           C
ATOM   4626  CD  LYS C 158      67.364  47.232 112.484 12.01 12.01           C
ATOM   4627  CE  LYS C 158      67.585  48.021 113.764 12.01 12.01           C
ATOM   4628  NZ  LYS C 158      67.696  47.123 114.945 14.01 14.01           N
ATOM   4629  N   TYR C 159      63.215  47.748 111.504 14.01 14.01           N
ATOM   4630  CA  TYR C 159      62.848  48.954 110.776 12.01 12.01           C
ATOM   4631  C   TYR C 159      61.362  49.227 110.625 12.01 12.01           C
ATOM   4632  O   TYR C 159      60.808  50.102 111.288 16.00 16.00           O
ATOM   4633  CB  TYR C 159      63.486  48.914 109.383 12.01 12.01           C
ATOM   4634  CG  TYR C 159      64.927  48.468 109.382 12.01 12.01           C
ATOM   4635  CD2 TYR C 159      65.676  48.472 110.547 12.01 12.01           C
ATOM   4636  CD1 TYR C 159      65.540  48.048 108.213 12.01 12.01           C
ATOM   4637  CE2 TYR C 159      66.996  48.072 110.545 12.01 12.01           C
ATOM   4638  CE1 TYR C 159      66.858  47.646 108.202 12.01 12.01           C
ATOM   4639  CZ  TYR C 159      67.578  47.660 109.370 12.01 12.01           C
ATOM   4640  OH  TYR C 159      68.888  47.260 109.362 16.00 16.00           O
END
"""

pdb_str_1="""
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
CRYST1   11.101   13.080   14.163  90.00  90.00  90.00 P 1
SCALE1      0.090082  0.000000  0.000000        0.00000
SCALE2      0.000000  0.076453  0.000000        0.00000
SCALE3      0.000000  0.000000  0.070607        0.00000
ATOM      1  N   THR A   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THR A   2       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THR A   3       5.075   5.694   8.249  1.00 20.00           C
TER
ATOM      4  O   THR B   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THR B   5       6.101   6.421   6.092  1.00 20.00           C
TER
ATOM      6  OG1 THR A   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THR A   7       5.964   5.000   5.565  1.00 20.00           C
TER
END
"""

pdb_str_1a="""
CRYST1   11.101   13.080   14.163  90.00  90.00  90.00 P 1
SCALE1      0.090082  0.000000  0.000000        0.00000
SCALE2      0.000000  0.076453  0.000000        0.00000
SCALE3      0.000000  0.000000  0.070607        0.00000
ATOM      1  N   THRA1   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THRA1   2       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THRA1   3       5.075   5.694   8.249  1.00 20.00           C
TER
ATOM      4  O   THRB1   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THRB1   5       6.101   6.421   6.092  1.00 20.00           C
TER
ATOM      6  OG1 THRA1   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THRA1   7       5.964   5.000   5.565  1.00 20.00           C
TER
ATOM      1  N   THRA2   1       1.790   4.994  11.032  1.00 20.00           N
ATOM      2  CA  THRA2   2       2.306   4.556   9.740  1.00 20.00           C
ATOM      3  C   THRA2   3       3.660   3.870   9.892  1.00 20.00           C
TER
ATOM      4  O   THRB2   4       4.517   4.327  10.648  1.00 20.00           O
ATOM      5  CB  THRB2   5       2.445   5.734   8.759  1.00 20.00           C
TER
ATOM      6  OG1 THRA2   6       1.166   6.348   8.560  1.00 20.00           O
ATOM      7  CG2 THRA2   7       2.985   5.251   7.420  1.00 20.00           C
TER
ATOM      1  N   THRA3   1       4.100  -0.147  11.534  1.00 20.00           N
ATOM      2  CA  THRA3   2       3.886   0.555  10.273  1.00 20.00           C
ATOM      3  C   THRA3   3       5.088   1.422   9.915  1.00 20.00           C
TER
ATOM      4  O   THRB3   4       5.659   2.093  10.775  1.00 20.00           O
ATOM      5  CB  THRB3   5       2.625   1.437  10.325  1.00 20.00           C
TER
ATOM      6  OG1 THRA3   6       1.479   0.622  10.600  1.00 20.00           O
ATOM      7  CG2 THRA3   7       2.424   2.158   9.000  1.00 20.00           C
TER
END
"""

pdb_str_2="""
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THR A   2       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THR A   3       5.075   5.694   8.249  1.00 20.00           C
TER
ATOM      4  O   THR B   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THR B   5       6.101   6.421   6.092  1.00 20.00           C
TER
ATOM      6  OG1 THR A   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THR A   7       5.964   5.000   5.565  1.00 20.00           C
TER
END
"""

pdb_str_2a="""
ATOM      1  N   THRA1   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THRA1   2       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THRA1   3       5.075   5.694   8.249  1.00 20.00           C
TER
ATOM      4  O   THRB1   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THRB1   5       6.101   6.421   6.092  1.00 20.00           C
TER
ATOM      6  OG1 THRA1   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THRA1   7       5.964   5.000   5.565  1.00 20.00           C
TER
ATOM      1  N   THRA2   1       1.790   4.994  11.032  1.00 20.00           N
ATOM      2  CA  THRA2   2       2.306   4.556   9.740  1.00 20.00           C
ATOM      3  C   THRA2   3       3.660   3.870   9.892  1.00 20.00           C
TER
ATOM      4  O   THRB2   4       4.517   4.327  10.648  1.00 20.00           O
ATOM      5  CB  THRB2   5       2.445   5.734   8.759  1.00 20.00           C
TER
ATOM      6  OG1 THRA2   6       1.166   6.348   8.560  1.00 20.00           O
ATOM      7  CG2 THRA2   7       2.985   5.251   7.420  1.00 20.00           C
TER
ATOM      1  N   THRA3   1       4.100  -0.147  11.534  1.00 20.00           N
ATOM      2  CA  THRA3   2       3.886   0.555  10.273  1.00 20.00           C
ATOM      3  C   THRA3   3       5.088   1.422   9.915  1.00 20.00           C
TER
ATOM      4  O   THRB3   4       5.659   2.093  10.775  1.00 20.00           O
ATOM      5  CB  THRB3   5       2.625   1.437  10.325  1.00 20.00           C
TER
ATOM      6  OG1 THRA3   6       1.479   0.622  10.600  1.00 20.00           O
ATOM      7  CG2 THRA3   7       2.424   2.158   9.000  1.00 20.00           C
TER
END
"""

pdb_str_2b="""
HELIX    1   1 THR A    1  THR A    2  1                                   6
SHEET    1   A 2 THR A   1  THR A   3  0
SHEET    2   A 2 THR B   4  THR B   5 -1  O  THR B   4   N  THR A   2
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THR A   2       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THR A   3       5.075   5.694   8.249  1.00 20.00           C
TER
ATOM      4  O   THR B   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THR B   5       6.101   6.421   6.092  1.00 20.00           C
TER
ATOM      6  OG1 THR A   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THR A   7       5.964   5.000   5.565  1.00 20.00           C
TER
END
"""

pdb_str_3="""
data_1A37
#
loop_
_database_PDB_rev_record.rev_num
_database_PDB_rev_record.type
_database_PDB_rev_record.details
2 SOURCE ?
2 COMPND ?
2 REMARK ?
2 SEQRES ?
2 KEYWDS ?
2 HEADER ?
3 VERSN  ?
4 MTRIX1 ?
4 MTRIX2 ?
4 MTRIX3 ?
#
_cell.entry_id           1A37
_cell.length_a           94.730
_cell.length_b           94.730
_cell.length_c           250.870
_cell.angle_alpha        90.00
_cell.angle_beta         90.00
_cell.angle_gamma        120.00
_cell.Z_PDB              12
#
_symmetry.entry_id                         1A37
_symmetry.space_group_name_H-M             'P 65'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                ?
_symmetry.space_group_name_Hall            ?
#
loop_
_struct_ncs_oper.id
_struct_ncs_oper.code
_struct_ncs_oper.details
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.vector[3]
1 given    ? 1.000000  0.000000 0.000000  0.000000  1.000000  0.000000  0.000000  0.000000  1.000000 0.00000  0.00000  0.00000
2 generate ? -0.997443 0.000760 -0.071468 -0.000162 -0.999965 -0.008376 -0.071472 -0.008343 0.997408 59.52120 80.32820 2.38680
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1    N N   . MET A 1 1   ? 10.710 38.460 14.825  1.00 89.21  ? ? ? ? ? ? 1   MET A N   1
ATOM 2    C CA  . MET A 1 1   ? 11.257 39.553 13.961  1.00 89.21  ? ? ? ? ? ? 1   MET A CA  1
ATOM 3    C C   . MET A 1 1   ? 11.385 40.985 14.516  1.00 89.21  ? ? ? ? ? ? 1   MET A C   1
ATOM 4    O O   . MET A 1 1   ? 12.376 41.648 14.218  1.00 89.21  ? ? ? ? ? ? 1   MET A O   1
ATOM 5    C CB  . MET A 1 1   ? 10.514 39.584 12.633  1.00 72.05  ? ? ? ? ? ? 1   MET A CB  1
ATOM 6    C CG  . MET A 1 1   ? 11.115 38.664 11.596  1.00 72.05  ? ? ? ? ? ? 1   MET A CG  1
ATOM 7    S SD  . MET A 1 1   ? 12.048 39.609 10.386  1.00 72.05  ? ? ? ? ? ? 1   MET A SD  1
ATOM 8    C CE  . MET A 1 1   ? 13.456 40.084 11.391  1.00 72.05  ? ? ? ? ? ? 1   MET A CE  1
ATOM 9    N N   . ASP A 1 2   ? 10.381 41.467 15.263  1.00 81.99  ? ? ? ? ? ? 2   ASP A N   1
ATOM 10   C CA  . ASP A 1 2   ? 10.350 42.822 15.886  1.00 81.99  ? ? ? ? ? ? 2   ASP A CA  1
ATOM 11   C C   . ASP A 1 2   ? 10.651 44.060 15.038  1.00 81.99  ? ? ? ? ? ? 2   ASP A C   1
ATOM 12   O O   . ASP A 1 2   ? 11.725 44.645 15.140  1.00 81.99  ? ? ? ? ? ? 2   ASP A O   1
ATOM 13   C CB  . ASP A 1 2   ? 11.208 42.882 17.167  1.00 70.41  ? ? ? ? ? ? 2   ASP A CB  1
ATOM 14   C CG  . ASP A 1 2   ? 11.000 44.178 17.963  1.00 70.41  ? ? ? ? ? ? 2   ASP A CG  1
ATOM 15   O OD1 . ASP A 1 2   ? 10.015 44.907 17.702  1.00 70.41  ? ? ? ? ? ? 2   ASP A OD1 1
ATOM 16   O OD2 . ASP A 1 2   ? 11.821 44.453 18.866  1.00 70.41  ? ? ? ? ? ? 2   ASP A OD2 1
#
"""

pdb_str_3a="""
data_default
_cell.angle_beta                  90.000
_cell.angle_gamma                 120.000
_cell.length_b                    94.730
_cell.length_c                    250.870
_cell.angle_alpha                 90.000
_cell.volume                      1949640.043
_cell.length_a                    94.730
_space_group.crystal_system       hexagonal
_space_group.name_H-M_alt         'P 65'
_space_group.IT_number            170
_space_group.name_Hall            ' P 65'
_symmetry.space_group_name_H-M    'P 65'
_symmetry.Int_Tables_number       170
_symmetry.space_group_name_Hall   ' P 65'
loop_
  _space_group_symop.id
  _space_group_symop.operation_xyz
   1 x,y,z
   2 x-y,x,z+5/6
   3 y,-x+y,z+1/6
   4 -y,x-y,z+2/3
   5 -x+y,-x,z+1/3
   6 -x,-y,z+1/2

loop_
  _atom_site.group_PDB
  _atom_site.id
  _atom_site.label_atom_id
  _atom_site.label_alt_id
  _atom_site.label_comp_id
  _atom_site.auth_asym_id
  _atom_site.auth_seq_id
  _atom_site.pdbx_PDB_ins_code
  _atom_site.Cartn_x
  _atom_site.Cartn_y
  _atom_site.Cartn_z
  _atom_site.occupancy
  _atom_site.B_iso_or_equiv
  _atom_site.type_symbol
  _atom_site.pdbx_formal_charge
  _atom_site.label_asym_id
  _atom_site.label_entity_id
  _atom_site.label_seq_id
  _atom_site.pdbx_PDB_model_num
   ATOM 1 N . MET A1 1 ? 10.71000 38.46000 14.82500 1.000 89.21000 N ? A ? 1 1
   ATOM 2 CA . MET A1 1 ? 11.25700 39.55300 13.96100 1.000 89.21000 C ? A ? 1 1
   ATOM 3 C . MET A1 1 ? 11.38500 40.98500 14.51600 1.000 89.21000 C ? A ? 1 1
   ATOM 4 O . MET A1 1 ? 12.37600 41.64800 14.21800 1.000 89.21000 O ? A ? 1 1
   ATOM 5 CB . MET A1 1 ? 10.51400 39.58400 12.63300 1.000 72.05000 C ? A ? 1 1
   ATOM 6 CG . MET A1 1 ? 11.11500 38.66400 11.59600 1.000 72.05000 C ? A ? 1 1
   ATOM 7 SD . MET A1 1 ? 12.04800 39.60900 10.38600 1.000 72.05000 S ? A ? 1 1
   ATOM 8 CE . MET A1 1 ? 13.45600 40.08400 11.39100 1.000 72.05000 C ? A ? 1 1
   ATOM 9 N . ASP A1 2 ? 10.38100 41.46700 15.26300 1.000 81.99000 N ? A ? 2 1
   ATOM 10 CA . ASP A1 2 ? 10.35000 42.82200 15.88600 1.000 81.99000 C ? A ? 2 1
   ATOM 11 C . ASP A1 2 ? 10.65100 44.06000 15.03800 1.000 81.99000 C ? A ? 2 1
   ATOM 12 O . ASP A1 2 ? 11.72500 44.64500 15.14000 1.000 81.99000 O ? A ? 2 1
   ATOM 13 CB . ASP A1 2 ? 11.20800 42.88200 17.16700 1.000 70.41000 C ? A ? 2 1
   ATOM 14 CG . ASP A1 2 ? 11.00000 44.17800 17.96300 1.000 70.41000 C ? A ? 2 1
   ATOM 15 OD1 . ASP A1 2 ? 10.01500 44.90700 17.70200 1.000 70.41000 O ? A ? 2 1
   ATOM 16 OD2 . ASP A1 2 ? 11.82100 44.45300 18.86600 1.000 70.41000 O ? A ? 2 1
   ATOM 1 N . MET A2 1 ? 47.80830 41.74364 16.08704 1.000 89.21000 N ? B ? 3 1
   ATOM 2 CA . MET A2 1 ? 47.32528 40.65782 15.17706 1.000 89.21000 C ? B ? 3 1
   ATOM 3 C . MET A2 1 ? 47.15903 39.22120 15.70953 1.000 89.21000 C ? B ? 3 1
   ATOM 4 O . MET A2 1 ? 46.19237 38.56056 15.33594 1.000 89.21000 O ? B ? 3 1
   ATOM 5 CB . MET A2 1 ? 48.16131 40.63807 13.90535 1.000 72.05000 C ? B ? 3 1
   ATOM 6 CG . MET A2 1 ? 47.63526 41.56662 12.83576 1.000 72.05000 C ? B ? 3 1
   ATOM 7 SD . MET A2 1 ? 46.79184 40.63164 11.55433 1.000 72.05000 S ? B ? 3 1
   ATOM 8 CE . MET A2 1 ? 45.31598 40.14801 12.45213 1.000 72.05000 C ? B ? 3 1
   ATOM 9 N . ASP A2 2 ? 48.10744 38.73313 16.52233 1.000 81.99000 N ? B ? 4 1
   ATOM 10 CA . ASP A2 2 ? 48.09487 37.37296 17.13462 1.000 81.99000 C ? B ? 4 1
   ATOM 11 C . ASP A2 2 ? 47.85618 36.14206 16.25698 1.000 81.99000 C ? B ? 4 1
   ATOM 12 O . ASP A2 2 ? 46.77809 35.55605 16.27707 1.000 81.99000 O ? B ? 4 1
   ATOM 13 CB . ASP A2 2 ? 47.14756 37.30209 18.35048 1.000 70.41000 C ? B ? 4 1
   ATOM 14 CG . ASP A2 2 ? 47.29912 35.99951 19.14847 1.000 70.41000 C ? B ? 4 1
   ATOM 15 OD1 . ASP A2 2 ? 48.30081 35.27288 18.95247 1.000 70.41000 O ? B ? 4 1
   ATOM 16 OD2 . ASP A2 2 ? 46.41590 35.71682 19.98816 1.000 70.41000 O ? B ? 4 1
"""
# Test a virus structure and cif biological assemblies
pdb_str_4 = """
data_1qgt

_cell.entry_id           1QGT
_cell.length_a           538.400
_cell.length_b           354.800
_cell.length_c           370.100
_cell.angle_alpha        90.00
_cell.angle_beta         132.30
_cell.angle_gamma        90.00
_cell.Z_PDB              480
_cell.pdbx_unique_axis   ?
#
_symmetry.entry_id                         1QGT
_symmetry.space_group_name_H-M             'C 1 2 1'
_symmetry.pdbx_full_space_group_name_H-M   ?
_symmetry.cell_setting                     ?
_symmetry.Int_Tables_number                5

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1    N N   . MET A 1 1   ? 94.766  87.569  60.599  1.00 82.49  ? 1   MET C N   1
ATOM 2    C CA  . MET A 1 1   ? 95.063  86.713  59.420  1.00 82.49  ? 1   MET C CA  1
ATOM 3    C C   . MET A 1 1   ? 96.298  85.843  59.619  1.00 82.49  ? 1   MET C C   1
ATOM 4    O O   . MET A 1 1   ? 97.405  86.233  59.238  1.00 82.49  ? 1   MET C O   1
ATOM 5    C CB  . MET A 1 1   ? 95.282  87.593  58.195  1.00 56.52  ? 1   MET C CB  1
ATOM 6    C CG  . MET A 1 1   ? 94.060  87.760  57.324  1.00 56.52  ? 1   MET C CG  1
ATOM 7    S SD  . MET A 1 1   ? 94.423  87.384  55.596  1.00 56.52  ? 1   MET C SD  1
ATOM 1137 N N   . MET B 1 1   ? 114.364 79.329  34.867  1.00 68.58  ? 1   MET D N   1
ATOM 1138 C CA  . MET B 1 1   ? 112.979 79.421  34.349  1.00 68.58  ? 1   MET D CA  1
ATOM 1139 C C   . MET B 1 1   ? 112.397 80.800  34.542  1.00 68.58  ? 1   MET D C   1
ATOM 1140 O O   . MET B 1 1   ? 112.555 81.420  35.588  1.00 68.58  ? 1   MET D O   1
ATOM 1141 C CB  . MET B 1 1   ? 112.091 78.418  35.056  1.00 58.39  ? 1   MET D CB  1
ATOM 1142 C CG  . MET B 1 1   ? 112.508 78.119  36.469  1.00 58.39  ? 1   MET D CG  1
ATOM 1143 S SD  . MET B 1 1   ? 111.735 76.594  36.985  1.00 58.39  ? 1   MET D SD  1
ATOM 1144 C CE  . MET B 1 1   ? 111.016 77.085  38.539  1.00 58.39  ? 1   MET D CE  1
ATOM 2281 N N   . MET C 1 1   ? 95.532  53.692  97.007  1.00 95.02  ? 1   MET B N   1
ATOM 2282 C CA  . MET C 1 1   ? 95.276  52.408  96.310  1.00 95.02  ? 1   MET B CA  1
ATOM 2283 C C   . MET C 1 1   ? 96.481  52.045  95.462  1.00 95.02  ? 1   MET B C   1
ATOM 2284 O O   . MET C 1 1   ? 97.428  52.822  95.372  1.00 95.02  ? 1   MET B O   1
ATOM 2285 C CB  . MET C 1 1   ? 94.045  52.564  95.433  1.00 60.11  ? 1   MET B CB  1
ATOM 2286 C CG  . MET C 1 1   ? 92.963  53.443  96.068  1.00 60.11  ? 1   MET B CG  1
ATOM 2287 S SD  . MET C 1 1   ? 92.386  54.911  95.149  1.00 60.11  ? 1   MET B SD  1
ATOM 2288 C CE  . MET C 1 1   ? 93.715  55.160  94.042  1.00 60.11  ? 1   MET B CE  1
ATOM 3425 N N   . MET D 1 1   ? 102.469 72.436  72.485  1.00 69.18  ? 1   MET A N   1
ATOM 3426 C CA  . MET D 1 1   ? 100.988 72.350  72.525  1.00 69.18  ? 1   MET A CA  1
ATOM 3427 C C   . MET D 1 1   ? 100.463 72.810  73.873  1.00 69.18  ? 1   MET A C   1
ATOM 3428 O O   . MET D 1 1   ? 101.137 72.689  74.892  1.00 69.18  ? 1   MET A O   1
ATOM 3429 C CB  . MET D 1 1   ? 100.535 70.914  72.256  1.00 55.96  ? 1   MET A CB  1
ATOM 3430 C CG  . MET D 1 1   ? 101.249 69.877  73.096  1.00 55.96  ? 1   MET A CG  1
ATOM 3431 S SD  . MET D 1 1   ? 101.376 68.290  72.255  1.00 55.96  ? 1   MET A SD  1
ATOM 3432 C CE  . MET D 1 1   ? 101.716 67.165  73.624  1.00 55.96  ? 1   MET A CE  1
ATOM 3433 N N   . ASP D 1 2   ? 99.252  73.349  73.862  1.00 56.63  ? 2   ASP A N   1
ATOM 3434 C CA  . ASP D 1 2   ? 98.608  73.838  75.074  1.00 56.63  ? 2   ASP A CA  1


loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 'complete icosahedral assembly'                ? 'complete icosahedral assembly' 240
2 'icosahedral asymmetric unit'                  ? tetrameric                      4
3 'icosahedral pentamer'                         ? eicosameric                     20
4 'icosahedral 23 hexamer'                       ? 24-meric                        24
5 'icosahedral asymmetric unit, std point frame' ? tetrameric                      4
6 'crystal asymmetric unit, crystal frame'       ? 120-meric                       120
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 '(1-60)'                       A,B,C,D
2 1                              A,B,C,D
3 '(1-5)'                        A,B,C,D
4 '(1,2,6,10,23,24)'             A,B,C,D
5 P                              A,B,C,D
6 '(X0)(1-10,16-25,36-40,46-50)' A,B,C,D

_pdbx_point_symmetry.entry_id             1QGT
_pdbx_point_symmetry.Schoenflies_symbol   I
_pdbx_point_symmetry.H-M_notation         532


loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.name
_pdbx_struct_oper_list.symmetry_operation
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
P  'transform to point frame'   ? ? -0.01520725 0.80901697  -0.58758853 0.00000 -0.50775026 0.49999999  0.70156231  0.00000
0.86137008  0.30901706  0.40317495  0.00000
X0 'transform to crystal frame' ? ? 1.00000000  0.00000000  0.00000000  0.00000 0.00000000  1.00000000  0.00000000  0.00000
0.00000000  0.00000000  1.00000000  0.00000
1  'point symmetry operation'   ? ? 1.00000000  0.00000000  0.00000000  0.00000 0.00000000  1.00000000  0.00000000  0.00000
0.00000000  0.00000000  1.00000000  0.00000
2  'point symmetry operation'   ? ? 0.45892970  -0.50775021 0.72909071  0.00000 0.84616288  0.50000004  -0.18441348 0.00000
-0.27090940 0.70156232  0.65910425  0.00000
3  'point symmetry operation'   ? ? -0.41654043 0.02460579  0.90878415  0.00000 0.86137010  -0.30901689 0.40317504  0.00000
0.29075009  0.95073820  0.10752333  0.00000
4  'point symmetry operation'   ? ? -0.41654043 0.86137010  0.29075009  0.00000 0.02460579  -0.30901689 0.95073820  0.00000
0.90878415  0.40317504  0.10752333  0.00000
5  'point symmetry operation'   ? ? 0.45892970  0.84616288  -0.27090940 0.00000 -0.50775021 0.50000004  0.70156232  0.00000
0.72909071  -0.18441348 0.65910425  0.00000
6  'point symmetry operation'   ? ? 0.48391682  0.53235611  0.69456568  0.00000 0.53235610  -0.80901691 0.24917588  0.00000
0.69456568  0.24917588  -0.67489991 0.00000
7  'point symmetry operation'   ? ? 0.48437941  0.50775032  0.71243681  0.00000 -0.50775014 -0.50000000 0.70156240  0.00000
0.71243693  -0.70156227 0.01562059  0.00000
8  'point symmetry operation'   ? ? 0.45892974  0.50775025  0.72909065  0.00000 -0.84616291 0.49999996  0.18441359  0.00000
-0.27090925 -0.70156234 0.65910429  0.00000
9  'point symmetry operation'   ? ? 0.44273840  0.53235600  0.72151216  0.00000 -0.01520726 0.80901701  -0.58758847 0.00000
-0.89652186 0.24917576  0.36627857  0.00000
10 'point symmetry operation'   ? ? 0.45818126  0.54756325  0.70017456  0.00000 0.83676435  0.00000010  -0.54756317 0.00000
-0.29982554 0.83676429  -0.45818136 0.00000
11 'point symmetry operation'   ? ? -0.48437934 -0.50775025 -0.71243690 0.00000 -0.50775025 -0.50000001 0.70156230  0.00000
-0.71243690 0.70156230  -0.01562064 0.00000
12 'point symmetry operation'   ? ? -0.45892964 -0.50775032 -0.72909067 0.00000 -0.84616295 0.49999995  0.18441345  0.00000
0.27090931  0.70156231  -0.65910430 0.00000
13 'point symmetry operation'   ? ? -0.44273840 -0.53235610 -0.72151209 0.00000 -0.01520725 0.80901693  -0.58758859 0.00000
0.89652186  -0.24917582 -0.36627853 0.00000
14 'point symmetry operation'   ? ? -0.45818137 -0.54756325 -0.70017449 0.00000 0.83676431  -0.00000003 -0.54756323 0.00000
0.29982548  -0.83676429 0.45818140  0.00000
15 'point symmetry operation'   ? ? -0.48391689 -0.53235600 -0.69456571 0.00000 0.53235599  -0.80901699 0.24917584  0.00000
-0.69456572 -0.24917582 0.67489990  0.00000
16 'point symmetry operation'   ? ? -0.99953748 -0.02460585 0.01787121  0.00000 -0.02460585 0.30901692  -0.95073818 0.00000
0.01787121  -0.95073818 -0.30947944 0.00000
17 'point symmetry operation'   ? ? -0.48437947 0.50775020  -0.71243684 0.00000 0.50775020  -0.49999999 -0.70156236 0.00000
-0.71243684 -0.70156236 -0.01562054 0.00000
18 'point symmetry operation'   ? ? 0.40034909  0.00000006  -0.91636271 0.00000 0.00000006  -1.00000000 -0.00000004 0.00000
-0.91636271 -0.00000004 -0.40034909 0.00000
19 'point symmetry operation'   ? ? 0.43198340  -0.84616285 -0.31208776 0.00000 -0.84616285 -0.50000009 0.18441350  0.00000
-0.31208776 0.18441350  -0.93198331 0.00000
20 'point symmetry operation'   ? ? -0.43319407 -0.86137014 0.26530055  0.00000 -0.86137014 0.30901686  -0.40317499 0.00000
0.26530055  -0.40317499 -0.87582278 0.00000
21 'point symmetry operation'   ? ? -0.44273847 0.01520715  0.89652183  0.00000 0.53235589  0.80901706  0.24917584  0.00000
-0.72151220 0.58758841  -0.36627859 0.00000
22 'point symmetry operation'   ? ? -0.43319429 0.86137006  0.26530044  0.00000 0.86137006  0.30901713  0.40317494  0.00000
0.26530044  0.40317494  -0.87582284 0.00000
23 'point symmetry operation'   ? ? 0.45818126  0.83676435  -0.29982554 0.00000 0.54756325  0.00000010  0.83676429  0.00000
0.70017456  -0.54756317 -0.45818136 0.00000
24 'point symmetry operation'   ? ? 0.99953748  -0.02460572 -0.01787121 0.00000 0.02460581  0.30901699  0.95073816  0.00000
-0.01787109 -0.95073816 0.30947951  0.00000
25 'point symmetry operation'   ? ? 0.44273847  -0.53235599 0.72151212  0.00000 0.01520714  0.80901697  0.58758853  0.00000
-0.89652183 -0.24917590 0.36627855  0.00000
26 'point symmetry operation'   ? ? 0.41654032  -0.02460576 -0.90878420 0.00000 0.86137013  -0.30901696 0.40317491  0.00000
-0.29075015 -0.95073817 -0.10752336 0.00000
27 'point symmetry operation'   ? ? 0.41654043  -0.86137006 -0.29075021 0.00000 0.02460578  -0.30901703 0.95073815  0.00000
-0.90878415 -0.40317501 -0.10752340 0.00000
28 'point symmetry operation'   ? ? -0.45892964 -0.84616295 0.27090931  0.00000 -0.50775032 0.49999995  0.70156231  0.00000
-0.72909067 0.18441345  -0.65910430 0.00000
29 'point symmetry operation'   ? ? -1.00000000 -0.00000013 0.00000000  0.00000 -0.00000013 1.00000000  -0.00000007 0.00000
0.00000000  -0.00000007 -1.00000000 0.00000
30 'point symmetry operation'   ? ? -0.45892981 0.50775014  -0.72909068 0.00000 0.84616285  0.50000005  -0.18441362 0.00000
0.27090934  -0.70156236 -0.65910423 0.00000
31 'point symmetry operation'   ? ? -0.43198330 0.84616291  0.31208774  0.00000 -0.84616288 -0.50000000 0.18441361  0.00000
0.31208782  -0.18441346 0.93198330  0.00000
32 'point symmetry operation'   ? ? 0.43319418  0.86137010  -0.26530050 0.00000 -0.86137010 0.30901699  -0.40317496 0.00000
-0.26530049 0.40317497  0.87582281  0.00000
33 'point symmetry operation'   ? ? 0.99953748  0.02460581  -0.01787109 0.00000 -0.02460572 0.30901699  -0.95073816 0.00000
-0.01787121 0.95073816  0.30947951  0.00000
34 'point symmetry operation'   ? ? 0.48437941  -0.50775014 0.71243693  0.00000 0.50775032  -0.50000000 -0.70156227 0.00000
0.71243681  0.70156240  0.01562059  0.00000
35 'point symmetry operation'   ? ? -0.40034909 0.00000007  0.91636271  0.00000 0.00000007  -1.00000000 0.00000011  0.00000
0.91636271  0.00000011  0.40034909  0.00000
36 'point symmetry operation'   ? ? 0.45818144  -0.83676431 -0.29982537 0.00000 -0.54756314 -0.00000010 -0.83676437 0.00000
0.70017453  0.54756323  -0.45818134 0.00000
37 'point symmetry operation'   ? ? -0.41654032 -0.86137009 0.29075027  0.00000 -0.02460574 -0.30901710 -0.95073813 0.00000
0.90878420  -0.40317489 0.10752343  0.00000
38 'point symmetry operation'   ? ? -0.99878911 -0.01520721 0.04678732  0.00000 -0.01520721 -0.80901704 -0.58758844 0.00000
0.04678732  -0.58758844 0.80780615  0.00000
39 'point symmetry operation'   ? ? -0.48391689 0.53235599  -0.69456572 0.00000 -0.53235600 -0.80901699 -0.24917582 0.00000
-0.69456571 0.24917584  0.67489990  0.00000
40 'point symmetry operation'   ? ? 0.41654043  0.02460578  -0.90878415 0.00000 -0.86137006 -0.30901703 -0.40317501 0.00000
-0.29075021 0.95073815  -0.10752340 0.00000
41 'point symmetry operation'   ? ? -0.44273847 0.53235589  -0.72151220 0.00000 0.01520715  0.80901706  0.58758841  0.00000
0.89652183  0.24917584  -0.36627859 0.00000
42 'point symmetry operation'   ? ? 0.44273840  -0.01520726 -0.89652186 0.00000 0.53235600  0.80901701  0.24917576  0.00000
0.72151216  -0.58758847 0.36627857  0.00000
43 'point symmetry operation'   ? ? 0.43319418  -0.86137010 -0.26530049 0.00000 0.86137010  0.30901699  0.40317497  0.00000
-0.26530050 -0.40317496 0.87582281  0.00000
44 'point symmetry operation'   ? ? -0.45818133 -0.83676435 0.29982543  0.00000 0.54756314  0.00000003  0.83676437  0.00000
-0.70017460 0.54756317  0.45818130  0.00000
45 'point symmetry operation'   ? ? -0.99953749 0.02460569  0.01787109  0.00000 0.02460569  0.30901707  0.95073814  0.00000
0.01787109  0.95073814  -0.30947958 0.00000
46 'point symmetry operation'   ? ? -0.43198330 -0.84616288 0.31208782  0.00000 0.84616291  -0.50000000 -0.18441346 0.00000
0.31208774  0.18441361  0.93198330  0.00000
47 'point symmetry operation'   ? ? -0.99878911 0.01520719  0.04678725  0.00000 0.01520719  -0.80901695 0.58758856  0.00000
0.04678725  0.58758856  0.80780606  0.00000
48 'point symmetry operation'   ? ? -0.45818133 0.54756314  -0.70017460 0.00000 -0.83676435 0.00000003  0.54756317  0.00000
0.29982543  0.83676437  0.45818130  0.00000
49 'point symmetry operation'   ? ? 0.44273847  0.01520714  -0.89652183 0.00000 -0.53235599 0.80901697  -0.24917590 0.00000
0.72151212  0.58758853  0.36627855  0.00000
50 'point symmetry operation'   ? ? 0.45892974  -0.84616291 -0.27090925 0.00000 0.50775025  0.49999996  -0.70156234 0.00000
0.72909065  0.18441359  0.65910429  0.00000
51 'point symmetry operation'   ? ? 0.45818144  -0.54756314 0.70017453  0.00000 -0.83676431 -0.00000010 0.54756323  0.00000
-0.29982537 -0.83676437 -0.45818134 0.00000
52 'point symmetry operation'   ? ? -0.44273840 -0.01520725 0.89652186  0.00000 -0.53235610 0.80901693  -0.24917582 0.00000
-0.72151209 -0.58758859 -0.36627853 0.00000
53 'point symmetry operation'   ? ? -0.45892981 0.84616285  0.27090934  0.00000 0.50775014  0.50000005  -0.70156236 0.00000
-0.72909068 -0.18441362 -0.65910423 0.00000
54 'point symmetry operation'   ? ? 0.43198319  0.84616294  -0.31208780 0.00000 0.84616294  -0.49999991 -0.18441357 0.00000
-0.31208780 -0.18441357 -0.93198328 0.00000
55 'point symmetry operation'   ? ? 0.99878911  -0.01520709 -0.04678732 0.00000 0.01520731  -0.80901699 0.58758850  0.00000
-0.04678725 -0.58758850 -0.80780611 0.00000
56 'point symmetry operation'   ? ? 0.41654032  0.86137013  -0.29075015 0.00000 -0.02460576 -0.30901696 -0.95073817 0.00000
-0.90878420 0.40317491  -0.10752336 0.00000
57 'point symmetry operation'   ? ? 0.99878911  0.01520731  -0.04678725 0.00000 -0.01520709 -0.80901699 -0.58758850 0.00000
-0.04678732 0.58758850  -0.80780611 0.00000
58 'point symmetry operation'   ? ? 0.48391696  -0.53235589 0.69456575  0.00000 -0.53235589 -0.80901708 -0.24917578 0.00000
0.69456575  -0.24917578 -0.67489988 0.00000
59 'point symmetry operation'   ? ? -0.41654032 -0.02460574 0.90878420  0.00000 -0.86137009 -0.30901710 -0.40317489 0.00000
0.29075027  -0.95073813 0.10752343  0.00000
60 'point symmetry operation'   ? ? -0.45818137 0.83676431  0.29982548  0.00000 -0.54756325 -0.00000003 -0.83676429 0.00000
-0.70017449 -0.54756323 0.45818140  0.00000

"""


def exercise_000(file_name="tst_mtrix_biomt_cmdl_000.pdb"):
  """
  Make sure SS gets populated by BIOMT
  """

  of = open(file_name,"w")
  print(pdb_str_0, file=of)
  of.close()
  # template when pdb_as_cif will handle BIOMT records as well...
  # easy_run.call("phenix.pdb_as_cif %s"%file_name)
  # for file_type in ['pdb', 'cif']:
  for file_type in ['pdb']:
    print("file_type:", file_type)
    easy_run.call("phenix.pdb.biomt_reconstruction %s.%s" % (file_name[:-4], file_type))
    pdb_inp = iotbx.pdb.input(
      file_name="%s_BIOMT_expanded.%s" % (file_name[:-4], file_type))
    a = pdb_inp.extract_secondary_structure()
    assert a.get_n_helices() == 27
    assert a.get_n_sheets() == 9, "%d" % a.get_n_sheets()
    # checking chain ids. If this part is failing, then something is changed in
    # chain expanding which made chain ids in hierarchy.py:join_roots()
    # not compatible with those used in secondary_structure.py:multiply_to_asu
    chain_ids = [h.start_chain_id for h in a.helices]
    assert chain_ids == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K',
        'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V',
        'W', 'X', 'Y', 'Z', '0'], chain_ids
    # checking sheets
    for i, sh in enumerate(a.sheets):
      assert sh.n_strands == 2
      assert sh.registrations[0] == None
      assert sh.registrations[1].cur_chain_id == chain_ids[i*3+1]
      assert sh.registrations[1].prev_chain_id == chain_ids[i*3]
      assert sh.strands[0].start_chain_id == chain_ids[i*3], sh.strands[0].start_chain_id
      assert sh.strands[0].end_chain_id == chain_ids[i*3]
      assert sh.strands[1].start_chain_id == chain_ids[i*3+1], sh.strands[1].start_chain_id
      assert sh.strands[1].end_chain_id == chain_ids[i*3+1]

def exercise_001(file_name="tst_mtrix_biomt_cmdl_001.pdb"):
  """
  Make sure SS gets populated by MTRIX
  """
  of = open(file_name,"w")
  print(pdb_str_2b, file=of)
  of.close()
  easy_run.call("phenix.pdb.mtrix_reconstruction %s"%file_name)
  pdb_inp = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_001_MTRIX_expanded.pdb")
  a = pdb_inp.extract_secondary_structure()
  assert a.get_n_helices() == 3, a.get_n_helices()
  assert a.get_n_sheets() == 3, "%d" % a.get_n_sheets()
  # checking chain ids. If this part is failing, then something is changed in
  # chain expanding which made chain ids in hierarchy.py:join_roots()
  # not compatible with those used in secondary_structure.py:multiply_to_asu
  chain_ids = [h.start_chain_id for h in a.helices]
  assert chain_ids == ['A', 'C', 'E'], chain_ids
  # checking sheets
  assert not show_diff (a.as_pdb_str(),"""\
HELIX    1   1 THR A    1  THR A    2  1                                   6
HELIX    2   2 THR C    1  THR C    2  1                                   6
HELIX    3   3 THR E    1  THR E    2  1                                   6
SHEET    1   1 2 THR A   1  THR A   3  0
SHEET    2   1 2 THR B   4  THR B   5 -1  O  THR B   4   N  THR A   2
SHEET    1   2 2 THR C   1  THR C   3  0
SHEET    2   2 2 THR D   4  THR D   5 -1  O  THR D   4   N  THR C   2
SHEET    1   3 2 THR E   1  THR E   3  0
SHEET    2   3 2 THR F   4  THR F   5 -1  O  THR F   4   N  THR E   2""")

def exercise_002(file_name="tst_mtrix_biomt_cmdl_002.pdb"):
  """
  Make sure SS is fine when expanding without MTRIX
  """
  of = open(file_name,"w")
  print(pdb_str_0, file=of)
  of.close()
  easy_run.call("phenix.pdb.mtrix_reconstruction %s"%file_name)
  pdb_inp = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_002_MTRIX_expanded.pdb")
  a = pdb_inp.extract_secondary_structure()
  assert a.get_n_helices() == 3, a.get_n_helices()
  assert a.get_n_sheets() == 1, "%d" % a.get_n_sheets()
  # checking chain ids. If this part is failing, then something is changed in
  # chain expanding which made chain ids in hierarchy.py:join_roots()
  # not compatible with those used in secondary_structure.py:multiply_to_asu
  chain_ids = [h.start_chain_id for h in a.helices]
  assert chain_ids == ['A', 'B', 'C'], chain_ids
  # checking sheets
  for i, sh in enumerate(a.sheets):
    assert sh.n_strands == 2
    assert sh.registrations[0] == None
    assert sh.registrations[1].cur_chain_id == 'B'
    assert sh.registrations[1].prev_chain_id == 'A'
    assert sh.strands[0].start_chain_id == 'A'
    assert sh.strands[0].end_chain_id == 'A'
    assert sh.strands[1].start_chain_id == 'B'
    assert sh.strands[1].end_chain_id == 'B'

def exercise_003(file_name="tst_mtrix_biomt_cmdl_003.pdb"):
  """
  Make sure SS is fine when expanding without BIOMT
  """
  of = open(file_name,"w")
  print(pdb_str_2b, file=of)
  of.close()
  easy_run.call("phenix.pdb.biomt_reconstruction %s"%file_name)
  pdb_inp = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_003_BIOMT_expanded.pdb")
  a = pdb_inp.extract_secondary_structure()
  assert a.get_n_helices() == 1, a.get_n_helices()
  assert a.get_n_sheets() == 1, "%d" % a.get_n_sheets()
  # checking chain ids. If this part is failing, then something is changed in
  # chain expanding which made chain ids in hierarchy.py:join_roots()
  # not compatible with those used in secondary_structure.py:multiply_to_asu
  chain_ids = [h.start_chain_id for h in a.helices]
  assert chain_ids == ['A'], chain_ids
  # checking sheets
  for i, sh in enumerate(a.sheets):
    assert sh.n_strands == 2
    assert sh.registrations[0] == None
    assert sh.registrations[1].cur_chain_id == 'B'
    assert sh.registrations[1].prev_chain_id == 'A'
    assert sh.strands[0].start_chain_id == 'A'
    assert sh.strands[0].end_chain_id == 'A'
    assert sh.strands[1].start_chain_id == 'B'
    assert sh.strands[1].end_chain_id == 'B'

def exercise_01(file_name="tst_mtrix_biomt_cmdl_01.pdb"):
  """
  Use MTRIX.
  """
  of = open(file_name,"w")
  print(pdb_str_1, file=of)
  of.close()
  easy_run.call("phenix.pdb.mtrix_reconstruction %s"%file_name)
  pdb_inp1 = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_01_MTRIX_expanded.pdb")
  pdb_inp2 = iotbx.pdb.input(source_info=None, lines=pdb_str_1a)
  dist = flex.sqrt((
    pdb_inp1.atoms().extract_xyz() - pdb_inp2.atoms().extract_xyz()).dot())
  assert approx_equal(dist.min_max_mean().as_tuple(), [0,0,0])
  assert pdb_inp1.crystal_symmetry().is_similar_symmetry(
         pdb_inp2.crystal_symmetry())
  assert pdb_inp1.crystal_symmetry() is not None
  assert pdb_inp2.crystal_symmetry() is not None
  assert pdb_inp1.crystal_symmetry().unit_cell() is not None
  assert pdb_inp2.crystal_symmetry().unit_cell() is not None

def exercise_02(file_name="tst_mtrix_biomt_cmdl_02.pdb"):
  """
  Make sure it's not going to make up CRYST1.
  """
  of = open(file_name,"w")
  print(pdb_str_2, file=of)
  of.close()
  easy_run.call("phenix.pdb.mtrix_reconstruction %s"%file_name)
  pdb_inp1 = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_02_MTRIX_expanded.pdb")
  pdb_inp2 = iotbx.pdb.input(source_info=None, lines=pdb_str_2a)
  dist = flex.sqrt((
    pdb_inp1.atoms().extract_xyz() - pdb_inp2.atoms().extract_xyz()).dot())
  assert approx_equal(dist.min_max_mean().as_tuple(), [0,0,0])
  assert [pdb_inp1.crystal_symmetry(),
          pdb_inp2.crystal_symmetry()].count(None)==2

def exercise_03(file_name="tst_mtrix_biomt_cmdl_03.cif"):
  """
  Handle mmCIF.
  """
  of = open(file_name,"w")
  print(pdb_str_3, file=of)
  of.close()
  easy_run.call("phenix.pdb.mtrix_reconstruction %s"%file_name)
  pdb_inp1 = iotbx.pdb.input(
    file_name="tst_mtrix_biomt_cmdl_03_MTRIX_expanded.cif")
  pdb_inp2 = iotbx.pdb.input(source_info=None, lines=pdb_str_3a)
  dist = flex.sqrt((
    pdb_inp1.atoms().extract_xyz() - pdb_inp2.atoms().extract_xyz()).dot())
  assert approx_equal(dist.min_max_mean().as_tuple(), [0,0,0])
  assert pdb_inp1.crystal_symmetry().is_similar_symmetry(
         pdb_inp2.crystal_symmetry())
  assert pdb_inp1.crystal_symmetry() is not None
  assert pdb_inp2.crystal_symmetry() is not None
  assert pdb_inp1.crystal_symmetry().unit_cell() is not None
  assert pdb_inp2.crystal_symmetry().unit_cell() is not None

def exercise_04(file_name="tst_mtrix_biomt_cmdl_04.cif"):
  """
  Test that the number of chains is correct after expansion
  """
  from iotbx.pdb.mmcif import cif_input
  from mmtbx.model.model import manager as model_manager

  with open(file_name,"w") as fh:
    print(pdb_str_4,file=fh)

  # open the expanded version
  easy_run.call("phenix.pdb.biomt_reconstruction %s"%file_name)
  inp_biomt = cif_input(
    file_name="tst_mtrix_biomt_cmdl_04_BIOMT_expanded.cif")

  # open the asu
  inp_asu= cif_input(file_name=file_name)

  # make models
  model_asu = model_manager(inp_asu)
  model_biomt = model_manager(inp_biomt)

  # check asu
  counts = model_asu.get_hierarchy().overall_counts()
  nres = counts.n_residues
  nchains = counts.n_chains
  assert nres == 5
  assert nchains == 4

  # check expanded from command line
  counts = model_biomt.get_hierarchy().overall_counts()
  nres = counts.n_residues
  nchains = counts.n_chains
  # 60 is number of expected point symmetry ops
  assert nres == 5*60
  assert nchains == 4*60

def exercise_05(file_name="tst_mtrix_biomt_cmdl_05.cif"):
  """
  Test that the number of chains is correct choosing different assemblies
  """
  from iotbx.pdb.mmcif import cif_input
  from mmtbx.model.model import manager as model_manager

  with open(file_name,"w") as fh:
    print(pdb_str_4,file=fh)

  # open the asu
  inp_asu= cif_input(file_name=file_name)
  model_asu = model_manager(inp_asu)
  assert model_asu._assembly_operators is not None

  expanded_sizes = [240,4,20,24,4,120] # n chains expected
  for i,(key,operator) in enumerate(model_asu._assembly_operators.items()):
    model_biomt = model_manager(model_asu._model_input,expand_with_mtrix=False)
    model_biomt.expand_with_BIOMT_records(assembly_id=key)
    with open("excercise_05_expand"+str(i)+".cif", "w") as fh:
      print(model_biomt.model_as_mmcif(), file=fh)
    counts = model_biomt.get_hierarchy().overall_counts()
    nchains = counts.n_chains
    #print(nchains,expanded_sizes[i])
    assert nchains == expanded_sizes[i]



if(__name__=='__main__'):
  if libtbx.env.has_module("phenix"):
    # SS records handling
    exercise_000()
    exercise_001()
    exercise_002()
    exercise_003()
    # # other
    exercise_01()
    exercise_02()
    exercise_03()
    exercise_04()
    exercise_05()
  else:
    print("Skipping tests: phenix not found")
