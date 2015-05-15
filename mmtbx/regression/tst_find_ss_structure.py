from __future__ import division

from libtbx.utils import null_out
from cStringIO import StringIO
from mmtbx.secondary_structure.find_ss_from_ca import find_secondary_structure

def remove_blank(text):
  return text.replace(" ","").replace("\n","")

def tst_01():
  text="""
ATOM      1  N   GLY U  11     -20.099   8.864  10.230  1.00 52.86           N
ATOM      2  CA  GLY U  11     -19.400   7.571   9.992  1.00 52.86           C
ATOM      3  C   GLY U  11     -17.867   8.054   9.511  1.00 52.86           C
ATOM      4  O   GLY U  11     -18.033   8.077   8.293  1.00 52.86           O
ATOM      5  N   GLY U  12     -16.673   8.172  10.087  1.00 52.86           N
ATOM      6  CA  GLY U  12     -15.517   7.978   9.123  1.00 52.86           C
ATOM      7  C   GLY U  12     -15.338   6.468   9.675  1.00 52.86           C
ATOM      8  O   GLY U  12     -15.070   5.904  10.738  1.00 52.86           O
ATOM      9  N   GLY U  13     -15.479   5.785   8.520  1.00 52.86           N
ATOM     10  CA  GLY U  13     -14.900   4.879   9.412  1.00 52.86           C
ATOM     11  C   GLY U  13     -14.424   3.514   9.380  1.00 52.86           C
ATOM     12  O   GLY U  13     -13.660   2.944   8.602  1.00 52.86           O
ATOM     13  N   GLY U  14     -14.286   3.425  10.715  1.00 52.86           N
ATOM     14  CA  GLY U  14     -13.338   2.264  10.393  1.00 52.86           C
ATOM     15  C   GLY U  14     -12.091   1.361   9.704  1.00 52.86           C
ATOM     16  O   GLY U  14     -11.343   0.974   8.806  1.00 52.86           O
ATOM     17  N   GLY U  15     -12.859   0.538  10.412  1.00 52.86           N
ATOM     18  CA  GLY U  15     -12.030  -0.556  10.373  1.00 52.86           C
ATOM     19  C   GLY U  15     -12.456  -1.530   9.072  1.00 52.86           C
ATOM     20  O   GLY U  15     -13.665  -1.682   8.859  1.00 52.86           O
ATOM     21  N   GLY U  16     -11.909  -3.009   7.572  1.00 42.55           N
ATOM     22  CA  GLY U  16     -11.031  -3.180   6.040  1.00 42.55           C
ATOM     23  C   GLY U  16     -10.265  -3.837   5.718  1.00 42.55           C
ATOM     24  O   GLY U  16      -9.515  -4.299   6.577  1.00 42.55           O
ATOM     25  N   GLY U  17     -11.247  -3.944   5.165  1.00 52.86           N
ATOM     26  CA  GLY U  17     -10.273  -5.222   4.456  1.00 52.86           C
ATOM     27  C   GLY U  17      -9.735  -4.585   3.332  1.00 52.86           C
ATOM     28  O   GLY U  17     -10.280  -3.765   2.596  1.00 52.86           O
ATOM     29  N   GLY U  18      -8.496  -4.996   3.135  1.00 52.86           N
ATOM     30  CA  GLY U  18      -8.015  -4.822   2.186  1.00 52.86           C
ATOM     31  C   GLY U  18      -7.002  -5.187   1.126  1.00 52.86           C
ATOM     32  O   GLY U  18      -6.107  -5.649   1.834  1.00 52.86           O
ATOM     33  N   GLY U  19      -7.061  -5.398  -0.184  1.00 52.86           N
ATOM     34  CA  GLY U  19      -6.009  -6.212  -0.948  1.00 52.86           C
ATOM     35  C   GLY U  19      -4.828  -5.254  -1.458  1.00 52.86           C
ATOM     36  O   GLY U  19      -5.131  -4.184  -1.991  1.00 52.86           O
ATOM     37  N   GLY U  20      -3.582  -5.721  -1.443  1.00 52.86           N
ATOM     38  CA  GLY U  20      -2.729  -4.680  -2.268  1.00 52.86           C
ATOM     39  C   GLY U  20      -2.086  -4.815  -3.421  1.00 52.86           C
ATOM     40  O   GLY U  20      -1.108  -5.401  -3.871  1.00 52.86           O
ATOM     41  N   GLY U  21      -2.755  -3.911  -4.145  1.00 52.86           N
ATOM     42  CA  GLY U  21      -2.425  -3.693  -5.384  1.00 52.86           C
ATOM     43  C   GLY U  21      -0.789  -3.646  -5.564  1.00 52.86           C
ATOM     44  O   GLY U  21      -0.196  -4.136  -6.539  1.00 52.86           O
"""
  print "Finding beta strands...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(args=['include_single=True'],hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 11  Start: 11 End: 21
Class:  Beta strand  N: 4 Start: 12 End: 15  Rise: 3.12 A Dot: 1.00
Class:  Beta strand  N: 5 Start: 17 End: 21  Rise: 3.35 A Dot: 0.96

PDB RECORDS:
SHEET    1   1 1 GLY U  12  GLY U  15  0
SHEET    1   2 1 GLY U  17  GLY U  21  0

PDB Selections:
" ( chain 'U' and resseq 12:15 and icode ' ')  or  ( chain 'U' and resseq 17:21 and icode ' ') "

"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"


def tst_02():
  text="""
ATOM      1  N   GLY A   1      42.375 -12.180  24.780  1.00 35.31
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57
ATOM      3  C   GLY A   1      43.288 -10.171  23.615  1.00 34.64
ATOM      4  O   GLY A   1      42.111  -9.896  23.277  1.00 35.82
ATOM      5  N   ILE A   2      44.323  -9.391  23.299  1.00 32.23
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55
ATOM      7  C   ILE A   2      43.750  -8.629  21.119  1.00 24.92
ATOM      8  O   ILE A   2      43.068  -7.904  20.409  1.00 23.73
ATOM      9  CB  ILE A   2      45.525  -7.320  22.425  1.00 30.10
ATOM     10  CG1 ILE A   2      45.924  -6.837  23.820  1.00 29.64
ATOM     11  CG2 ILE A   2      45.555  -6.173  21.386  1.00 30.54
ATOM     12  CD1 ILE A   2      44.837  -6.338  24.762  1.00 32.44
ATOM     13  N   GLY A   3      44.161  -9.867  20.749  1.00 22.69
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05
ATOM     15  C   GLY A   3      42.433 -10.405  19.166  1.00 22.08
ATOM     16  O   GLY A   3      41.912 -10.061  18.096  1.00 22.86
ATOM     17  N   ALA A   4      41.862 -10.961  20.191  1.00 21.60
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80
ATOM     19  C   ALA A   4      39.584  -9.950  20.087  1.00 21.67
ATOM     20  O   ALA A   4      38.676  -9.747  19.278  1.00 22.21
ATOM     21  CB  ALA A   4      40.061 -12.080  21.350  1.00 22.97
ATOM     22  N   VAL A   5      39.936  -9.001  20.956  1.00 22.13
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34
ATOM     24  C   VAL A   5      39.536  -6.896  19.795  1.00 18.81
ATOM     25  O   VAL A   5      38.626  -6.314  19.126  1.00 17.22
ATOM     26  CB  VAL A   5      39.843  -6.933  22.338  1.00 20.40
ATOM     27  CG1 VAL A   5      39.237  -5.519  22.413  1.00 23.06
ATOM     28  CG2 VAL A   5      39.745  -7.587  23.653  1.00 21.67
ATOM     29  N   LEU A   6      40.752  -7.021  19.222  1.00 16.05
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59
ATOM     31  C   LEU A   6      40.230  -6.935  16.870  1.00 20.01
ATOM     32  O   LEU A   6      39.649  -6.121  16.029  1.00 21.90
ATOM     33  CB  LEU A   6      42.627  -6.461  17.880  1.00 24.58
ATOM     34  CG  LEU A   6      43.125  -6.023  16.524  1.00 23.91
ATOM     35  CD1 LEU A   6      42.706  -4.584  16.210  1.00 27.44
ATOM     36  CD2 LEU A   6      44.669  -6.152  16.638  1.00 29.31
ATOM     37  N   LYS A   7      39.981  -8.229  16.721  1.00 19.83
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55
ATOM     39  C   LYS A   7      37.648  -8.063  15.784  1.00 19.04
ATOM     40  O   LYS A   7      37.031  -7.839  14.731  1.00 21.18
ATOM     41  CB  LYS A   7      38.854 -10.176  15.616  1.00 27.62
ATOM     42  CG  LYS A   7      40.011 -10.993  15.144  1.00 40.15
ATOM     43  CD  LYS A   7      39.691 -12.487  15.325  1.00 47.84
ATOM     44  CE  LYS A   7      40.599 -13.394  14.493  1.00 53.11
ATOM     45  NZ  LYS A   7      39.966 -14.755  14.319  1.00 55.47
ATOM     46  N   VAL A   8      37.111  -7.988  16.981  1.00 19.69
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52
ATOM     48  C   VAL A   8      35.776  -5.881  16.885  1.00 20.31
ATOM     49  O   VAL A   8      34.775  -5.402  16.257  1.00 20.12
ATOM     50  CB  VAL A   8      35.113  -7.600  18.562  1.00 23.09
ATOM     51  CG1 VAL A   8      34.774  -9.045  18.851  1.00 24.22
ATOM     52  CG2 VAL A   8      35.769  -6.970  19.726  1.00 27.95
ATOM     68  N   THR A  11      35.315  -5.870  13.113  1.00 21.12
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14
ATOM     70  C   THR A  11      32.874  -6.299  13.140  1.00 24.13
ATOM     71  O   THR A  11      31.823  -5.889  12.630  1.00 27.78
ATOM     72  CB  THR A  11      34.462  -7.925  11.935  1.00 27.78
ATOM     73  OG1 THR A  11      34.702  -8.591  13.199  1.00 31.22
ATOM     74  CG2 THR A  11      35.695  -7.959  11.047  1.00 30.31
ATOM     75  N   GLY A  12      32.918  -6.606  14.420  1.00 23.00
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17
ATOM     77  C   GLY A  12      31.264  -5.168  15.584  1.00 24.45
ATOM     78  O   GLY A  12      30.060  -4.914  15.603  1.00 23.79
ATOM     79  N   LEU A  13      32.195  -4.276  15.948  1.00 22.05
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24
ATOM     81  C   LEU A  13      31.251  -2.004  15.324  1.00 21.24
ATOM     82  O   LEU A  13      30.232  -1.308  15.679  1.00 23.02
ATOM     83  CB  LEU A  13      33.146  -2.183  17.008  1.00 25.55
ATOM     84  CG  LEU A  13      32.913  -1.523  18.351  1.00 27.01
ATOM     85  CD1 LEU A  13      33.999  -0.517  18.760  1.00 27.53
ATOM     86  CD2 LEU A  13      31.587  -0.842  18.531  1.00 24.54
ATOM     87  N   PRO A  14      31.689  -1.979  14.106  1.00 17.38
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52
ATOM     89  C   PRO A  14      29.521  -1.670  12.857  1.00 18.98
ATOM     90  O   PRO A  14      28.657  -0.801  12.744  1.00 16.75
ATOM     91  CB  PRO A  14      31.845  -1.614  11.816  1.00 17.80
ATOM     92  CG  PRO A  14      33.118  -2.205  12.303  1.00 18.05
ATOM     93  CD  PRO A  14      33.098  -2.363  13.769  1.00 17.96
ATOM     94  N   ALA A  15      29.207  -2.952  12.868  1.00 16.39
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10
ATOM     96  C   ALA A  15      27.023  -3.016  13.951  1.00 16.98
ATOM     97  O   ALA A  15      25.872  -2.551  13.769  1.00 16.78
ATOM     98  CB  ALA A  15      27.741  -4.906  12.502  1.00 19.58
ATOM     99  N   LEU A  16      27.570  -3.117  15.127  1.00 15.97
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20
ATOM    101  C   LEU A  16      26.614  -1.169  16.344  1.00 20.28
ATOM    102  O   LEU A  16      25.599  -0.734  16.933  1.00 18.32
ATOM    103  CB  LEU A  16      27.811  -3.027  17.542  1.00 19.70
ATOM    104  CG  LEU A  16      27.384  -2.550  18.921  1.00 22.23
ATOM    105  CD1 LEU A  16      26.031  -3.234  19.257  1.00 27.80
ATOM    106  CD2 LEU A  16      28.445  -2.970  19.933  1.00 21.91
ATOM    107  N   ILE A  17      27.514  -0.365  15.791  1.00 20.97
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41
ATOM    109  C   ILE A  17      26.081   1.392  14.758  1.00 18.17
ATOM    110  O   ILE A  17      25.380   2.240  15.282  1.00 16.46
ATOM    111  CB  ILE A  17      28.579   1.847  15.132  1.00 21.10
ATOM    112  CG1 ILE A  17      29.586   1.858  16.352  1.00 25.66
ATOM    113  CG2 ILE A  17      28.268   3.288  14.691  1.00 22.04
ATOM    114  CD1 ILE A  17      30.856   2.696  16.161  1.00 27.00
ATOM    115  N   SER A  18      25.930   0.759  13.657  1.00 16.97
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98
ATOM    117  C   SER A  18      23.499   0.405  13.438  1.00 18.89
ATOM    118  O   SER A  18      22.557   1.165  13.352  1.00 18.37
ATOM    119  CB  SER A  18      25.076   0.039  11.491  1.00 20.39
ATOM    120  OG  SER A  18      23.902   0.046  10.670  1.00 23.87
ATOM    121  N   TRP A  19      23.512  -0.661  14.161  1.00 17.71
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72
ATOM    123  C   TRP A  19      22.083   0.004  16.012  1.00 18.02
ATOM    124  O   TRP A  19      20.820   0.244  16.160  1.00 16.93
ATOM    125  CB  TRP A  19      22.854  -2.410  15.767  1.00 15.59
ATOM    126  CG  TRP A  19      21.803  -2.993  16.678  1.00 17.94
ATOM    127  CD1 TRP A  19      20.917  -3.950  16.210  1.00 18.03
ATOM    128  CD2 TRP A  19      21.448  -2.745  18.041  1.00 16.03
ATOM    129  NE1 TRP A  19      20.060  -4.304  17.222  1.00 21.28
ATOM    130  CE2 TRP A  19      20.357  -3.624  18.372  1.00 20.45
ATOM    131  CE3 TRP A  19      21.879  -1.892  19.048  1.00 14.41
ATOM    132  CZ2 TRP A  19      19.784  -3.690  19.612  1.00 17.15
ATOM    133  CZ3 TRP A  19      21.292  -1.950  20.288  1.00 17.24
ATOM    134  CH2 TRP A  19      20.230  -2.805  20.601  1.00 15.13
ATOM    135  N   ILE A  20      22.930   0.594  16.823  1.00 14.82
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67
ATOM    137  C   ILE A  20      21.917   2.819  17.080  1.00 17.51
ATOM    138  O   ILE A  20      20.942   3.365  17.655  1.00 17.70
ATOM    139  CB  ILE A  20      23.902   2.192  18.499  1.00 16.25
ATOM    140  CG1 ILE A  20      24.481   0.986  19.363  1.00 15.50
ATOM    141  CG2 ILE A  20      23.599   3.421  19.390  1.00 14.54
ATOM    142  CD1 ILE A  20      26.033   1.304  19.637  1.00 18.18
ATOM    143  N   LYS A  21      22.464   3.177  15.957  1.00 16.61
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84
ATOM    145  C   LYS A  21      20.436   3.910  14.752  1.00 21.02
ATOM    146  O   LYS A  21      19.685   4.899  14.971  1.00 22.80
ATOM    147  CB  LYS A  21      22.699   4.646  13.935  1.00 16.73
ATOM    148  CG  LYS A  21      23.944   5.416  14.471  1.00 22.19
ATOM    149  CD  LYS A  21      24.919   5.669  13.300  1.00 25.86
ATOM    150  CE  LYS A  21      26.173   6.287  13.908  1.00 32.91
ATOM    151  NZ  LYS A  21      27.199   6.564  12.863  1.00 39.11
ATOM    152  N   ARG A  22      20.075   2.728  14.351  1.00 19.18
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38
ATOM    154  C   ARG A  22      17.807   2.365  15.177  1.00 22.46
ATOM    155  O   ARG A  22      16.648   2.899  15.096  1.00 23.24
ATOM    156  CB  ARG A  22      18.733   0.860  13.418  1.00 22.73
ATOM    157  CG  ARG A  22      19.309   0.956  12.004  1.00 22.29
ATOM    158  CD  ARG A  22      19.117  -0.300  11.254  1.00 25.78
ATOM    159  NE  ARG A  22      19.382  -1.562  11.991  1.00 28.94
ATOM    160  CZ  ARG A  22      20.624  -2.104  11.889  1.00 33.77
ATOM    161  NH1 ARG A  22      21.664  -1.554  11.252  1.00 31.21
ATOM    162  NH2 ARG A  22      20.870  -3.271  12.498  1.00 33.11
ATOM    163  N   LYS A  23      18.257   1.937  16.323  1.00 20.49
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62
ATOM    165  C   LYS A  23      17.216   3.361  18.100  1.00 25.47
ATOM    166  O   LYS A  23      16.204   3.554  18.811  1.00 24.62
ATOM    167  CB  LYS A  23      18.257   1.128  18.589  1.00 24.16
ATOM    168  CG  LYS A  23      17.979  -0.388  18.463  1.00 31.03
ATOM    169  CD  LYS A  23      16.858  -0.657  19.514  1.00 38.52
ATOM    170  CE  LYS A  23      16.197  -1.986  19.278  1.00 44.05
ATOM    171  NZ  LYS A  23      15.412  -2.493  20.477  1.00 48.30
ATOM    172  N   ARG A  24      18.155   4.268  17.844  1.00 23.99
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11
ATOM    174  C   ARG A  24      16.996   6.355  17.441  1.00 28.34
ATOM    175  O   ARG A  24      16.224   7.166  18.029  1.00 31.82
ATOM    176  CB  ARG A  24      19.446   6.379  18.188  1.00 22.24
ATOM    177  CG  ARG A  24      20.182   6.236  19.542  1.00 20.01
ATOM    178  CD  ARG A  24      21.577   6.742  19.433  1.00 25.00
ATOM    179  NE  ARG A  24      21.715   8.115  18.950  1.00 24.24
ATOM    180  CZ  ARG A  24      21.745   9.136  19.837  1.00 22.72
ATOM    181  NH1 ARG A  24      21.508   8.848  21.086  1.00 24.10
ATOM    182  NH2 ARG A  24      22.134  10.375  19.512  1.00 24.09
ATOM    183  N   GLN A  25      16.889   6.080  16.192  1.00 29.90
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50
ATOM    185  C   GLN A  25      14.463   6.162  15.554  1.00 39.92
ATOM    186  O   GLN A  25      13.435   6.799  15.227  1.00 42.39
ATOM    187  CB  GLN A  25      16.283   6.713  13.874  1.00 40.82
ATOM    188  CG  GLN A  25      17.607   7.505  13.800  1.00 46.91
ATOM    189  CD  GLN A  25      18.413   7.205  12.564  1.00 50.89
ATOM    190  OE1 GLN A  25      19.596   7.559  12.477  1.00 54.23
ATOM    191  NE2 GLN A  25      17.757   6.512  11.624  1.00 51.61
ATOM    192  N   GLN A  26      14.363   4.945  16.064  1.00 43.06
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66
ATOM    194  C   GLN A  26      12.618   5.071  17.837  1.00 48.00
ATOM    195  O   GLN A  26      11.409   5.344  17.597  1.00 50.85
ATOM    196  CB  GLN A  26      13.144   2.874  16.862  1.00 49.14
ATOM    197  CG  GLN A  26      13.166   1.972  15.625  1.00 54.49
ATOM    198  CD  GLN A  26      13.506   0.551  16.116  1.00 58.10
ATOM    199  OE1 GLN A  26      13.083   0.180  17.239  1.00 59.27
ATOM    200  NE2 GLN A  26      14.304  -0.117  15.290  1.00 57.61
"""
  print "Finding helices...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 8  Start: 1 End: 8
Class:  Alpha helix  N: 8 Start: 1 End: 8  Rise: 1.56 A Dot: 0.98

Model 2  N: 16  Start: 11 End: 26
Class:  Alpha helix  N: 16 Start: 11 End: 26  Rise: 1.58 A Dot: 0.98

PDB RECORDS:
HELIX    1   1 GLY A    1  VAL A    8  1                                   8
HELIX    2   2 THR A   11  GLN A   26  1                                  16

PDB Selections:
" ( chain 'A' and resseq 1:8 and icode ' ')  or  ( chain 'A' and resseq 11:26 and icode ' ') "

"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  print "OK"

if __name__=="__main__":
  tst_01()
  tst_02()
