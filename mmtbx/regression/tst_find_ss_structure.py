from __future__ import division

from libtbx.utils import null_out
from cStringIO import StringIO
from mmtbx.secondary_structure.find_ss_from_ca import find_secondary_structure
from libtbx import test_utils

# def remove_blank(text):
#   return text.replace(" ","").replace("\n","")

def tst_01():
  text="""
ATOM      2  CA  GLY U  11     -19.400   7.571   9.992  1.00 52.86           C
ATOM      6  CA  GLY U  12     -15.517   7.978   9.123  1.00 52.86           C
ATOM     10  CA  GLY U  13     -14.900   4.879   9.412  1.00 52.86           C
ATOM     14  CA  GLY U  14     -13.338   2.264  10.393  1.00 52.86           C
ATOM     18  CA  GLY U  15     -12.030  -0.556  10.373  1.00 52.86           C
ATOM     22  CA  GLY U  16     -11.031  -3.180   6.040  1.00 42.55           C
ATOM     26  CA  GLY U  17     -10.273  -5.222   4.456  1.00 52.86           C
ATOM     30  CA  GLY U  18      -8.015  -4.822   2.186  1.00 52.86           C
ATOM     34  CA  GLY U  19      -6.009  -6.212  -0.948  1.00 52.86           C
ATOM     38  CA  GLY U  20      -2.729  -4.680  -2.268  1.00 52.86           C
ATOM     42  CA  GLY U  21      -2.425  -3.693  -5.384  1.00 52.86           C
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
" ( chain 'U' and resid 12  through 15 )  or  ( chain 'U' and resid 17  through 21 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  assert not test_utils.show_diff(found_text, expected_text)
  # if remove_blank(found_text)!=remove_blank(expected_text):
  #   print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
  #   raise AssertionError, "FAILED"
  print "OK"


def tst_02():
  text="""
ATOM      2  CA  GLY A   1      43.603 -11.488  24.325  1.00 35.57
ATOM      6  CA  ILE A   2      44.200  -8.183  22.475  1.00 27.55
ATOM     14  CA  GLY A   3      43.999 -10.264  19.329  1.00 21.05
ATOM     18  CA  ALA A   4      40.378 -11.260  20.106  1.00 21.80
ATOM     23  CA  VAL A   5      39.355  -7.658  21.083  1.00 19.34
ATOM     30  CA  LEU A   6      41.062  -6.432  17.957  1.00 17.59
ATOM     38  CA  LYS A   7      39.079  -8.646  15.636  1.00 22.55
ATOM     47  CA  VAL A   8      35.792  -7.369  17.211  1.00 20.52
ATOM     69  CA  THR A  11      34.132  -6.405  12.343  1.00 24.14
ATOM     76  CA  GLY A  12      31.584  -6.595  15.140  1.00 24.17
ATOM     80  CA  LEU A  13      31.923  -2.919  16.364  1.00 23.24
ATOM     88  CA  PRO A  14      31.026  -1.278  13.030  1.00 17.52
ATOM     95  CA  ALA A  15      27.822  -3.418  12.724  1.00 17.10
ATOM    100  CA  LEU A  16      26.958  -2.649  16.351  1.00 18.20
ATOM    108  CA  ILE A  17      27.343   1.056  15.618  1.00 20.41
ATOM    116  CA  SER A  18      24.825   0.827  12.744  1.00 19.98
ATOM    122  CA  TRP A  19      22.492  -1.085  15.081  1.00 15.72
ATOM    136  CA  ILE A  20      22.628   1.633  17.766  1.00 15.67
ATOM    144  CA  LYS A  21      21.888   4.236  15.157  1.00 19.84
ATOM    153  CA  ARG A  22      18.740   2.273  14.020  1.00 20.38
ATOM    164  CA  LYS A  23      17.500   1.928  17.550  1.00 22.62
ATOM    173  CA  ARG A  24      18.059   5.674  18.276  1.00 27.11
ATOM    184  CA  GLN A  25      15.836   6.730  15.339  1.00 37.50
ATOM    193  CA  GLN A  26      13.132   4.360  16.583  1.00 46.66
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
" ( chain 'A' and resid 1  through 8 )  or  ( chain 'A' and resid 11  through 26 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  assert not test_utils.show_diff(found_text, expected_text)
  # if remove_blank(found_text)!=remove_blank(expected_text):
  #   print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
  #   raise AssertionError, "FAILED"
  print "OK"

def tst_03():
  text="""
ATOM      2  CA  ALA A   1      11.323  32.055  11.635  1.00 40.00           C
ATOM      7  CA  ALA A   2       8.288  29.768  10.916  1.00 40.00           C
ATOM     12  CA  ALA A   3      10.313  27.854   8.231  1.00 40.00           C
ATOM     17  CA  ALA A   4      13.089  27.116  10.822  1.00 40.00           C
ATOM     22  CA  ALA A   5      10.573  25.488  13.298  1.00 40.00           C
ATOM     27  CA  ALA A   6       9.258  23.514  10.260  1.00 40.00           C
ATOM     32  CA  ALA A   7      12.788  22.543   8.962  1.00 40.00           C
ATOM     37  CA  ALA A   8      13.846  21.459  12.515  1.00 40.00           C
ATOM     42  CA  ALA A   9      10.716  19.261  12.994  1.00 40.00           C
ATOM     47  CA  ALA A  10      11.063  17.985   9.357  1.00 40.00           C
ATOM     52  CA  ALA A  11      14.754  17.018   9.967  1.00 40.00           C
ATOM     57  CA  ALA A  12      13.721  15.483  13.371  1.00 40.00           C
ATOM     62  CA  ALA A  13      10.821  13.516  11.708  1.00 40.00           C
ATOM     67  CA  ALA A  14      13.246  12.367   8.939  1.00 40.00           C
ATOM     72  CA  ALA A  15      15.847  11.407  11.629  1.00 40.00           C
ATOM     77  CA  ALA A  16      13.099   9.317  13.370  1.00 40.00           C
ATOM      2  CA  ALA B   2       1.733  -3.620  -2.296  1.00  1.00
ATOM      7  CA  ALA B   3      -1.902  -4.065  -1.341  1.00  1.00
ATOM     12  CA  ALA B   4      -2.941  -0.441  -1.685  1.00  1.00
ATOM     17  CA  ALA B   5      -0.320   0.578  -4.218  1.00  1.00
ATOM     22  CA  ALA B   6       0.221  -2.836  -5.759  1.00  1.00
ATOM     27  CA  ALA B   7      -3.192  -4.271  -4.973  1.00  1.00
ATOM     32  CA  ALA B   8      -5.081  -0.993  -4.849  1.00  1.00
ATOM     37  CA  ALA B   9      -2.802   0.969  -7.148  1.00  1.00
ATOM     42  CA  ALA B  10      -1.460  -1.967  -9.123  1.00  1.00
ATOM     47  CA  ALA B  11      -4.418  -4.277  -8.632  1.00  1.00
ATOM     52  CA  ALA B  12      -7.044  -1.601  -8.116  1.00  1.00
ATOM     57  CA  ALA B  13      -5.323   1.151 -10.064  1.00  1.00
ATOM     62  CA  ALA B  14      -3.322  -1.073 -12.383  1.00  1.00
ATOM     67  CA  ALA B  15      -5.629  -4.072 -12.291  1.00  1.00
ATOM     72  CA  ALA B  16      -8.822  -2.205 -11.488  1.00  1.00
ATOM     77  CA  ALA B  17      -7.833   1.122 -12.996  1.00  1.00
ATOM     82  CA  ALA B  18      -5.368  -0.211 -15.540  1.00  1.00
ATOM     87  CA  ALA B  19      -6.878  -3.661 -15.920  1.00  1.00
ATOM     92  CA  ALA B  20     -10.423  -2.748 -14.958  1.00  1.00
ATOM     97  CA  ALA B  21     -10.280   0.896 -15.972  1.00  1.00
ATOM    102  CA  ALA B  22      -7.582   0.562 -18.606  1.00  1.00
ATOM      2  CA  ALA C   2       1.202  -3.661  -1.646  1.00  1.00           C
ATOM      7  CA  ALA C   3      -1.466  -2.408  -4.020  1.00  1.00           C
ATOM     12  CA  ALA C   4       1.288  -2.503  -6.614  1.00  1.00           C
ATOM     17  CA  ALA C   5       0.312  -6.139  -7.010  1.00  1.00           C
ATOM     22  CA  ALA C   6      -2.284  -4.816  -9.426  1.00  1.00           C
ATOM     27  CA  ALA C   7       0.502  -5.008 -11.981  1.00  1.00           C
ATOM     32  CA  ALA C   8      -0.579  -8.614 -12.375  1.00  1.00           C
ATOM     37  CA  ALA C   9      -3.100  -7.225 -14.833  1.00  1.00           C
ATOM     42  CA  ALA C  10      -0.285  -7.514 -17.347  1.00  1.00           C
ATOM     47  CA  ALA C  11      -1.470 -11.087 -17.740  1.00  1.00           C
ATOM     52  CA  ALA C  12      -3.913  -9.634 -20.239  1.00  1.00           C
ATOM     57  CA  ALA C  13      -1.074 -10.021 -22.713  1.00  1.00           C
ATOM     62  CA  ALA C  14      -2.362 -13.558 -23.106  1.00  1.00           C
ATOM     67  CA  ALA C  15      -4.725 -12.045 -25.646  1.00  1.00           C
ATOM     72  CA  ALA C  16      -1.865 -12.529 -28.077  1.00  1.00           C
ATOM     77  CA  ALA C  17      -3.254 -16.028 -28.473  1.00  1.00           C
ATOM     82  CA  ALA C  18      -5.534 -14.456 -31.052  1.00  1.00           C
ATOM     87  CA  ALA C  19      -2.657 -15.038 -33.442  1.00  1.00           C
ATOM     92  CA  ALA C  20      -4.146 -18.495 -33.840  1.00  1.00           C
ATOM     97  CA  ALA C  21      -6.342 -16.867 -36.458  1.00  1.00           C
ATOM    102  CA  ALA C  22      -3.451 -17.549 -38.805  1.00  1.00           C
"""
  print "Finding alpha,3-10 and pi helices...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 16  Start: 1 End: 16
Class:  Alpha helix  N: 16 Start: 1 End: 16  Rise: 1.51 A Dot: 0.98

Model 2  N: 21  Start: 2 End: 22
Class:     Pi helix  N: 21 Start: 2 End: 22  Rise: 0.96 A Dot: 0.98

Model 3  N: 21  Start: 2 End: 22
Class:   3-10 helix  N: 20 Start: 2 End: 21  Rise: 1.99 A Dot: 1.00

PDB RECORDS:
HELIX    1   1 ALA A    1  ALA A   16  1                                  16
HELIX    1   1 ALA C    2  ALA C   21  5                                  20
HELIX    1   1 ALA B    2  ALA B   22  3                                  21



PDB Selections:
" ( chain 'A' and resid 1  through 16 )  or  ( chain 'C' and resid 2  through 21 )  or  ( chain 'B' and resid 2  through 22 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  assert not test_utils.show_diff(found_text, expected_text)
  # if remove_blank(found_text)!=remove_blank(expected_text):
  #   print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
  #   raise AssertionError, "FAILED"
  print "OK"

def tst_04():
  text="""
ATOM      2  CA  THRAa   3     186.743 125.884 251.259  1.00100.00           C
ATOM      5  CA  ASNAa   4     189.629 123.742 252.763  1.00100.00           C
ATOM      8  CA  SERAa   5     191.072 126.112 255.320  1.00100.00           C
ATOM     11  CA  ASPAa   6     192.080 124.928 258.848  1.00100.00           C
ATOM     14  CA  PHEAa   7     189.384 124.585 261.530  1.00100.00           C
ATOM     17  CA  VALAa   8     189.248 124.466 265.315  1.00100.00           C
ATOM     20  CA  VALAa   9     187.059 122.294 267.547  1.00100.00           C
ATOM     23  CA  ILEAa  10     185.534 123.893 270.679  1.00100.00           C
ATOM     26  CA  LYSAa  11     183.570 122.134 273.450  1.00100.00           C
ATOM     29  CA  ALAAa  12     181.897 124.298 276.085  1.00100.00           C
ATOM     32  CA  LEUAa  13     182.733 123.145 279.601  1.00100.00           C
ATOM     35  CA  GLUAa  14     180.241 125.609 281.156  1.00100.00           C
ATOM     38  CA  ASPAa  15     177.155 127.540 279.985  1.00100.00           C
ATOM     41  CA  GLYAa  16     177.637 130.843 278.162  1.00100.00           C
ATOM     44  CA  VALAa  17     180.958 130.212 276.395  1.00100.00           C
ATOM     47  CA  ASNAa  18     181.477 132.715 273.547  1.00100.00           C
ATOM     50  CA  VALAa  19     183.320 131.753 270.320  1.00100.00           C
ATOM     53  CA  ILEAa  20     184.043 135.156 268.674  1.00100.00           C
ATOM     56  CA  GLYAa  21     185.054 135.558 264.994  1.00100.00           C
ATOM     59  CA  LEUAa  22     187.345 138.529 264.419  1.00100.00           C
ATOM     62  CA  THRAa  23     187.310 140.218 261.033  1.00100.00           C
ATOM     65  CA  ARGAa  24     189.831 139.523 258.335  1.00100.00           C
ATOM     68  CA  GLYAa  25     191.359 142.673 256.805  1.00100.00           C
ATOM     71  CA  ALAAa  26     192.794 146.041 257.837  1.00100.00           C
ATOM     74  CA  ASPAa  27     190.126 146.289 260.564  1.00100.00           C
ATOM     77  CA  THRAa  28     189.912 143.928 263.570  1.00100.00           C
ATOM     80  CA  ARGAa  29     186.413 143.856 265.033  1.00100.00           C
ATOM     83  CA  PHEAa  30     183.873 141.240 266.091  1.00100.00           C
ATOM     86  CA  HISAa  31     181.625 140.079 263.343  1.00100.00           C
ATOM     89  CA  HISAa  32     179.931 137.209 265.203  1.00100.00           C
ATOM     92  CA  SERAa  33     179.805 135.702 268.677  1.00100.00           C
ATOM     95  CA  GLUAa  34     178.501 132.109 268.857  1.00100.00           C
ATOM     98  CA  CYSAa  35     177.222 131.284 272.342  1.00100.00           C
ATOM    101  CA  LEUAa  36     177.646 127.700 273.502  1.00100.00           C
ATOM    104  CA  ASPAa  37     175.969 125.990 276.438  1.00100.00           C
ATOM    107  CA  LYSAa  38     177.682 123.298 278.488  1.00100.00           C
ATOM    110  CA  GLYAa  39     178.623 120.300 276.385  1.00100.00           C
ATOM    113  CA  GLUAa  40     177.892 121.761 272.941  1.00100.00           C
ATOM    116  CA  VALAa  41     180.597 121.439 270.276  1.00100.00           C
ATOM    119  CA  LEUAa  42     181.492 123.998 267.594  1.00100.00           C
ATOM    122  CA  ILEAa  43     183.793 123.155 264.645  1.00100.00           C
ATOM    125  CA  ALAAa  44     184.701 126.388 262.889  1.00100.00           C
ATOM    128  CA  GLNAa  45     186.987 127.209 259.959  1.00100.00           C
ATOM    131  CA  PHEAa  46     189.115 130.161 259.157  1.00100.00           C
ATOM    134  CA  THRAa  47     187.356 131.901 256.203  1.00100.00           C
ATOM    137  CA  GLUAa  48     187.180 134.953 253.965  1.00100.00           C
ATOM    140  CA  HISAa  49     185.578 136.805 256.905  1.00100.00           C
ATOM    143  CA  THRAa  50     187.343 135.292 259.938  1.00100.00           C
ATOM    146  CA  SERAa  51     191.129 135.327 260.339  1.00100.00           C
ATOM    149  CA  ALAAa  52     191.231 135.094 264.170  1.00100.00           C
ATOM    152  CA  ILEAa  53     188.989 133.390 266.744  1.00100.00           C
ATOM    155  CA  LYSAa  54     188.770 134.368 270.428  1.00100.00           C
ATOM    158  CA  VALAa  55     187.303 131.970 273.016  1.00100.00           C
ATOM    161  CA  ARGAa  56     185.817 133.382 276.214  1.00100.00           C
ATOM    164  CA  GLYAa  57     184.672 131.065 278.997  1.00100.00           C
ATOM    167  CA  LYSAa  58     185.698 127.553 280.004  1.00100.00           C
ATOM    170  CA  ALAAa  59     186.172 125.294 276.966  1.00100.00           C
ATOM    173  CA  TYRAa  60     188.258 122.444 275.620  1.00100.00           C
ATOM    176  CA  ILEAa  61     189.863 123.277 272.265  1.00100.00           C
ATOM    179  CA  GLNAa  62     191.492 121.098 269.577  1.00100.00           C
ATOM    182  CA  THRAa  63     193.550 122.431 266.653  1.00100.00           C
ATOM    185  CA  ARGAa  64     196.271 121.116 264.358  1.00100.00           C
ATOM    188  CA  HISAa  65     198.826 122.305 266.995  1.00100.00           C
ATOM    191  CA  GLYAa  66     197.443 120.330 269.914  1.00100.00           C
ATOM    194  CA  VALAa  67     194.865 120.679 272.646  1.00100.00           C
ATOM    197  CA  ILEAa  68     194.232 123.486 275.120  1.00100.00           C
ATOM    200  CA  GLUAa  69     191.576 124.693 277.564  1.00100.00           C
ATOM    203  CA  SERAa  70     190.301 128.219 277.907  1.00100.00           C
ATOM    206  CA  GLUAa  71     189.167 129.249 281.377  1.00100.00           C
ATOM    209  CA  GLYAa  72     186.003 131.073 282.428  1.00100.00           C
"""
  print "Finding sheets...",
  import iotbx.pdb
  from cctbx.array_family import flex
  hierarchy=iotbx.pdb.input(source_info='text',
       lines=flex.split_lines(text)).construct_hierarchy()
  fss=find_secondary_structure(hierarchy=hierarchy,out=null_out())

  expected_text="""
Model 1  N: 70  Start: 3 End: 72
Class:  Beta strand  N: 10 Start: 3 End: 12  Rise: 3.32 A Dot: 0.88
Class:  Beta strand  N: 9 Start: 16 End: 24  Rise: 3.24 A Dot: 0.97
Class:  Beta strand  N: 4 Start: 27 End: 30  Rise: 3.34 A Dot: 0.95
Class:  Beta strand  N: 6 Start: 31 End: 36  Rise: 3.29 A Dot: 0.99
Class:  Beta strand  N: 8 Start: 40 End: 47  Rise: 3.30 A Dot: 0.96
Class:  Beta strand  N: 5 Start: 51 End: 55  Rise: 3.41 A Dot: 1.00
Class:  Beta strand  N: 6 Start: 58 End: 63  Rise: 3.41 A Dot: 0.96
Class:  Beta strand  N: 7 Start: 66 End: 72  Rise: 3.41 A Dot: 0.98

PDB RECORDS:
SHEET    1   1 3 HISAa  32  LEUAa  36  0
SHEET    2   1 3 VALAa  17  LEUAa  22 -1  N  GLYAa  21   O  HISAa  32
SHEET    3   1 3 ALAAa  52  VALAa  55 -1  N  LYSAa  54   O  ILEAa  20
SHEET    1   2 4 GLUAa  40  GLNAa  45  0
SHEET    2   2 4 PHEAa   7  ALAAa  12 -1  N  ALAAa  12   O  GLUAa  40
SHEET    3   2 4 LYSAa  58  THRAa  63 -1  N  GLNAa  62   O  VALAa   9
SHEET    4   2 4 GLYAa  66  GLUAa  71 -1  N  SERAa  70   O  ALAAa  59



PDB Selections:
" ( chain 'Aa' and resid 32  through 36 )  or  ( chain 'Aa' and resid 17  through 22 )  or  ( chain 'Aa' and resid 52  through 55 )  or  ( chain 'Aa' and resid 40  through 45 )  or  ( chain 'Aa' and resid 7  through 12 )  or  ( chain 'Aa' and resid 58  through 63 )  or  ( chain 'Aa' and resid 66  through 71 ) "
"""
  f=StringIO()
  fss.show_summary(out=f,verbose=True)
  found_text=f.getvalue()
  assert not test_utils.show_diff(found_text, expected_text)
  # if remove_blank(found_text)!=remove_blank(expected_text):
  #   print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
  #   raise AssertionError, "FAILED"
  print "OK"

if __name__=="__main__":
  tst_01()
  tst_02()
  tst_03()
  tst_04()
