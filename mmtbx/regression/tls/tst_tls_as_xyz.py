from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx import easy_run
import iotbx.pdb

pdb_str = """
REMARK   3  TLS DETAILS.
REMARK   3   NUMBER OF TLS GROUPS: 2
REMARK   3   ORIGIN: CENTER OF MASS
REMARK   3   TLS GROUP : 1
REMARK   3    SELECTION: CHAIN A
REMARK   3    ORIGIN FOR THE GROUP (A):   8.7338  28.3021  16.6793
REMARK   3    T TENSOR
REMARK   3      T11:   0.0695 T22:   0.0842
REMARK   3      T33:   0.0992 T12:   0.0039
REMARK   3      T13:  -0.0041 T23:   0.0014
REMARK   3    L TENSOR
REMARK   3      L11:   2.1725 L22:   0.5881
REMARK   3      L33:   1.1562 L12:   0.4879
REMARK   3      L13:  -0.0881 L23:   0.2404
REMARK   3    S TENSOR
REMARK   3      S11:   0.0512 S12:   0.0543 S13:   0.0992
REMARK   3      S21:  -0.0892 S22:  -0.0659 S23:  -0.0441
REMARK   3      S31:  -0.0042 S32:  -0.0302 S33:   0.0147
REMARK   3   TLS GROUP : 2
REMARK   3    SELECTION: CHAIN B
REMARK   3    ORIGIN FOR THE GROUP (A):  18.7338  48.3021  46.6793
REMARK   3    T TENSOR
REMARK   3      T11:   0.0695 T22:   0.0842
REMARK   3      T33:   0.0992 T12:   0.0039
REMARK   3      T13:  -0.0041 T23:   0.0014
REMARK   3    L TENSOR
REMARK   3      L11:   2.1725 L22:   0.5881
REMARK   3      L33:   1.1562 L12:   0.4879
REMARK   3      L13:  -0.0881 L23:   0.2404
REMARK   3    S TENSOR
REMARK   3      S11:   0.0512 S12:   0.0543 S13:   0.0992
REMARK   3      S21:  -0.0892 S22:  -0.0659 S23:  -0.0441
REMARK   3      S31:  -0.0042 S32:  -0.0302 S33:   0.0147
REMARK   3
CRYST1   35.050   40.500   42.370  90.00  90.00  90.00 P 21 21 21
SCALE1      0.028531  0.000000  0.000000        0.00000
SCALE2      0.000000  0.024691  0.000000        0.00000
SCALE3      0.000000  0.000000  0.023602        0.00000
ATOM      1  CA  THR A   6       6.096  14.546  15.382  1.00  0.00           C
ATOM      2  CA  THR A   7       4.643  17.379  17.472  1.00  0.00           C
ATOM      3  CA  TYR A   8       7.031  20.308  17.626  1.00  0.00           C
ATOM      4  CA  LYS A   9       6.715  23.191  20.108  1.00  0.00           C
ATOM      5  CA  LEU A  10       7.484  26.896  19.736  1.00  0.00           C
ATOM      6  CA  VAL A  11       8.313  28.973  22.814  1.00  0.00           C
ATOM      7  CA  ILE A  12       7.745  32.666  21.951  1.00  0.00           C
ATOM      8  CA  ASN A  13       9.367  35.454  23.951  1.00  0.00           C
ATOM      9  CA  GLY A  14       8.223  38.436  21.869  1.00  0.00           C
ATOM     10  CA  LYS A  15       7.732  42.071  22.635  1.00  0.00           C
ATOM     11  CA  THR A  16       3.910  41.794  22.427  1.00  0.00           C
ATOM     12  CA  LEU A  17       3.328  38.080  21.651  1.00  0.00           C
ATOM     13  CA  LYS A  18       4.268  35.598  24.351  1.00  0.00           C
ATOM     14  CA  GLY A  19       3.567  31.957  25.117  1.00  0.00           C
ATOM     15  CA  GLU A  20       3.694  28.633  23.293  1.00  0.00           C
ATOM     16  CA  THR A  21       2.207  27.042  20.183  1.00  0.00           C
ATOM     17  CA  THR A  22       2.609  23.707  18.385  1.00  0.00           C
ATOM     18  CA  THR A  23       2.682  22.183  14.916  1.00  0.00           C
ATOM     19  CA  LYS A  24       2.823  18.625  13.556  1.00  0.00           C
ATOM     20  CA  ALA A  25       5.312  17.719  10.843  1.00  0.00           C
ATOM     21  CA  VAL A  26       7.260  14.854   9.299  1.00  0.00           C
ATOM     22  CA  ASP A  27      10.612  16.526  10.001  1.00  0.00           C
ATOM     23  CA  ALA A  28      12.142  19.541  11.794  1.00  0.00           C
ATOM     24  CA  GLU A  29      12.712  21.542   8.632  1.00  0.00           C
ATOM     25  CA  THR A  30       9.005  21.472   7.817  1.00  0.00           C
ATOM     26  CA  ALA A  31       8.121  22.548  11.353  1.00  0.00           C
ATOM     27  CA  GLU A  32      10.719  25.317  11.217  1.00  0.00           C
ATOM     28  CA  LYS A  33       9.139  26.736   8.061  1.00  0.00           C
ATOM     29  CA  ALA A  34       5.654  26.682   9.658  1.00  0.00           C
ATOM     30  CA  PHE A  35       6.940  28.471  12.774  1.00  0.00           C
ATOM     31  CA  LYS A  36       8.915  31.065  10.827  1.00  0.00           C
ATOM     32  CA  GLN A  37       5.784  31.817   8.802  1.00  0.00           C
ATOM     33  CA  TYR A  38       3.757  32.125  12.022  1.00  0.00           C
ATOM     34  CA  ALA A  39       6.322  34.523  13.495  1.00  0.00           C
ATOM     35  CA  ASN A  40       6.375  36.639  10.306  1.00  0.00           C
ATOM     36  CA  ASP A  41       2.536  36.702  10.145  1.00  0.00           C
ATOM     37  CA  ASN A  42       2.519  38.093  13.718  1.00  0.00           C
ATOM     38  CA  GLY A  43       5.296  40.638  13.264  1.00  0.00           C
ATOM     39  CA  VAL A  44       7.921  38.823  15.338  1.00  0.00           C
ATOM     40  CA  ASP A  45      11.526  39.426  14.290  1.00  0.00           C
ATOM     41  CA  GLY A  46      13.865  37.643  16.623  1.00  0.00           C
ATOM     42  CA  VAL A  47      16.614  35.096  17.081  1.00  0.00           C
ATOM     43  CA  TRP A  48      16.026  31.364  17.129  1.00  0.00           C
ATOM     44  CA  THR A  49      17.146  28.133  18.719  1.00  0.00           C
ATOM     45  CA  TYR A  50      16.165  24.543  17.999  1.00  0.00           C
ATOM     46  CA  ASP A  51      16.716  21.634  20.410  1.00  0.00           C
ATOM     47  CA  ASP A  52      16.413  18.270  18.691  1.00  0.00           C
ATOM     48  CA  ALA A  53      16.219  16.432  22.049  1.00  0.00           C
ATOM     49  CA  THR A  54      12.971  18.177  23.062  1.00  0.00           C
ATOM     50  CA  LYS A  55      11.609  19.022  19.577  1.00  0.00           C
ATOM     51  CA  THR A  56      11.365  22.651  20.684  1.00  0.00           C
ATOM     52  CA  PHE A  57      12.052  25.886  18.846  1.00  0.00           C
ATOM     53  CA  THR A  58      12.461  29.177  20.755  1.00  0.00           C
ATOM     54  CA  VAL A  59      12.275  32.720  19.421  1.00  0.00           C
ATOM     55  CA  THR A  60      13.202  35.851  21.380  1.00  0.00           C
ATOM     56  CA  GLU A  61      12.989  39.521  20.535  1.00  0.00           C
TER
ATOM      1  CA  THR B   6      18.094  36.747  38.683  1.00  0.00           C
ATOM      2  CA  THR B   7      17.002  37.634  42.222  1.00  0.00           C
ATOM      3  CA  TYR B   8      18.786  40.734  43.451  1.00  0.00           C
ATOM      4  CA  LYS B   9      18.872  41.840  47.104  1.00  0.00           C
ATOM      5  CA  LEU B  10      18.852  45.310  48.658  1.00  0.00           C
ATOM      6  CA  VAL B  11      20.333  45.838  52.123  1.00  0.00           C
ATOM      7  CA  ILE B  12      18.910  49.102  53.545  1.00  0.00           C
ATOM      8  CA  ASN B  13      20.640  50.975  56.357  1.00  0.00           C
ATOM      9  CA  GLY B  14      18.382  54.043  56.519  1.00  0.00           C
ATOM     10  CA  LYS B  15      17.597  56.519  59.219  1.00  0.00           C
ATOM     11  CA  THR B  16      14.034  55.170  59.682  1.00  0.00           C
ATOM     12  CA  LEU B  17      13.836  52.292  57.150  1.00  0.00           C
ATOM     13  CA  LYS B  18      16.034  49.280  57.805  1.00  0.00           C
ATOM     14  CA  GLY B  19      16.242  45.699  56.592  1.00  0.00           C
ATOM     15  CA  GLU B  20      16.278  43.860  53.274  1.00  0.00           C
ATOM     16  CA  THR B  21      14.098  43.538  50.182  1.00  0.00           C
ATOM     17  CA  THR B  22      14.399  41.766  46.823  1.00  0.00           C
ATOM     18  CA  THR B  23      13.529  42.164  43.156  1.00  0.00           C
ATOM     19  CA  LYS B  24      13.775  39.919  40.086  1.00  0.00           C
ATOM     20  CA  ALA B  25      15.298  41.242  36.875  1.00  0.00           C
ATOM     21  CA  VAL B  26      17.040  40.230  33.662  1.00  0.00           C
ATOM     22  CA  ASP B  27      20.109  42.345  34.455  1.00  0.00           C
ATOM     23  CA  ALA B  28      21.647  44.472  37.241  1.00  0.00           C
ATOM     24  CA  GLU B  29      20.766  47.787  35.639  1.00  0.00           C
ATOM     25  CA  THR B  30      17.068  46.930  35.697  1.00  0.00           C
ATOM     26  CA  ALA B  31      17.284  45.873  39.341  1.00  0.00           C
ATOM     27  CA  GLU B  32      19.190  49.044  40.204  1.00  0.00           C
ATOM     28  CA  LYS B  33      16.417  51.191  38.731  1.00  0.00           C
ATOM     29  CA  ALA B  34      13.747  49.286  40.715  1.00  0.00           C
ATOM     30  CA  PHE B  35      15.711  49.704  43.961  1.00  0.00           C
ATOM     31  CA  LYS B  36      16.449  53.384  43.382  1.00  0.00           C
ATOM     32  CA  GLN B  37      12.736  53.956  42.785  1.00  0.00           C
ATOM     33  CA  TYR B  38      11.911  52.051  45.988  1.00  0.00           C
ATOM     34  CA  ALA B  39      14.398  54.150  47.965  1.00  0.00           C
ATOM     35  CA  ASN B  40      13.011  57.407  46.510  1.00  0.00           C
ATOM     36  CA  ASP B  41       9.393  56.311  47.200  1.00  0.00           C
ATOM     37  CA  ASN B  42      10.372  55.772  50.867  1.00  0.00           C
ATOM     38  CA  GLY B  43      12.371  58.965  51.313  1.00  0.00           C
ATOM     39  CA  VAL B  44      15.806  57.334  51.476  1.00  0.00           C
ATOM     40  CA  ASP B  45      18.685  59.472  50.212  1.00  0.00           C
ATOM     41  CA  GLY B  46      21.939  57.654  50.662  1.00  0.00           C
ATOM     42  CA  VAL B  47      25.055  56.218  49.087  1.00  0.00           C
ATOM     43  CA  TRP B  48      25.136  52.936  47.216  1.00  0.00           C
ATOM     44  CA  THR B  49      27.244  49.886  46.524  1.00  0.00           C
ATOM     45  CA  TYR B  50      26.676  46.957  44.186  1.00  0.00           C
ATOM     46  CA  ASP B  51      28.485  43.605  44.454  1.00  0.00           C
ATOM     47  CA  ASP B  52      28.165  41.547  41.287  1.00  0.00           C
ATOM     48  CA  ALA B  53      29.434  38.394  43.060  1.00  0.00           C
ATOM     49  CA  THR B  54      26.490  38.319  45.499  1.00  0.00           C
ATOM     50  CA  LYS B  55      23.900  40.218  43.402  1.00  0.00           C
ATOM     51  CA  THR B  56      23.461  42.607  46.326  1.00  0.00           C
ATOM     52  CA  PHE B  57      22.940  46.353  46.449  1.00  0.00           C
ATOM     53  CA  THR B  58      23.434  48.296  49.708  1.00  0.00           C
ATOM     54  CA  VAL B  59      22.228  51.780  50.588  1.00  0.00           C
ATOM     55  CA  THR B  60      23.245  53.732  53.695  1.00  0.00           C
ATOM     56  CA  GLU B  61      22.160  57.082  55.047  1.00  0.00           C
TER
END
"""

def run(file_name = "tst_tls_as_xyz.pdb"):
  of = open(file_name,"w")
  print(pdb_str, file=of)
  of.close()
  uc = iotbx.pdb.input(file_name=file_name).crystal_symmetry().unit_cell()
  #for n in range(10,100,10)+range(100,1000,100)+range(1000,10001,1000)+[15000,20000]:
  for n in [1000,]:
    assert not easy_run.call("phenix.tls_as_xyz %s n_models=%s > tst_tls_as_xyz.log"%(
      file_name,str(n)))
    for i in [0,1]:
      u1 = iotbx.pdb.input(file_name=
        "tst_tls_as_xyz_u_from_ensemble_%s.pdb"%str(i)
        ).xray_structure_simple().scatterers().extract_u_cart(uc)
      u2 = iotbx.pdb.input(file_name=
        "tst_tls_as_xyz_u_from_tls_%s.pdb"%str(i)
        ).xray_structure_simple().scatterers().extract_u_cart(uc)

    u1, u2 = u1.as_double(), u2.as_double()
    cc = flex.linear_correlation(x=u1, y=u2).coefficient()
    r = flex.sum(flex.abs(flex.abs(u1)-flex.abs(u2)))/\
        flex.sum(flex.abs(flex.abs(u1)+flex.abs(u2)))*2
    print("%5d %6.4f %6.4f"%(n, cc, r))
  assert cc>0.99, cc
  assert r<0.06, r

if (__name__ == "__main__"):
  run()
