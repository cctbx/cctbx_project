from __future__ import absolute_import, division, print_function
import time

from mmtbx.programs import holton_geometry_validation
from iotbx.cli_parser import run_program



pdb_str = """
ATOM    111  N  AASP A   8       7.955  -2.855   0.018  0.50  7.81           N
ATOM    112  CA AASP A   8       8.002  -4.252  -0.397  0.50  7.24           C
ATOM    113  C  AASP A   8       8.554  -5.084   0.753  0.50  8.19           C
ATOM    114  O  AASP A   8       8.480  -4.659   1.910  0.50  7.28           O
ATOM    115  CB AASP A   8       8.775  -4.435  -1.715  0.50  7.27           C
ATOM    116  CG AASP A   8      10.242  -4.030  -1.622  0.50  7.87           C
ATOM    117  OD1AASP A   8      10.746  -3.749  -0.513  0.50  7.71           O
ATOM    118  OD2AASP A   8      10.901  -4.004  -2.686  0.50  7.34           O
ATOM    119  N  BASP A   8       7.795  -2.908   0.329  0.50  6.89           N
ATOM    120  CA BASP A   8       7.797  -4.330   0.033  0.50  7.76           C
ATOM    121  C  BASP A   8       8.506  -5.085   1.155  0.50  7.71           C
ATOM    122  O  BASP A   8       8.982  -4.499   2.132  0.50  8.33           O
ATOM    123  CB BASP A   8       8.401  -4.597  -1.349  0.50  8.40           C
ATOM    124  CG BASP A   8       9.896  -4.295  -1.427  0.50  8.54           C
ATOM    125  OD1BASP A   8      10.507  -3.890  -0.409  0.50  8.52           O
ATOM    126  OD2BASP A   8      10.467  -4.478  -2.524  0.50  8.59           O
ATOM    143  N  AVAL A  10      11.698  -5.862   1.272  0.50  7.04           N
ATOM    144  CA AVAL A  10      13.074  -5.505   1.619  0.50  6.68           C
ATOM    145  C  AVAL A  10      13.196  -4.022   1.969  0.50  6.97           C
ATOM    146  O  AVAL A  10      14.268  -3.427   1.822  0.50  7.19           O
ATOM    147  CB AVAL A  10      14.068  -5.906   0.511  0.50  8.14           C
ATOM    148  CG1AVAL A  10      14.094  -7.426   0.333  0.50  7.45           C
ATOM    149  CG2AVAL A  10      13.711  -5.221  -0.792  0.50  7.14           C
ATOM    150  N  BVAL A  10      11.294  -6.157   1.579  0.50  9.47           N
ATOM    151  CA BVAL A  10      12.679  -5.802   1.861  0.50  8.32           C
ATOM    152  C  BVAL A  10      12.838  -4.313   2.188  0.50  7.84           C
ATOM    153  O  BVAL A  10      13.927  -3.767   2.057  0.50  7.90           O
ATOM    154  CB BVAL A  10      13.639  -6.267   0.749  0.50  8.37           C
ATOM    155  CG1BVAL A  10      13.659  -7.780   0.662  0.50  8.88           C
ATOM    156  CG2BVAL A  10      13.239  -5.655  -0.585  0.50  8.99           C
ATOM    157  N  AASN A  11      12.095  -3.409   2.406  0.50  7.18           N
ATOM    158  CA AASN A  11      12.072  -2.046   2.939  0.50  6.24           C
ATOM    159  C  AASN A  11      12.158  -0.957   1.869  0.50  6.76           C
ATOM    160  O  AASN A  11      12.518   0.183   2.185  0.50  6.62           O
ATOM    161  CB AASN A  11      13.136  -1.843   4.030  0.50  7.10           C
ATOM    162  CG AASN A  11      12.842  -0.668   4.936  0.50  6.86           C
ATOM    163  OD1AASN A  11      11.710  -0.464   5.371  0.50  6.49           O
ATOM    164  ND2AASN A  11      13.876   0.109   5.233  0.50  6.63           N
ATOM    165  N  BASN A  11      11.747  -3.647   2.585  0.50  7.49           N
ATOM    166  CA BASN A  11      11.781  -2.275   3.101  0.50  7.01           C
ATOM    167  C  BASN A  11      12.007  -1.219   2.016  0.50  7.03           C
ATOM    168  O  BASN A  11      12.615  -0.177   2.280  0.50  6.89           O
ATOM    169  CB BASN A  11      12.824  -2.145   4.230  0.50  7.54           C
ATOM    170  CG BASN A  11      12.593  -0.931   5.117  0.50  7.89           C
ATOM    171  OD1BASN A  11      11.461  -0.620   5.482  0.50  6.89           O
ATOM    172  ND2BASN A  11      13.673  -0.244   5.471  0.50  7.51           N
ATOM    173  N  ACYS A  12      11.818  -1.258   0.612  0.50  6.26           N
ATOM    174  CA ACYS A  12      11.832  -0.272  -0.466  0.50  5.29           C
ATOM    175  C  ACYS A  12      10.434   0.290  -0.679  0.50  5.74           C
ATOM    176  O  ACYS A  12       9.456  -0.464  -0.737  0.50  5.92           O
ATOM    177  CB ACYS A  12      12.327  -0.891  -1.776  0.50  5.75           C
ATOM    178  SG ACYS A  12      13.966  -0.380  -2.285  0.50  5.62           S
ATOM    179  N  BCYS A  12      11.494  -1.447   0.808  0.50  6.57           N
ATOM    180  CA BCYS A  12      11.622  -0.502  -0.294  0.50  7.46           C
ATOM    181  C  BCYS A  12      10.280   0.166  -0.578  0.50  6.92           C
ATOM    182  O  BCYS A  12       9.254  -0.511  -0.685  0.50  6.31           O
ATOM    183  CB BCYS A  12      12.112  -1.223  -1.547  0.50  7.67           C
ATOM    184  SG BCYS A  12      13.695  -2.037  -1.328  0.50  8.82           S
ATOM    199  N  ATYR A  14       7.322   2.584  -2.627  0.50  6.23           N
ATOM    200  CA ATYR A  14       6.663   2.581  -3.933  0.50  5.87           C
ATOM    201  C  ATYR A  14       6.927   3.906  -4.639  0.50  6.21           C
ATOM    202  O  ATYR A  14       6.425   4.948  -4.212  0.50  5.65           O
ATOM    203  CB ATYR A  14       5.154   2.434  -3.753  0.50  5.82           C
ATOM    204  CG ATYR A  14       4.624   1.049  -3.457  0.50  5.69           C
ATOM    205  CD1ATYR A  14       4.972   0.371  -2.294  0.50  6.79           C
ATOM    206  CD2ATYR A  14       3.702   0.452  -4.312  0.50  6.13           C
ATOM    207  CE1ATYR A  14       4.446  -0.878  -2.018  0.50  6.09           C
ATOM    208  CE2ATYR A  14       3.177  -0.791  -4.044  0.50  6.42           C
ATOM    209  CZ ATYR A  14       3.548  -1.451  -2.897  0.50  5.51           C
ATOM    210  OH ATYR A  14       3.023  -2.697  -2.637  0.50  6.12           O
ATOM    211  N  BTYR A  14       7.325   2.527  -2.608  0.50  6.17           N
ATOM    212  CA BTYR A  14       6.638   2.433  -3.897  0.50  6.72           C
ATOM    213  C  BTYR A  14       6.782   3.752  -4.647  0.50  5.79           C
ATOM    214  O  BTYR A  14       6.167   4.754  -4.270  0.50  6.16           O
ATOM    215  CB BTYR A  14       5.148   2.159  -3.689  0.50  6.80           C
ATOM    216  CG BTYR A  14       4.761   0.731  -3.368  0.50  6.89           C
ATOM    217  CD1BTYR A  14       5.150   0.132  -2.175  0.50  6.18           C
ATOM    218  CD2BTYR A  14       3.957  -0.005  -4.236  0.50  6.90           C
ATOM    219  CE1BTYR A  14       4.783  -1.175  -1.868  0.50  8.20           C
ATOM    220  CE2BTYR A  14       3.582  -1.312  -3.933  0.50  6.96           C
ATOM    221  CZ BTYR A  14       3.997  -1.890  -2.748  0.50  6.78           C
ATOM    222  OH BTYR A  14       3.625  -3.181  -2.445  0.50  7.07           O
ATOM    337  N  ACYS A  22      -1.297   5.808  -6.189  0.50  6.56           N
ATOM    338  CA ACYS A  22      -1.078   5.393  -4.808  0.50  6.42           C
ATOM    339  C  ACYS A  22      -2.383   5.027  -4.114  0.50  6.35           C
ATOM    340  O  ACYS A  22      -2.398   4.124  -3.271  0.50  6.60           O
ATOM    341  CB ACYS A  22      -0.353   6.485  -4.020  0.50  6.14           C
ATOM    342  SG ACYS A  22       1.349   6.836  -4.545  0.50  5.40           S
ATOM    343  N  BCYS A  22      -1.301   5.667  -6.228  0.50  6.21           N
ATOM    344  CA BCYS A  22      -1.071   5.339  -4.825  0.50  6.37           C
ATOM    345  C  BCYS A  22      -2.368   5.018  -4.110  0.50  6.36           C
ATOM    346  O  BCYS A  22      -2.415   4.102  -3.282  0.50  6.57           O
ATOM    347  CB BCYS A  22      -0.368   6.485  -4.117  0.50  6.21           C
ATOM    348  SG BCYS A  22       1.310   6.624  -4.618  0.50  5.77           S
ATOM    461  N  ALYS A  30      -5.499  -2.976   1.528  0.50  9.91           N
ATOM    462  CA ALYS A  30      -6.635  -2.839   2.439  0.50 10.71           C
ATOM    463  C  ALYS A  30      -6.628  -1.571   3.288  0.50  9.70           C
ATOM    464  O  ALYS A  30      -7.236  -1.549   4.362  0.50  8.81           O
ATOM    465  CB ALYS A  30      -6.820  -4.083   3.317  0.50 10.52           C
ATOM    466  CG ALYS A  30      -7.599  -5.189   2.640  0.50 12.51           C
ATOM    467  CD ALYS A  30      -7.803  -6.365   3.568  0.50 10.25           C
ATOM    468  CE ALYS A  30      -7.876  -7.644   2.765  0.50 10.80           C
ATOM    469  NZ ALYS A  30      -6.655  -7.773   1.930  0.50 10.53           N
ATOM    470  N  BLYS A  30      -5.855  -3.072   1.681  0.50 10.64           N
ATOM    471  CA BLYS A  30      -6.946  -2.847   2.624  0.50  9.92           C
ATOM    472  C  BLYS A  30      -6.791  -1.543   3.390  0.50  9.35           C
ATOM    473  O  BLYS A  30      -7.488  -1.330   4.398  0.50  9.45           O
ATOM    474  CB BLYS A  30      -7.162  -4.032   3.575  0.50 11.35           C
ATOM    475  CG BLYS A  30      -7.228  -5.374   2.860  0.50 11.31           C
ATOM    476  CD BLYS A  30      -8.402  -5.359   1.913  0.50  9.67           C
ATOM    477  CE BLYS A  30      -8.733  -6.731   1.364  0.50  9.32           C
ATOM    478  NZ BLYS A  30      -9.986  -6.662   0.557  0.50 10.04           N
ATOM    487  N  AGLU A  32      -7.365   2.739   3.840  0.50  8.10           N
ATOM    488  CA AGLU A  32      -8.473   3.648   3.543  0.50  7.98           C
ATOM    489  C  AGLU A  32      -8.288   4.335   2.195  0.50  6.88           C
ATOM    490  O  AGLU A  32      -9.216   4.387   1.376  0.50  6.45           O
ATOM    491  CB AGLU A  32      -8.595   4.711   4.634  0.50  7.29           C
ATOM    492  CG AGLU A  32      -9.613   5.805   4.308  0.50  7.60           C
ATOM    493  CD AGLU A  32      -9.391   7.074   5.105  0.50  7.23           C
ATOM    494  OE1AGLU A  32      -8.512   7.868   4.713  0.50  8.01           O
ATOM    495  OE2AGLU A  32     -10.081   7.283   6.122  0.50  7.79           O
ATOM    496  N  BGLU A  32      -7.300   2.584   3.633  0.50  6.66           N
ATOM    497  CA BGLU A  32      -8.443   3.417   3.269  0.50  6.75           C
ATOM    498  C  BGLU A  32      -8.189   4.186   1.980  0.50  6.70           C
ATOM    499  O  BGLU A  32      -8.982   4.121   1.029  0.50  6.92           O
ATOM    500  CB BGLU A  32      -8.718   4.419   4.384  0.50  7.54           C
ATOM    501  CG BGLU A  32      -9.756   5.449   3.995  0.50  6.77           C
ATOM    502  CD BGLU A  32      -9.626   6.748   4.747  0.50  7.18           C
ATOM    503  OE1BGLU A  32      -8.610   7.447   4.548  0.50  7.36           O
ATOM    504  OE2BGLU A  32     -10.546   7.076   5.526  0.50  7.04           O
ATOM    667  N  AGLY A  43      12.501   9.989  -2.091  0.50  7.99           N
ATOM    668  CA AGLY A  43      11.943  10.332  -3.372  0.50  8.26           C
ATOM    669  C  AGLY A  43      10.431  10.232  -3.337  0.50  7.19           C
ATOM    670  O  AGLY A  43       9.797  10.095  -2.281  0.50  6.81           O
ATOM    671  N  BGLY A  43      12.410   9.938  -1.840  0.50  7.47           N
ATOM    672  CA BGLY A  43      12.031  10.247  -3.203  0.50  7.28           C
ATOM    673  C  BGLY A  43      10.526  10.359  -3.311  0.50  6.96           C
ATOM    674  O  BGLY A  43       9.836  10.572  -2.316  0.50  7.68           O
ATOM    701  N  ACYS A  46       4.585   6.467  -2.227  0.50  5.36           N
ATOM    702  CA ACYS A  46       3.189   6.387  -1.826  0.50  5.02           C
ATOM    703  C  ACYS A  46       3.063   6.397  -0.308  0.50  5.06           C
ATOM    704  O  ACYS A  46       3.836   5.742   0.403  0.50  5.83           O
ATOM    705  CB ACYS A  46       2.532   5.126  -2.400  0.50  5.32           C
ATOM    706  SG ACYS A  46       2.338   5.087  -4.206  0.50  4.49           S
ATOM    707  N  BCYS A  46       4.534   6.419  -2.240  0.50  5.44           N
ATOM    708  CA BCYS A  46       3.142   6.346  -1.807  0.50  5.05           C
ATOM    709  C  BCYS A  46       3.037   6.367  -0.290  0.50  5.09           C
ATOM    710  O  BCYS A  46       3.725   5.604   0.403  0.50  5.43           O
ATOM    711  CB BCYS A  46       2.463   5.062  -2.291  0.50  5.52           C
ATOM    712  SG BCYS A  46       2.136   4.845  -4.042  0.50  5.43           S
ATOM    773  N  ALYS A  50      -4.766   4.450   7.185  0.50  5.92           N
ATOM    774  CA ALYS A  50      -5.598   3.384   7.750  0.50  6.92           C
ATOM    775  C  ALYS A  50      -5.256   1.999   7.200  0.50  6.62           C
ATOM    776  O  ALYS A  50      -6.133   1.212   6.829  0.50  6.24           O
ATOM    777  CB ALYS A  50      -7.089   3.703   7.631  0.50  6.66           C
ATOM    778  CG ALYS A  50      -7.476   4.993   8.327  0.50  7.19           C
ATOM    779  CD ALYS A  50      -8.989   5.149   8.428  0.50  7.18           C
ATOM    780  CE ALYS A  50      -9.356   6.389   9.224  0.50  7.39           C
ATOM    781  NZ ALYS A  50     -10.788   6.413   9.617  0.50  7.71           N
ATOM    782  N  BLYS A  50      -4.605   4.533   7.317  0.50  6.77           N
ATOM    783  CA BLYS A  50      -5.393   3.489   7.979  0.50  6.78           C
ATOM    784  C  BLYS A  50      -5.189   2.100   7.365  0.50  6.66           C
ATOM    785  O  BLYS A  50      -6.143   1.382   7.056  0.50  7.37           O
ATOM    786  CB BLYS A  50      -6.873   3.875   8.060  0.50  7.42           C
ATOM    787  CG BLYS A  50      -7.111   5.205   8.763  0.50  6.97           C
ATOM    788  CD BLYS A  50      -8.584   5.462   9.036  0.50  7.73           C
ATOM    789  CE BLYS A  50      -8.793   6.865   9.613  0.50  8.24           C
ATOM    790  NZ BLYS A  50      -8.648   7.937   8.584  0.50  7.95           N
"""
def run():
  filename = 'tst_holton_geometry_validation.pdb'
  f = open(filename,'w')
  print(pdb_str, file = f)
  f.close()
  result = run_program(program_class=holton_geometry_validation.Program,
   args = [filename])
  # print (result)
  expected_worst_table = {'CBETADEV': ['CBETADEV 8.1796   10 | B VAL A   10'], 'CLASH': ['CLASH 1.0814 -0.413 |  HB2 A GLU  A   32 -  HB2 A LYS  A   50 '], 'OMEGA': ['OMEGA 4.472084 10 0 171.49 |  CA  B VAL  A   10'], 'RAMA': ['RAMA 3.4214   11 | 0.89 B ASN A   11'], 'ROTA': ['ROTA 5.4138   12 | 0.08 A CYS A   12'], 'ANGLE': ['ANGLE 4.70 -4.99 119.39 114.40 2.3 |  CA  B CYS  A   46  -  CB  B CYS  A   46  -  SG  B CYS  A   46 '], 'BOND': ['BOND 2.42 0.051 1.757 1.808 0.03 |  CB  B CYS  A   22 -   SG  B CYS  A   22 '], 'CHIR': ['CHIR 1.315 0.23 2.28 2.51 0.2 |  CA  A LYS  A   30  -  N   A LYS  A   30  -  C   A LYS  A   30  -  CB  A LYS  A   30'], 'TORSION': ['TORSION 8.9570 -80.78 80.78 0.00 30.0 |  CB  A GLU  A   32  -  CG  A GLU  A   32  -  CD  A GLU  A   32  -  OE1 A GLU  A   32'], 'NONBOND': ['NONBOND 1.743667 0.513 2.927 3.440 1 |  O   B VAL  A   10  -  CB  B ASN  A   11'], 'PLANE': ['PLANE 1.4400 -0.024 0.024 0 0.02 |  CG  A TYR  A   14']}

  expected_result_table = {'CBETADEV': ['CBETADEV', '20', '3.25', '8.18', '1.0000', '0.7974', '6.5225', '3.2475'], 'CLASH': ['CLASH', '1', '1.08', '2.08', '0.7016', '0.6177', '1.9033', '1.4603'], 'OMEGA': ['OMEGA', '4', '2.16', '4.47', '0.9289', '0.7053', '3.1541', '2.0039'], 'RAMA': ['RAMA', '2', '3.23', '3.42', '0.9604', '0.6856', '2.3459', '3.1020'], 'ROTA': ['ROTA', '20', '0.82', '5.41', '0.3079', '0.6033', '3.2660', '0.2523'], 'ANGLE': ['ANGLE', '176', '0.35', '4.70', '0.0000', '0.2904', '1.3655', '0.0000'], 'BOND': ['BOND', '154', '0.08', '2.42', '0.0000', '0.1301', '0.3144', '0.0000'], 'CHIR': ['CHIR', '22', '0.32', '1.31', '0.0010', '0.1638', '0.2154', '0.0003'], 'TORSION': ['TORSION', '94', '1.57', '8.96', '0.9996', '0.6876', '6.1589', '1.5671'], 'FULL_NONBOND': ['FULL_NONBOND', '531', '-0.45', 'None', '0.0000', '0.0000', '0.0000', '0.0000'], 'NONBOND': ['NONBOND', '531', '1.05', '1.74', '0.6123', '0.3741', '0.6523', '0.6412'], 'PLANE': ['PLANE', '56', '0.07', '1.44', '0.0000', '0.1139', '0.1640', '0.0000']}

  if 0:
    print(result.worst_table)
    print(expected_worst_table)
  compare_tables(result.worst_table, expected_worst_table)
  if 0:
    print(result.result_table)
    print(expected_result_table)
  compare_tables(result.result_table, expected_result_table)


def compare_tables(a,b):
  a_keys = list(a.keys())
  b_keys = list(b.keys())
  a_keys.sort()
  b_keys.sort()
  assert a_keys == b_keys, (a_keys, b_keys)
  for key in a_keys:
    t1 = a[key]
    t2 = b[key]
    assert len(t1) == len(t2), (t1, t2)
    for x1, x2 in zip(t1, t2):
      try:
        xx1 = float(x1)
        xx2 = float(x2)
      except Exception as e:
        xx1 = 0
        xx2 = 0
      assert abs(xx1 - xx2) < 0.02, (xx1, xx2, x1,x2)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.3f"%(time.time()-t0))
  print("OK")
