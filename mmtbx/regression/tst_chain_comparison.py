from __future__ import division

from cStringIO import StringIO
from mmtbx.validation.chain_comparison import run

def remove_blank(text):
  return text.replace(" ","").replace("\n","")

model="""
CRYST1  113.949  113.949   32.474  90.00  90.00  90.00 I 4
ATOM      9  CA  LYS     4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP     5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL     6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET     7     108.936  26.802  39.218  1.00 45.06      P9
ATOM     40  CA  SER     8     105.419  28.211  39.036  1.00 41.81      P9
ATOM     46  CA  THR     9     103.906  27.616  35.600  1.00 37.71      P9
ATOM     53  CA  LYS    10     103.277  30.858  33.696  1.00 33.20      P9
ATOM     62  CA  TYR    11     100.215  31.269  31.467  1.00 29.03      P9
ATOM     74  CA  VAL    12      98.628  33.636  28.957  1.00 24.81      P9
ATOM     81  CA  GLU    13      95.372  33.451  27.024  1.00 21.62      P9
ATOM     90  CA  ALA    14      95.648  31.766  23.627  1.00 21.34      P9
ATOM     95  CA  GLY    15      94.222  34.881  22.004  1.00 22.33      P9
ATOM     99  CA  GLU    16      97.221  36.890  23.214  1.00 23.17      P9
ATOM    108  CA  LEU    17      99.807  34.753  21.414  1.00 23.42      P9
ATOM    116  CA  LYS    18     101.429  36.062  18.237  1.00 24.53      P9
ATOM    125  CA  GLU    19     103.948  34.832  15.691  1.00 24.34      P9
ATOM    134  CA  GLY    20     107.227  34.699  17.578  1.00 22.47      P9
ATOM    138  CA  SER    21     105.632  33.977  20.956  1.00 20.20      P9
ATOM    144  CA  TYR    22     106.668  30.890  22.912  1.00 19.09      P9
ATOM    156  CA  VAL    23     104.283  28.225  24.176  1.00 18.63      P9
ATOM    163  CA  VAL    24     104.403  24.711  25.602  1.00 16.73      P9
ATOM    170  CA  ILE    25     102.485  22.073  23.645  1.00 16.21      P9
ATOM    178  CA  ASP    26     102.424  18.531  25.042  1.00 17.02      P9
ATOM    186  CA  GLY    27     105.535  19.131  27.135  1.00 19.18      P9
ATOM    190  CA  GLU    28     107.647  20.739  24.398  1.00 18.42      P9
ATOM    200  CA  PRO    29     108.490  24.457  24.454  1.00 17.13      P9
ATOM    206  CA  CYS    30     107.738  25.793  20.958  1.00 18.27      P9
ATOM    212  CA  ARG    31     108.013  29.010  18.950  1.00 18.81      P9
ATOM    223  CA  VAL    32     104.720  30.074  17.348  1.00 18.68      P9
ATOM    230  CA  VAL    33     104.990  30.324  13.554  1.00 21.62      P9
ATOM    237  CA  GLU    34     101.337  30.903  12.676  1.00 22.97      P9
ATOM    246  CA  ILE    35      97.990  31.567  14.325  1.00 23.13      P9
ATOM    254  CA  GLU    36      94.516  31.390  12.807  1.00 22.64      P9
ATOM    263  CA  LYS    37      91.314  32.656  14.422  1.00 21.25      P9
ATOM    272  CA  SER    38      87.754  31.489  13.758  1.00 20.23      P9
ATOM    278  CA  LYS    39      84.259  31.708  15.215  1.00 22.33      P9
ATOM    287  CA  THR    40      81.571  29.194  14.236  1.00 25.52      P9
ATOM    294  CA  GLY    41      79.184  29.756  17.114  0.93 29.60      P9
ATOM    298  CA  LYS    42      76.787  32.652  16.597  0.89 32.20      P9
ATOM    307  CA  HIS    43      77.170  33.509  20.288  1.00 31.87      P9
ATOM    317  CA  GLY    44      80.141  31.264  20.999  0.65 28.63      P9
ATOM    321  CA  SER    45      83.784  31.911  21.813  1.00 23.69      P9
ATOM    327  CA  ALA    46      86.373  32.582  19.132  1.00 18.89      P9
ATOM    332  CA  LYS    47      88.902  29.776  18.670  1.00 17.38      P9
ATOM    341  CA  ALA    48      92.611  29.914  17.899  1.00 17.02      P9
ATOM    346  CA  ARG    49      94.612  27.386  15.893  1.00 17.11      P9
ATOM    357  CA  ILE    50      98.284  27.635  16.768  1.00 17.68      P9
ATOM    365  CA  VAL    51     101.124  26.174  14.718  1.00 17.34      P9
ATOM    372  CA  ALA    52     104.518  26.094  16.413  1.00 17.71      P9
ATOM    377  CA  VAL    53     107.935  24.470  16.254  1.00 17.54      P9
ATOM    384  CA  GLY    54     109.704  22.649  19.055  1.00 18.35      P9
ATOM    388  CA  VAL    55     112.688  24.550  20.407  1.00 19.01      P9
ATOM    395  CA  PHE    56     114.509  21.294  21.108  1.00 19.90      P9
ATOM    406  CA  ASP    57     113.198  18.603  18.749  1.00 21.17      P9
ATOM    414  CA  GLY    58     112.293  20.976  15.923  1.00 21.69      P9
ATOM    418  CA  GLY    59     108.988  19.187  15.479  1.00 21.62      P9
ATOM    422  CA  LYS    60     105.864  20.872  14.151  1.00 20.30      P9
ATOM    431  CA  ARG    61     103.024  20.971  16.667  1.00 18.60      P9
ATOM    442  CA  THR    62      99.517  22.398  16.607  1.00 17.44      P9
ATOM    449  CA  LEU    63      96.998  23.440  19.242  1.00 16.19      P9
ATOM    457  CA  SER    64      93.370  24.446  18.716  1.00 15.39      P9
ATOM    463  CA  LEU    65      91.433  25.924  21.644  1.00 14.67      P9
ATOM    472  CA  PRO    66      89.260  28.905  22.730  1.00 14.79      P9
ATOM    478  CA  VAL    67      90.940  32.325  22.634  1.00 15.74      P9
ATOM    485  CA  ASP    68      90.325  32.738  26.369  1.00 17.10      P9
ATOM    493  CA  ALA    69      91.793  29.371  27.372  1.00 16.26      P9
ATOM    498  CA  GLN    70      95.115  29.425  29.223  1.00 18.00      P9
ATOM    507  CA  VAL    71      98.237  28.190  27.463  1.00 19.41      P9
ATOM    514  CA  GLU    72     101.563  27.596  29.201  1.00 21.29      P9
ATOM    523  CA  VAL    73     104.492  29.822  28.262  1.00 21.83      P9
ATOM    531  CA  PRO    74     108.047  28.603  28.963  0.94 22.19      P9
ATOM    537  CA  ILE    75     110.436  30.660  31.072  1.00 23.67      P9
ATOM    545  CA  ILE    76     113.508  31.336  28.943  0.96 24.46      P9
ATOM    553  CA  GLU    77     116.735  32.457  30.611  1.00 24.69      P9
ATOM    562  CA  LYS    78     119.294  33.935  28.213  1.00 24.12      P9
ATOM    571  CA  PHE    79     123.013  34.160  28.959  1.00 21.39      P9
ATOM    582  CA  THR    80     126.479  34.196  27.404  1.00 19.30      P9
ATOM    589  CA  ALA    81     128.958  31.320  27.436  1.00 17.03      P9
ATOM    594  CA  GLN    82     132.421  30.668  25.993  1.00 15.83      P9
ATOM    603  CA  ILE    83     133.128  27.536  23.964  1.00 14.35      P9
ATOM    611  CA  LEU    84     135.770  25.349  25.639  1.00 14.01      P9
ATOM    619  CA  SER    85     135.720  22.442  23.163  1.00 14.77      P9
ATOM    625  CA  VAL    86     133.697  21.017  20.292  1.00 15.76      P9
ATOM    632  CA  SER    87     133.246  17.344  19.410  1.00 15.97      P9
ATOM    638  CA  GLY    88     130.905  15.871  16.836  1.00 16.65      P9
ATOM    642  CA  ASP    89     128.352  15.243  19.587  1.00 17.36      P9
ATOM    650  CA  VAL    90     128.894  17.746  22.389  1.00 16.73      P9
ATOM    657  CA  ILE    91     129.780  21.401  22.870  1.00 16.05      P9
ATOM    665  CA  GLN    92     131.584  22.059  26.178  1.00 16.32      P9
ATOM    674  CA  LEU    93     130.870  25.590  27.465  1.00 15.69      P9
ATOM    682  CA  MET    94     131.982  27.904  30.270  1.00 17.32      P9
ATOM    690  CA  ASP    95     129.006  29.830  31.697  1.00 19.94      P9
ATOM    698  CA  MET    96     130.333  33.412  31.817  1.00 22.68      P9
ATOM    706  CA  ARG    97     128.065  34.204  34.767  1.00 24.59      P9
ATOM    717  CA  ASP    98     129.544  31.755  37.283  0.96 24.84      P9
ATOM    725  CA  TYR    99     132.285  30.135  35.175  0.99 24.06      P9
ATOM    737  CA  LYS   100     130.841  26.625  35.550  1.00 22.14      P9
ATOM    746  CA  THR   101     131.149  24.064  32.761  1.00 20.84      P9
ATOM    753  CA  ILE   102     127.988  23.035  30.916  1.00 19.78      P9
ATOM    761  CA  GLU   103     127.678  20.418  28.177  1.00 19.97      P9
ATOM    770  CA  VAL   104     125.251  20.963  25.299  1.00 18.00      P9
ATOM    778  CA  PRO   105     124.362  18.242  22.754  1.00 18.24      P9
ATOM    784  CA  MET   106     125.166  19.138  19.128  1.00 18.91      P9
ATOM    792  CA  LYS   107     121.491  18.558  18.303  1.00 19.94      P9
ATOM    801  CA  TYR   108     120.616  21.558  20.496  1.00 20.72      P9
ATOM    813  CA  VAL   109     122.621  24.018  18.411  1.00 23.20      P9
ATOM    820  CA  GLU   110     120.751  26.260  15.961  1.00 28.03      P9
ATOM    829  CA  GLU   111     121.704  25.099  12.460  0.74 31.85      P9
ATOM    838  CA  GLU   112     122.954  28.554  11.512  1.00 33.44      P9
ATOM    847  CA  ALA   113     125.368  28.615  14.458  1.00 31.90      P9
ATOM    852  CA  LYS   114     126.911  25.180  13.906  0.92 31.27      P9
ATOM    861  CA  GLY   115     129.136  26.579  11.174  1.00 30.17      P9
ATOM    865  CA  ARG   116     130.956  29.004  13.463  1.00 27.89      P9
ATOM    876  CA  LEU   117     131.365  26.777  16.514  1.00 25.28      P9
ATOM    884  CA  ALA   118     134.998  27.063  17.578  1.00 23.01      P9
ATOM    890  CA  PRO   119     136.912  26.997  20.869  1.00 21.21      P9
ATOM    896  CA  GLY   120     137.158  30.500  22.308  1.00 21.69      P9
ATOM    900  CA  ALA   121     134.046  31.875  20.620  1.00 21.58      P9
ATOM    905  CA  GLU   122     131.240  33.298  22.745  1.00 22.41      P9
ATOM    914  CA  VAL   123     127.678  32.132  22.205  1.00 20.97      P9
ATOM    921  CA  GLU   124     124.171  33.150  23.219  1.00 20.32      P9
ATOM    930  CA  VAL   125     122.627  30.323  25.229  1.00 19.82      P9
ATOM    937  CA  TRP   126     119.012  29.745  26.239  1.00 21.08      P9
ATOM    951  CA  GLN   127     118.027  27.643  29.216  1.00 22.22      P9
ATOM    960  CA  ILE   128     114.538  26.200  29.585  1.00 22.80      P9
ATOM    968  CA  LEU   129     114.153  24.029  32.671  1.00 24.59      P9
ATOM    976  CA  ASP   130     117.058  21.551  32.616  0.99 24.90      P9
ATOM    984  CA  ARG   131     117.952  21.918  28.938  1.00 22.57      P9
ATOM    995  CA  TYR   132     120.128  24.357  27.016  1.00 21.18      P9
ATOM   1007  CA  LYS   133     120.153  25.537  23.432  1.00 21.22      P9
ATOM   1016  CA  ILE   134     122.957  27.373  21.662  1.00 20.82      P9
ATOM   1024  CA  ILE   135     121.298  30.133  19.637  1.00 24.27      P9
ATOM   1032  CA  ARG   136     124.175  31.872  17.866  1.00 28.80      P9
ATOM   1043  CA  VAL   137     127.877  32.680  17.955  1.00 32.92      P9
ATOM   1050  CA  LYS   138     128.895  35.931  19.718  1.00 36.89      P9
ATOM   1059  CA  GLY   139     126.242  36.878  22.264  0.67 38.21      P9
"""

query="""
ATOM      2  CA  HIS A   7      47.968  37.350   1.343  1.00 30.00           C
ATOM     12  CA  ILE A   8      46.579  37.684   4.402  1.00 30.00           C
ATOM     20  CA  LYS A   9      47.784  38.702   7.227  1.00 30.00           C
ATOM     29  CA  LEU A  10      45.288  36.534   8.774  1.00 30.00           C
ANISOU   29  CA  LEU A  10     4726   8440   6836    240  -1242    212       C
ATOM     37  CA  MET A  11      43.852  36.569  12.293  1.00 30.00           C
ATOM     45  CA  ASN A  12      42.415  33.671  14.884  1.00 30.00           C
ATOM     53  CA  ALA U  31      39.285  33.668  16.254  1.00 30.00      UNK  C
ATOM     58  CA  CYS U  32      40.481  30.632  18.193  1.00 30.00      UNK  C
ATOM     64  CA  LYS U  33      38.694  27.705  19.208  1.00 30.00      UNK  C
ATOM     73  CA  GLY U  34      41.530  25.709  20.615  1.00 30.00      UNK  C
ATOM     77  CA  LEU U  35      42.598  23.279  18.933  1.00 30.00      UNK  C
ATOM     85  CA  GLY U  36      45.899  23.083  17.741  1.00 30.00      UNK  C
ATOM     89  CA  HIS U  37      47.964  25.323  17.214  1.00 30.00      UNK  C
ATOM     99  CA  TRP U  41      49.522  32.736   0.687  1.00 30.00      UNK  C
ATOM    113  CA  PRO U  42      46.846  33.404   2.158  1.00 30.00      UNK  C
ATOM    120  CA  GLY U  43      49.716  33.742   4.809  1.00 30.00      UNK  C
ATOM    124  CA  MET U  44      50.275  36.496   3.408  1.00 30.00      UNK  C
ATOM    132  CA  ARG U  51      66.471  31.618   0.165  1.00 30.00      UNK  C
ATOM    143  CA  THR U  52      62.673  32.286   0.000  1.00 30.00      UNK  C
ATOM    150  CA  TRP U  53      60.396  33.791   2.514  1.00 30.00      UNK  C
ATOM    164  CA  PRO U  54      56.975  33.518   4.059  1.00 30.00      UNK  C
ATOM    171  CA  GLY U  55      56.742  37.050   3.683  1.00 30.00      UNK  C
ATOM    175  CA  HIS U  56      58.241  37.385   0.677  1.00 30.00      UNK  C
ATOM    185  CA  GLY U  57      60.259  36.521  -2.061  1.00 30.00      UNK  C
ATOM    189  CA  TRP U  58      63.306  36.585  -0.658  1.00 30.00      UNK  C
ATOM    203  CA  HIS U  59      66.371  36.917  -1.353  1.00 30.00      UNK  C
ATOM    213  CA  ARG U  71      33.885  27.955   5.046  1.00 30.00      UNK  C
ATOM    224  CA  PRO U  72      34.112  31.066   3.070  1.00 30.00      UNK  C
ATOM    231  CA  GLY U  73      37.023  30.068   0.698  1.00 30.00      UNK  C
ATOM    235  CA  ASN U  74      40.330  31.142   0.443  1.00 30.00      UNK  C
ATOM    243  CA  GLY U  75      39.883  28.291  -1.998  1.00 30.00      UNK  C
ATOM    247  CA  ALA U  76      40.699  29.754  -4.736  1.00 30.00      UNK  C
ATOM    252  CA  TYR U  77      42.908  31.262  -3.328  1.00 30.00      UNK  C
ATOM    264  CA  GLY U  78      44.842  29.020  -1.084  1.00 30.00      UNK  C
ATOM    268  CA  ALA U  79      46.698  27.855  -3.383  1.00 30.00      UNK  C
ATOM    273  CA  VAL U  91      51.081  32.321   8.150  1.00 30.00      UNK  C
ATOM    280  CA  ARG U  92      51.757  33.221  11.457  1.00 30.00      UNK  C
ATOM    291  CA  TYR U  93      53.177  35.451  14.207  1.00 30.00      UNK  C
ATOM    303  CA  GLY U  94      55.930  33.240  16.341  1.00 30.00      UNK  C
ATOM    307  CA  LYS U  95      56.538  31.836  13.500  1.00 30.00      UNK  C
ATOM    316  CA  GLU U  96      53.605  29.086  13.128  1.00 30.00      UNK  C
ATOM    325  CA  TRP U  97      52.331  28.161   9.983  1.00 30.00      UNK  C
ATOM    339  CA  ASP U  98      49.378  26.104   9.729  1.00 30.00      UNK  C
ATOM    347  CA  GLU U  99      47.812  25.689   7.165  1.00 30.00      UNK  C
ATOM    356  CA  GLY U 111      34.962  41.965  -2.040  1.00 30.00      UNK  C
ATOM    360  CA  GLY U 112      34.552  41.382   1.896  1.00 30.00      UNK  C
ATOM    364  CA  PHE U 113      32.652  38.617   2.536  1.00 30.00      UNK  C
ATOM    375  CA  GLY U 114      33.035  35.430   1.690  1.00 30.00      UNK  C
ATOM    379  CA  LEU U 115      31.413  34.883   5.090  1.00 30.00      UNK  C
ATOM    387  CA  SER U 116      29.190  32.919   6.131  1.00 30.00      UNK  C
ATOM    393  CA  LEU U 117      25.231  33.701   6.955  1.00 30.00      UNK  C
ATOM    401  CA  PRO U 118      24.358  32.286  10.112  1.00 30.00      UNK  C
ATOM    408  CA  TYR U 119      26.935  34.425  10.677  1.00 30.00      UNK  C
ATOM    420  CA  GLY U 131      38.770  41.005   4.746  1.00 30.00      UNK  C
ATOM    424  CA  GLN U 132      39.116  39.652   2.001  1.00 30.00      UNK  C
ATOM    433  CA  TRP U 133      36.781  37.434   3.332  1.00 30.00      UNK  C
ATOM    447  CA  TYR U 134      37.051  34.918   6.022  1.00 30.00      UNK  C
ATOM    459  CA  ALA U 135      34.718  33.052   8.152  1.00 30.00      UNK  C
ATOM    464  CA  CYS U 136      35.953  32.045  10.889  1.00 30.00      UNK  C
ATOM    470  CA  GLY U 137      37.106  29.541   8.056  1.00 30.00      UNK  C
ATOM    474  CA  TYR U 138      39.250  29.269   5.412  1.00 30.00      UNK  C
"""

def tst_01():
  print "Comparing mixed model with target..."
  import iotbx.pdb
  from cctbx.array_family import flex
  model_pdb_inp=iotbx.pdb.input(source_info='model',
       lines=flex.split_lines(model))
  crystal_symmetry=model_pdb_inp.crystal_symmetry()
  model_hierarchy=model_pdb_inp.construct_hierarchy()
  query_hierarchy=iotbx.pdb.input(source_info='query',
       lines=flex.split_lines(query)).construct_hierarchy()

  f=StringIO()
  r=run(crystal_symmetry=crystal_symmetry,
    chain_hierarchy=query_hierarchy,target_hierarchy=model_hierarchy,out=f)
  expected_text="""
Residues matching in forward direction:     16  RMSD:   1.45
Residues matching in reverse direction:     31  RMSD:   1.40
Residues near but not matching one-to-one:  12  RMSD:   1.87

All residues near target:                   59  RMSD:   1.52
Residues far from target:                    2  RMSD:   2.04
"""

  found_text="\n".join(f.getvalue().splitlines()[-6:])
  if remove_blank(found_text)!=remove_blank(expected_text):
    print "Expected: \n%s \nFound: \n%s" %(expected_text,found_text)
    raise AssertionError, "FAILED"
  from libtbx.test_utils import approx_equal
  print r.get_values("forward")
  assert approx_equal(r.get_values("forward"),(1.4473857036049544, 16))
  print r.get_values("reverse")
  assert approx_equal(r.get_values("reverse"),(1.3969610738798282, 31))
  print r.get_values("close")
  assert approx_equal(r.get_values("close"),(1.5184018499613678, 59))
  print r.get_values("all_far")
  assert approx_equal(r.get_values("all_far"),(0,0))
  print "OK"


if __name__=="__main__":
  tst_01()
