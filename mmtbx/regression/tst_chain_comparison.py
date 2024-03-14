from __future__ import absolute_import, division, print_function

from six.moves import cStringIO as StringIO
from mmtbx.validation.chain_comparison import run

from iotbx.cli_parser import run_program



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

model1="""
CRYST1  113.949  113.949   32.474  90.00  90.00  90.00 I 4
ATOM      1  CA  LYS A   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM      2  CA  TRP A   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM      3  CA  VAL A   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM      4  CA  MET A   7     108.936  26.802  39.218  1.00 45.06      P9
ATOM      5  CA  SER A   8     105.419  28.211  39.036  1.00 41.81      P9
ATOM      6  CA  THR A   9     103.906  27.616  35.600  1.00 37.71      P9
ATOM      7  CA  LYS A  10     103.277  30.858  33.696  1.00 33.20      P9
ATOM      8  CA  TYR A  11     100.215  31.269  31.467  1.00 29.03      P9
ATOM      9  CA  VAL A  12      98.628  33.636  28.957  1.00 24.81      P9
ATOM     10  CA  GLU A  13      95.372  33.451  27.024  1.00 21.62      P9
ATOM     11  CA  ALA A  14      95.648  31.766  23.627  1.00 21.34      P9
ATOM     12  CA  GLY A  15      94.222  34.881  22.004  1.00 22.33      P9
ATOM     13  CA  GLU A  16      97.221  36.890  23.214  1.00 23.17      P9
ATOM     14  CA  LEU A  17      99.807  34.753  21.414  1.00 23.42      P9
ATOM     15  CA  LYS A  18     101.429  36.062  18.237  1.00 24.53      P9
ATOM     16  CA  GLU A  19     103.948  34.832  15.691  1.00 24.34      P9
ATOM     17  CA  GLY A  20     107.227  34.699  17.578  1.00 22.47      P9
ATOM     18  CA  SER A  21     105.632  33.977  20.956  1.00 20.20      P9
ATOM     19  CA  TYR A  22     106.668  30.890  22.912  1.00 19.09      P9
ATOM     20  CA  VAL A  23     104.283  28.225  24.176  1.00 18.63      P9
ATOM     21  CA  VAL A  24     104.403  24.711  25.602  1.00 16.73      P9
ATOM     22  CA  ILE A  25     102.485  22.073  23.645  1.00 16.21      P9
ATOM     23  CA  ASP A  26     102.424  18.531  25.042  1.00 17.02      P9
ATOM     24  CA  GLY A  27     105.535  19.131  27.135  1.00 19.18      P9
ATOM     25  CA  GLU A  28     107.647  20.739  24.398  1.00 18.42      P9
ATOM     26  CA  PRO A  29     108.490  24.457  24.454  1.00 17.13      P9
ATOM     27  CA  CYS A  30     107.738  25.793  20.958  1.00 18.27      P9
ATOM     28  CA  ARG A  31     108.013  29.010  18.950  1.00 18.81      P9
ATOM     29  CA  VAL A  32     104.720  30.074  17.348  1.00 18.68      P9
ATOM     30  CA  VAL A  33     104.990  30.324  13.554  1.00 21.62      P9
ATOM     31  CA  GLU A  34     101.337  30.903  12.676  1.00 22.97      P9
ATOM     32  CA  ILE A  35      97.990  31.567  14.325  1.00 23.13      P9
ATOM     33  CA  GLU A  36      94.516  31.390  12.807  1.00 22.64      P9
ATOM     34  CA  LYS A  37      91.314  32.656  14.422  1.00 21.25      P9
ATOM     35  CA  SER A  38      87.754  31.489  13.758  1.00 20.23      P9
ATOM     36  CA  LYS A  39      84.259  31.708  15.215  1.00 22.33      P9
ATOM     37  CA  THR A  40      81.571  29.194  14.236  1.00 25.52      P9
ATOM     38  CA  GLY A  41      79.184  29.756  17.114  0.93 29.60      P9
ATOM     39  CA  LYS A  42      76.787  32.652  16.597  0.89 32.20      P9
ATOM     40  CA  HIS A  43      77.170  33.509  20.288  1.00 31.87      P9
ATOM     41  CA  GLY A  44      80.141  31.264  20.999  0.65 28.63      P9
ATOM     42  CA  SER A  45      83.784  31.911  21.813  1.00 23.69      P9
ATOM     43  CA  ALA A  46      86.373  32.582  19.132  1.00 18.89      P9
ATOM     44  CA  LYS A  47      88.902  29.776  18.670  1.00 17.38      P9
ATOM     45  CA  ALA A  48      92.611  29.914  17.899  1.00 17.02      P9
ATOM     46  CA  ARG A  49      94.612  27.386  15.893  1.00 17.11      P9
ATOM     47  CA  ILE A  50      98.284  27.635  16.768  1.00 17.68      P9
ATOM     48  CA  VAL A  51     101.124  26.174  14.718  1.00 17.34      P9
ATOM     49  CA  ALA A  52     104.518  26.094  16.413  1.00 17.71      P9
ATOM     50  CA  VAL A  53     107.935  24.470  16.254  1.00 17.54      P9
ATOM     51  CA  GLY A  54     109.704  22.649  19.055  1.00 18.35      P9
ATOM     52  CA  VAL A  55     112.688  24.550  20.407  1.00 19.01      P9
ATOM     53  CA  PHE A  56     114.509  21.294  21.108  1.00 19.90      P9
ATOM     54  CA  ASP A  57     113.198  18.603  18.749  1.00 21.17      P9
ATOM     55  CA  GLY A  58     112.293  20.976  15.923  1.00 21.69      P9
ATOM     56  CA  GLY A  59     108.988  19.187  15.479  1.00 21.62      P9
ATOM     57  CA  LYS A  60     105.864  20.872  14.151  1.00 20.30      P9
ATOM     58  CA  ARG A  61     103.024  20.971  16.667  1.00 18.60      P9
ATOM     59  CA  THR A  62      99.517  22.398  16.607  1.00 17.44      P9
ATOM     60  CA  LEU A  63      96.998  23.440  19.242  1.00 16.19      P9
ATOM     61  CA  SER A  64      93.370  24.446  18.716  1.00 15.39      P9
ATOM     62  CA  LEU A  65      91.433  25.924  21.644  1.00 14.67      P9
ATOM     63  CA  PRO A  66      89.260  28.905  22.730  1.00 14.79      P9
ATOM     64  CA  VAL A  67      90.940  32.325  22.634  1.00 15.74      P9
ATOM     65  CA  ASP A  68      90.325  32.738  26.369  1.00 17.10      P9
ATOM     66  CA  ALA A  69      91.793  29.371  27.372  1.00 16.26      P9
ATOM     67  CA  GLN A  70      95.115  29.425  29.223  1.00 18.00      P9
ATOM     68  CA  VAL A  71      98.237  28.190  27.463  1.00 19.41      P9
ATOM     69  CA  GLU A  72     101.563  27.596  29.201  1.00 21.29      P9
ATOM     70  CA  VAL A  73     104.492  29.822  28.262  1.00 21.83      P9
ATOM     71  CA  PRO A  74     108.047  28.603  28.963  0.94 22.19      P9
ATOM     72  CA  ILE A  75     110.436  30.660  31.072  1.00 23.67      P9
ATOM     73  CA  ILE A  76     113.508  31.336  28.943  0.96 24.46      P9
ATOM     74  CA  GLU A  77     116.735  32.457  30.611  1.00 24.69      P9
ATOM     75  CA  LYS A  78     119.294  33.935  28.213  1.00 24.12      P9
ATOM     76  CA  PHE A  79     123.013  34.160  28.959  1.00 21.39      P9
ATOM     77  CA  THR A  80     126.479  34.196  27.404  1.00 19.30      P9
ATOM     78  CA  ALA A  81     128.958  31.320  27.436  1.00 17.03      P9
ATOM     79  CA  GLN A  82     132.421  30.668  25.993  1.00 15.83      P9
ATOM     80  CA  ILE A  83     133.128  27.536  23.964  1.00 14.35      P9
ATOM     81  CA  LEU A  84     135.770  25.349  25.639  1.00 14.01      P9
ATOM     82  CA  SER A  85     135.720  22.442  23.163  1.00 14.77      P9
ATOM     83  CA  VAL A  86     133.697  21.017  20.292  1.00 15.76      P9
ATOM     84  CA  SER A  87     133.246  17.344  19.410  1.00 15.97      P9
ATOM     85  CA  GLY A  88     130.905  15.871  16.836  1.00 16.65      P9
ATOM     86  CA  ASP A  89     128.352  15.243  19.587  1.00 17.36      P9
ATOM     87  CA  VAL A  90     128.894  17.746  22.389  1.00 16.73      P9
ATOM     88  CA  ILE A  91     129.780  21.401  22.870  1.00 16.05      P9
ATOM     89  CA  GLN A  92     131.584  22.059  26.178  1.00 16.32      P9
ATOM     90  CA  LEU A  93     130.870  25.590  27.465  1.00 15.69      P9
ATOM     91  CA  MET A  94     131.982  27.904  30.270  1.00 17.32      P9
ATOM     92  CA  ASP A  95     129.006  29.830  31.697  1.00 19.94      P9
ATOM     93  CA  MET A  96     130.333  33.412  31.817  1.00 22.68      P9
ATOM     94  CA  ARG A  97     128.065  34.204  34.767  1.00 24.59      P9
ATOM     95  CA  ASP A  98     129.544  31.755  37.283  0.96 24.84      P9
ATOM     96  CA  TYR A  99     132.285  30.135  35.175  0.99 24.06      P9
ATOM     97  CA  LYS A 100     130.841  26.625  35.550  1.00 22.14      P9
ATOM     98  CA  THR A 101     131.149  24.064  32.761  1.00 20.84      P9
ATOM     99  CA  ILE A 102     127.988  23.035  30.916  1.00 19.78      P9
ATOM    100  CA  GLU A 103     127.678  20.418  28.177  1.00 19.97      P9
ATOM    101  CA  VAL A 104     125.251  20.963  25.299  1.00 18.00      P9
ATOM    102  CA  PRO A 105     124.362  18.242  22.754  1.00 18.24      P9
ATOM    103  CA  MET A 106     125.166  19.138  19.128  1.00 18.91      P9
ATOM    104  CA  LYS A 107     121.491  18.558  18.303  1.00 19.94      P9
ATOM    105  CA  TYR A 108     120.616  21.558  20.496  1.00 20.72      P9
ATOM    106  CA  VAL A 109     122.621  24.018  18.411  1.00 23.20      P9
ATOM    107  CA  GLU A 110     120.751  26.260  15.961  1.00 28.03      P9
ATOM    108  CA  GLU A 111     121.704  25.099  12.460  0.74 31.85      P9
ATOM    109  CA  GLU A 112     122.954  28.554  11.512  1.00 33.44      P9
ATOM    110  CA  ALA A 113     125.368  28.615  14.458  1.00 31.90      P9
ATOM    111  CA  LYS A 114     126.911  25.180  13.906  0.92 31.27      P9
ATOM    112  CA  GLY A 115     129.136  26.579  11.174  1.00 30.17      P9
ATOM    113  CA  ARG A 116     130.956  29.004  13.463  1.00 27.89      P9
ATOM    114  CA  LEU A 117     131.365  26.777  16.514  1.00 25.28      P9
ATOM    115  CA  ALA A 118     134.998  27.063  17.578  1.00 23.01      P9
ATOM    116  CA  PRO A 119     136.912  26.997  20.869  1.00 21.21      P9
ATOM    117  CA  GLY A 120     137.158  30.500  22.308  1.00 21.69      P9
ATOM    118  CA  ALA A 121     134.046  31.875  20.620  1.00 21.58      P9
ATOM    119  CA  GLU A 122     131.240  33.298  22.745  1.00 22.41      P9
ATOM    120  CA  VAL A 123     127.678  32.132  22.205  1.00 20.97      P9
ATOM    121  CA  GLU A 124     124.171  33.150  23.219  1.00 20.32      P9
ATOM    122  CA  VAL A 125     122.627  30.323  25.229  1.00 19.82      P9
ATOM    123  CA  TRP A 126     119.012  29.745  26.239  1.00 21.08      P9
ATOM    124  CA  GLN A 127     118.027  27.643  29.216  1.00 22.22      P9
ATOM    125  CA  ILE A 128     114.538  26.200  29.585  1.00 22.80      P9
ATOM    126  CA  LEU A 129     114.153  24.029  32.671  1.00 24.59      P9
ATOM    127  CA  ASP A 130     117.058  21.551  32.616  0.99 24.90      P9
ATOM    128  CA  ARG A 131     117.952  21.918  28.938  1.00 22.57      P9
ATOM    129  CA  TYR A 132     120.128  24.357  27.016  1.00 21.18      P9
ATOM    130  CA  LYS A 133     120.153  25.537  23.432  1.00 21.22      P9
ATOM    131  CA  ILE A 134     122.957  27.373  21.662  1.00 20.82      P9
ATOM    132  CA  ILE A 135     121.298  30.133  19.637  1.00 24.27      P9
ATOM    133  CA  ARG A 136     124.175  31.872  17.866  1.00 28.80      P9
ATOM    134  CA  VAL A 137     127.877  32.680  17.955  1.00 32.92      P9
ATOM    135  CA  LYS A 138     128.895  35.931  19.718  1.00 36.89      P9
ATOM    136  CA  GLY A 139     126.242  36.878  22.264  0.67 38.21      P9
ATOM      1  CA  LYS B   4     514.976 423.221 449.266  1.00 48.61      P9
ATOM      2  CA  TRP B   5     514.182 426.755 448.110  1.00 47.90      P9
ATOM      3  CA  VAL B   6     515.823 428.654 445.250  1.00 46.89      P9
ATOM      4  CA  MET B   7     513.936 431.802 444.218  1.00 45.06      P9
ATOM      5  CA  SER B   8     510.419 433.211 444.036  1.00 41.81      P9
ATOM      6  CA  THR B   9     508.906 432.616 440.600  1.00 37.71      P9
ATOM      7  CA  LYS B  10     508.277 435.858 438.696  1.00 33.20      P9
ATOM      8  CA  TYR B  11     505.215 436.269 436.467  1.00 29.03      P9
ATOM      9  CA  VAL B  12     503.628 438.636 433.957  1.00 24.81      P9
ATOM     10  CA  GLU B  13     500.372 438.451 432.024  1.00 21.62      P9
ATOM     11  CA  ALA B  14     500.648 436.766 428.627  1.00 21.34      P9
ATOM     12  CA  GLY B  15     499.222 439.881 427.004  1.00 22.33      P9
ATOM     13  CA  GLU B  16     502.221 441.890 428.214  1.00 23.17      P9
ATOM     14  CA  LEU B  17     504.807 439.753 426.414  1.00 23.42      P9
ATOM     15  CA  LYS B  18     506.429 441.062 423.237  1.00 24.53      P9
ATOM     16  CA  GLU B  19     508.948 439.832 420.691  1.00 24.34      P9
ATOM     17  CA  GLY B  20     512.227 439.699 422.578  1.00 22.47      P9
ATOM     18  CA  SER B  21     510.632 438.977 425.956  1.00 20.20      P9
ATOM     19  CA  TYR B  22     511.668 435.890 427.912  1.00 19.09      P9
ATOM     20  CA  VAL B  23     509.283 433.225 429.176  1.00 18.63      P9
ATOM     21  CA  VAL B  24     509.403 429.711 430.602  1.00 16.73      P9
ATOM     22  CA  ILE B  25     507.485 427.073 428.645  1.00 16.21      P9
ATOM     23  CA  ASP B  26     507.424 423.531 430.042  1.00 17.02      P9
ATOM     24  CA  GLY B  27     510.535 424.131 432.135  1.00 19.18      P9
ATOM     25  CA  GLU B  28     512.647 425.739 429.398  1.00 18.42      P9
ATOM     26  CA  PRO B  29     513.490 429.457 429.454  1.00 17.13      P9
ATOM     27  CA  CYS B  30     512.738 430.793 425.958  1.00 18.27      P9
ATOM     28  CA  ARG B  31     513.013 434.010 423.950  1.00 18.81      P9
ATOM     29  CA  VAL B  32     509.720 435.074 422.348  1.00 18.68      P9
ATOM     30  CA  VAL B  33     509.990 435.324 418.554  1.00 21.62      P9
ATOM     31  CA  GLU B  34     506.337 435.903 417.676  1.00 22.97      P9
ATOM     32  CA  ILE B  35     502.990 436.567 419.325  1.00 23.13      P9
ATOM     33  CA  GLU B  36     499.516 436.390 417.807  1.00 22.64      P9
ATOM     34  CA  LYS B  37     496.314 437.656 419.422  1.00 21.25      P9
ATOM     35  CA  SER B  38     492.754 436.489 418.758  1.00 20.23      P9
ATOM     36  CA  LYS B  39     489.259 436.708 420.215  1.00 22.33      P9
ATOM     37  CA  THR B  40     486.571 434.194 419.236  1.00 25.52      P9
ATOM     38  CA  GLY B  41     484.184 434.756 422.114  0.93 29.60      P9
ATOM     39  CA  LYS B  42     481.787 437.652 421.597  0.89 32.20      P9
ATOM     40  CA  HIS B  43     482.170 438.509 425.288  1.00 31.87      P9
ATOM     41  CA  GLY B  44     485.141 436.264 425.999  0.65 28.63      P9
ATOM     42  CA  SER B  45     488.784 436.911 426.813  1.00 23.69      P9
ATOM     43  CA  ALA B  46     491.373 437.582 424.132  1.00 18.89      P9
ATOM     44  CA  LYS B  47     493.902 434.776 423.670  1.00 17.38      P9
ATOM     45  CA  ALA B  48     497.611 434.914 422.899  1.00 17.02      P9
ATOM     46  CA  ARG B  49     499.612 432.386 420.893  1.00 17.11      P9
ATOM     47  CA  ILE B  50     503.284 432.635 421.768  1.00 17.68      P9
ATOM     48  CA  VAL B  51     506.124 431.174 419.718  1.00 17.34      P9
ATOM     49  CA  ALA B  52     509.518 431.094 421.413  1.00 17.71      P9
ATOM     50  CA  VAL B  53     512.935 429.470 421.254  1.00 17.54      P9
ATOM     51  CA  GLY B  54     514.704 427.649 424.055  1.00 18.35      P9
ATOM     52  CA  VAL B  55     517.688 429.550 425.407  1.00 19.01      P9
ATOM     53  CA  PHE B  56     519.509 426.294 426.108  1.00 19.90      P9
ATOM     54  CA  ASP B  57     518.198 423.603 423.749  1.00 21.17      P9
ATOM     55  CA  GLY B  58     517.293 425.976 420.923  1.00 21.69      P9
ATOM     56  CA  GLY B  59     513.988 424.187 420.479  1.00 21.62      P9
ATOM     57  CA  LYS B  60     510.864 425.872 419.151  1.00 20.30      P9
ATOM     58  CA  ARG B  61     508.024 425.971 421.667  1.00 18.60      P9
ATOM     59  CA  THR B  62     504.517 427.398 421.607  1.00 17.44      P9
ATOM     60  CA  LEU B  63     501.998 428.440 424.242  1.00 16.19      P9
ATOM     61  CA  SER B  64     498.370 429.446 423.716  1.00 15.39      P9
ATOM     62  CA  LEU B  65     496.433 430.924 426.644  1.00 14.67      P9
ATOM     63  CA  PRO B  66     494.260 433.905 427.730  1.00 14.79      P9
ATOM     64  CA  VAL B  67     495.940 437.325 427.634  1.00 15.74      P9
ATOM     65  CA  ASP B  68     495.325 437.738 431.369  1.00 17.10      P9
ATOM     66  CA  ALA B  69     496.793 434.371 432.372  1.00 16.26      P9
ATOM     67  CA  GLN B  70     500.115 434.425 434.223  1.00 18.00      P9
ATOM     68  CA  VAL B  71     503.237 433.190 432.463  1.00 19.41      P9
ATOM     69  CA  GLU B  72     506.563 432.596 434.201  1.00 21.29      P9
ATOM     70  CA  VAL B  73     509.492 434.822 433.262  1.00 21.83      P9
ATOM     71  CA  PRO B  74     513.047 433.603 433.963  0.94 22.19      P9
ATOM     72  CA  ILE B  75     515.436 435.660 436.072  1.00 23.67      P9
ATOM     73  CA  ILE B  76     518.508 436.336 433.943  0.96 24.46      P9
ATOM     74  CA  GLU B  77     521.735 437.457 435.611  1.00 24.69      P9
ATOM     75  CA  LYS B  78     524.294 438.935 433.213  1.00 24.12      P9
ATOM     76  CA  PHE B  79     528.013 439.160 433.959  1.00 21.39      P9
ATOM     77  CA  THR B  80     531.479 439.196 432.404  1.00 19.30      P9
ATOM     78  CA  ALA B  81     533.958 436.320 432.436  1.00 17.03      P9
ATOM     79  CA  GLN B  82     537.421 435.668 430.993  1.00 15.83      P9
ATOM     80  CA  ILE B  83     538.128 432.536 428.964  1.00 14.35      P9
ATOM     81  CA  LEU B  84     540.770 430.349 430.639  1.00 14.01      P9
ATOM     82  CA  SER B  85     540.720 427.442 428.163  1.00 14.77      P9
ATOM     83  CA  VAL B  86     538.697 426.017 425.292  1.00 15.76      P9
ATOM     84  CA  SER B  87     538.246 422.344 424.410  1.00 15.97      P9
ATOM     85  CA  GLY B  88     535.905 420.871 421.836  1.00 16.65      P9
ATOM     86  CA  ASP B  89     533.352 420.243 424.587  1.00 17.36      P9
ATOM     87  CA  VAL B  90     533.894 422.746 427.389  1.00 16.73      P9
ATOM     88  CA  ILE B  91     534.780 426.401 427.870  1.00 16.05      P9
ATOM     89  CA  GLN B  92     536.584 427.059 431.178  1.00 16.32      P9
ATOM     90  CA  LEU B  93     535.870 430.590 432.465  1.00 15.69      P9
ATOM     91  CA  MET B  94     536.982 432.904 435.270  1.00 17.32      P9
ATOM     92  CA  ASP B  95     534.006 434.830 436.697  1.00 19.94      P9
ATOM     93  CA  MET B  96     535.333 438.412 436.817  1.00 22.68      P9
ATOM     94  CA  ARG B  97     533.065 439.204 439.767  1.00 24.59      P9
ATOM     95  CA  ASP B  98     534.544 436.755 442.283  0.96 24.84      P9
ATOM     96  CA  TYR B  99     537.285 435.135 440.175  0.99 24.06      P9
ATOM     97  CA  LYS B 100     535.841 431.625 440.550  1.00 22.14      P9
ATOM     98  CA  THR B 101     536.149 429.064 437.761  1.00 20.84      P9
ATOM     99  CA  ILE B 102     532.988 428.035 435.916  1.00 19.78      P9
ATOM    100  CA  GLU B 103     532.678 425.418 433.177  1.00 19.97      P9
ATOM    101  CA  VAL B 104     530.251 425.963 430.299  1.00 18.00      P9
ATOM    102  CA  PRO B 105     529.362 423.242 427.754  1.00 18.24      P9
ATOM    103  CA  MET B 106     530.166 424.138 424.128  1.00 18.91      P9
ATOM    104  CA  LYS B 107     526.491 423.558 423.303  1.00 19.94      P9
ATOM    105  CA  TYR B 108     525.616 426.558 425.496  1.00 20.72      P9
ATOM    106  CA  VAL B 109     527.621 429.018 423.411  1.00 23.20      P9
ATOM    107  CA  GLU B 110     525.751 431.260 420.961  1.00 28.03      P9
ATOM    108  CA  GLU B 111     526.704 430.099 417.460  0.74 31.85      P9
ATOM    109  CA  GLU B 112     527.954 433.554 416.512  1.00 33.44      P9
ATOM    110  CA  ALA B 113     530.368 433.615 419.458  1.00 31.90      P9
ATOM    111  CA  LYS B 114     531.911 430.180 418.906  0.92 31.27      P9
ATOM    112  CA  GLY B 115     534.136 431.579 416.174  1.00 30.17      P9
ATOM    113  CA  ARG B 116     535.956 434.004 418.463  1.00 27.89      P9
ATOM    114  CA  LEU B 117     536.365 431.777 421.514  1.00 25.28      P9
ATOM    115  CA  ALA B 118     539.998 432.063 422.578  1.00 23.01      P9
ATOM    116  CA  PRO B 119     541.912 431.997 425.869  1.00 21.21      P9
ATOM    117  CA  GLY B 120     542.158 435.500 427.308  1.00 21.69      P9
ATOM    118  CA  ALA B 121     539.046 436.875 425.620  1.00 21.58      P9
ATOM    119  CA  GLU B 122     536.240 438.298 427.745  1.00 22.41      P9
ATOM    120  CA  VAL B 123     532.678 437.132 427.205  1.00 20.97      P9
ATOM    121  CA  GLU B 124     529.171 438.150 428.219  1.00 20.32      P9
ATOM    122  CA  VAL B 125     527.627 435.323 430.229  1.00 19.82      P9
ATOM    123  CA  TRP B 126     524.012 434.745 431.239  1.00 21.08      P9
ATOM    124  CA  GLN B 127     523.027 432.643 434.216  1.00 22.22      P9
ATOM    125  CA  ILE B 128     519.538 431.200 434.585  1.00 22.80      P9
ATOM    126  CA  LEU B 129     519.153 429.029 437.671  1.00 24.59      P9
ATOM    127  CA  ASP B 130     522.058 426.551 437.616  0.99 24.90      P9
ATOM    128  CA  ARG B 131     522.952 426.918 433.938  1.00 22.57      P9
ATOM    129  CA  TYR B 132     525.128 429.357 432.016  1.00 21.18      P9
ATOM    130  CA  LYS B 133     525.153 430.537 428.432  1.00 21.22      P9
ATOM    131  CA  ILE B 134     527.957 432.373 426.662  1.00 20.82      P9
ATOM    132  CA  ILE B 135     526.298 435.133 424.637  1.00 24.27      P9
ATOM    133  CA  ARG B 136     529.175 436.872 422.866  1.00 28.80      P9
ATOM    134  CA  VAL B 137     532.877 437.680 422.955  1.00 32.92      P9
ATOM    135  CA  LYS B 138     533.895 440.931 424.718  1.00 36.89      P9
ATOM    136  CA  GLY B 139     531.242 441.878 427.264  0.67 38.21      P9
"""
def tst_01():
  print("Comparing mixed model with target...")
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
  expected_text=\
"""SEQ SCORE is fraction (close and matching target sequence).
MEAN LENGTH is the mean length of contiguous segments in the match with target sequence. (Each gap/reverse of direction starts new segment).



               ----ALL RESIDUES---  CLOSE RESIDUES ONLY    %
     MODEL     --CLOSE-    --FAR-- FORWARD REVERSE MIXED FOUND  CA                  SEQ
               RMSD   N      N       N       N      N          SCORE  SEQ MATCH(%)  SCORE  MEAN LENGTH  FRAGMENTS BAD CONNECTIONS

 Unique_target 1.55   54      7     14      29      11   39.7   0.26     9.3        0.04    6.0          7               6"""
  found_text="\n".join(f.getvalue().splitlines()[-10:])
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("Expected: \n%s \nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")
  from libtbx.test_utils import approx_equal
  print(r.get_values("forward"))
  assert approx_equal(r.get_values("forward"),(1.6751069901864204, 14))
  print(r.get_values("reverse"))
  assert approx_equal(r.get_values("reverse"),(1.388466550576198, 29))
  print(r.get_values("close"))
  assert approx_equal(r.get_values("close"),(1.545835235099158, 54))
  print(r.get_values("all_far"))
  assert approx_equal(r.get_values("all_far"),(0,0))
  print("OK")

def tst_02():
  print("Comparing mixed model with target with 2 chains...")
  import iotbx.pdb
  from cctbx.array_family import flex
  model_pdb_inp=iotbx.pdb.input(source_info='model',
       lines=flex.split_lines(model1))
  crystal_symmetry=model_pdb_inp.crystal_symmetry()
  model_hierarchy=model_pdb_inp.construct_hierarchy()
  query_hierarchy=iotbx.pdb.input(source_info='query',
       lines=flex.split_lines(query)).construct_hierarchy()

  f=StringIO()
  r=run(crystal_symmetry=crystal_symmetry,
    chain_hierarchy=query_hierarchy,target_hierarchy=model_hierarchy,out=f)
  expected_text="""
  SEQ SCORE is fraction (close and matching target sequence).
MEAN LENGTH is the mean length of contiguous segments in the match with target sequence. (Each gap/reverse of direction starts new segment).



               ----ALL RESIDUES---  CLOSE RESIDUES ONLY    %
     MODEL     --CLOSE-    --FAR-- FORWARD REVERSE MIXED FOUND  CA                  SEQ
               RMSD   N      N       N       N      N          SCORE  SEQ MATCH(%)  SCORE  MEAN LENGTH  FRAGMENTS BAD CONNECTIONS

 Unique_target 1.55   54      7     14      29      11   39.7   0.26     9.3        0.04    6.0          7               6"""
  found_text="\n".join(f.getvalue().splitlines()[-10:])
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("\n\nExpected: \n%s \n\nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")

  from libtbx.test_utils import approx_equal
  print(r.get_values("forward"))
  assert approx_equal(r.get_values("forward"),(1.6751069901864204, 14))
  print(r.get_values("reverse"))
  assert approx_equal(r.get_values("reverse"),(1.388466550576198, 29))
  print(r.get_values("close"))
  assert approx_equal(r.get_values("close"),(1.545835235099158, 54))
  print(r.get_values("all_far"))
  assert approx_equal(r.get_values("all_far"),(0,0))
  print("OK")

def tst_03():
  print("Comparing mixed model with target with 2 chains...as group")
  import iotbx.pdb
  from cctbx.array_family import flex
  model_pdb_inp=iotbx.pdb.input(source_info='model',
       lines=flex.split_lines(model1))
  crystal_symmetry=model_pdb_inp.crystal_symmetry()
  model_hierarchy=model_pdb_inp.construct_hierarchy()
  query_hierarchy=iotbx.pdb.input(source_info='query',
       lines=flex.split_lines(query)).construct_hierarchy()
  import os
  if not os.path.isdir("files"):
    os.mkdir("files")

  ff=open(os.path.join("files","query.pdb"),'w')
  print("CRYST1  113.949  113.949   32.474  90.00  90.00  90.00 I 4", file=ff)
  print(query_hierarchy.as_pdb_string(), file=ff)
  ff.close()
  ff=open("model.pdb",'w')
  print(model_hierarchy.as_pdb_string(), file=ff)
  ff.close()

  f=StringIO()
  args=["query_dir=files","model.pdb"]
  r=run(args,out=f)
  expected_text="""
  SEQ SCORE is fraction (close and matching target sequence).
MEAN LENGTH is the mean length of contiguous segments in the match with target sequence. (Each gap/reverse of direction starts new segment).



               ----ALL RESIDUES---  CLOSE RESIDUES ONLY    %
     MODEL     --CLOSE-    --FAR-- FORWARD REVERSE MIXED FOUND  CA                  SEQ
               RMSD   N      N       N       N      N          SCORE  SEQ MATCH(%)  SCORE  MEAN LENGTH  FRAGMENTS BAD CONNECTIONS

     query.pdb 1.55   54      7     14      29      11   39.7   0.26     9.3        0.04    6.0          7               6"""
  found_text="\n".join(f.getvalue().splitlines()[-10:])
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("Expected at tst_03: \n%s \nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")
  print("OK")

def tst_04():
    print("Testing choosing unique sequences")

    from mmtbx.validation.chain_comparison import \
       extract_unique_part_of_sequences as eups

    seqs=[ "abcdefgh","klafmalsd"]
    print()
    print(seqs)
    copies_in_unique,base_copies,unique_sequence_dict=eups(seqs)
    for seq in copies_in_unique.keys():
      print(copies_in_unique[seq],base_copies,seq)

    seqs=[ "abcdefgh",
           "klafmalsd",
           "klafmalsd"]
    print()
    print(seqs)
    copies_in_unique,base_copies,unique_sequence_dict=eups(seqs)
    for seq in copies_in_unique.keys():
      print(copies_in_unique[seq],base_copies,seq)

    seqs=[
        "abcdefgh",
        "klafmalsd",
         ]
    print()
    print(seqs)
    copies_in_unique,base_copies,unique_sequence_dict=eups(seqs)
    for seq in copies_in_unique.keys():
      print(copies_in_unique[seq],base_copies,seq)

    seqs=[
        "abcdefgh",
        "klafmalsd",
        "abcdefgh",
        "klafmalsd",
         ]
    print()
    print(seqs)
    copies_in_unique,base_copies,unique_sequence_dict=eups(seqs)
    for seq in copies_in_unique.keys():
      print(copies_in_unique[seq],base_copies,seq)

    seqs=[
        "abcdefgh",
        "klafmalsd",
        "klafmalsd",
        "abcdefgh",
        "klafmalsd",
        "klafmalsd",
         ]
    print()
    print(seqs)
    copies_in_unique,base_copies,unique_sequence_dict=eups(seqs)
    for seq in copies_in_unique.keys():
      print(copies_in_unique[seq],base_copies,seq)

    print("OK")

target="""
CRYST1  113.949  113.949   32.474  90.00  90.00  90.00 I 4
ATOM      9  CA  LYS U   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP U   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL U   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET U   7     108.936  26.802  39.218  1.00 45.06      P9
"""

modela="""
CRYST1  113.949  113.949   32.474  90.00  90.00  90.00 I 4
ATOM      9  CA  LYS A   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP A   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL A   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET A   7     108.936  26.802  39.218  1.00 45.06      P9
"""
modelb="""
ATOM     40  CA  SER B   8     205.419  28.211  39.036  1.00 41.81      P9
ATOM     46  CA  THR B   9     203.906  27.616  35.600  1.00 37.71      P9
ATOM     53  CA  LYS B  10     203.277  30.858  33.696  1.00 33.20      P9
ATOM     62  CA  TYR B  11     200.215  31.269  31.467  1.00 29.03      P9
"""

modelaa="""
ATOM      9  CA  LYS A   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP A   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL A   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET A   7     108.936  26.802  39.218  1.00 45.06      P9
ATOM      9  CA  LYS C  14     209.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP C  15     209.182  21.755  43.110  1.00 47.90      P9
ATOM     15  CA  VAL C  16     210.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET C  17     208.936  26.802  39.218  1.00 45.06      P9
"""

modelaab="""
ATOM      9  CA  LYS A   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP A   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL A   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET A   7     108.936  26.802  39.218  1.00 45.06      P9
ATOM      9  CA  LYS C  14     209.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP C  15     209.182  21.755  43.110  1.00 47.90      P9
ATOM     15  CA  VAL C  16     210.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET C  17     208.936  26.802  39.218  1.00 45.06      P9
ATOM     40  CA  SER B   8     205.419  28.211  39.036  1.00 41.81      P9
ATOM     46  CA  THR B   9     203.906  27.616  35.600  1.00 37.71      P9
ATOM     53  CA  LYS B  10     203.277  30.858  33.696  1.00 33.20      P9
ATOM     62  CA  TYR B  11     200.215  31.269  31.467  1.00 29.03      P9
"""

modelaabaab="""
ATOM      9  CA  LYS A   4     109.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP A   5     109.182  21.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL A   6     110.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET A   7     108.936  26.802  39.218  1.00 45.06      P9
ATOM      9  CA  LYS C  14     209.976  18.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP C  15     209.182  21.755  43.110  1.00 47.90      P9
ATOM     15  CA  VAL C  16     210.823  23.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET C  17     208.936  26.802  39.218  1.00 45.06      P9
ATOM     40  CA  SER B   8     205.419  28.211  39.036  1.00 41.81      P9
ATOM     46  CA  THR B   9     203.906  27.616  35.600  1.00 37.71      P9
ATOM     53  CA  LYS B  10     203.277  30.858  33.696  1.00 33.20      P9
ATOM     62  CA  TYR B  11     200.215  31.269  31.467  1.00 29.03      P9
ATOM      9  CA  LYS F   4     109.976 118.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP F   5     109.182 121.755  43.110  1.00 47.90      P9
ATOM     25  CA  VAL F   6     110.823 123.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET F   7     108.936 126.802  39.218  1.00 45.06      P9
ATOM      9  CA  LYS H  14     209.976 118.221  44.266  1.00 48.61      P9
ATOM     11  CA  TRP H  15     209.182 121.755  43.110  1.00 47.90      P9
ATOM     15  CA  VAL H  16     210.823 123.654  40.250  1.00 46.89      P9
ATOM     32  CA  MET H  17     208.936 126.802  39.218  1.00 45.06      P9
ATOM     40  CA  SER G   8     205.419 128.211  39.036  1.00 41.81      P9
ATOM     46  CA  THR G   9     203.906 127.616  35.600  1.00 37.71      P9
ATOM     53  CA  LYS G  10     203.277 130.858  33.696  1.00 33.20      P9
ATOM     62  CA  TYR G  11     200.215 131.269  31.467  1.00 29.03      P9
"""

ncs_spec="""
Summary of NCS information
Fri Jan 12 11:58:49 2018
/net/anaconda/raid1/terwill/misc/junk





new_ncs_group
new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth     0.0000    0.0000    0.0000

center_orth  112.2897   27.5224   23.1305
CHAIN A
RMSD 0
MATCHING 136
  RESSEQ 4:139

new_operator

rota_matrix    1.0000    0.0000    0.0000
rota_matrix    0.0000    1.0000    0.0000
rota_matrix    0.0000    0.0000    1.0000
tran_orth   -405.0000  -405.0000  -405.0000

center_orth  517.2897  432.5224  428.1305
CHAIN B
RMSD 4.32467304409e-13
MATCHING 136
  RESSEQ 4:139

"""
def tst_05():
  from mmtbx.validation.chain_comparison import \
       extract_unique_part_of_hierarchy as euph

  print("Testing extraction of unique part and unique matching")
  for m in [modela,modelb,modelaa,modelaab,modelaabaab]:
    import iotbx.pdb
    from cctbx.array_family import flex
    model_pdb_inp=iotbx.pdb.input(source_info='model',
         lines=flex.split_lines(m))
    crystal_symmetry=model_pdb_inp.crystal_symmetry()
    model_hierarchy=model_pdb_inp.construct_hierarchy()

    print("\nExtraction of unique MODEL with %s residues" %(
       model_hierarchy.overall_counts().n_residues))
    query_hierarchy=iotbx.pdb.input(source_info='query',
         lines=flex.split_lines(target)).construct_hierarchy()
    unique_hierarchy=euph(model_hierarchy,target_ph=query_hierarchy)
    print("FINAL chain ids: %s \n" %(" ".join(unique_hierarchy.chain_ids())))
  print("OK")

def tst_06():
  try:
    from phenix.programs import chain_comparison
  except Exception as e:
    return # no phenix

  print("Comparing mixed model with target with 2 chains...using ncs")
  import iotbx.pdb
  from cctbx.array_family import flex
  model_pdb_inp=iotbx.pdb.input(source_info='model',
       lines=flex.split_lines(model1))
  crystal_symmetry=model_pdb_inp.crystal_symmetry()
  f=open('ncs.ncs_spec','w')
  print(ncs_spec, file=f)
  f.close()
  f=open('model.pdb','w')
  print(model1, file=f) #model_hierarchy.as_pdb_string()
  f.close()
  f=open('query.pdb','w')
  print(query, file=f) #query_hierarchy.as_pdb_string()
  f.close()

  f=StringIO()
  args=['model.pdb','query.pdb','ncs_file=ncs.ncs_spec']
  run_program(program_class=chain_comparison.Program, args = args, logger = f)
  f.flush()
  expected_text="""
  SEQ SCORE is fraction (close and matching target sequence).
MEAN LENGTH is the mean length of contiguous segments in the match with target sequence. (Each gap/reverse of direction starts new segment).



               ----ALL RESIDUES---  CLOSE RESIDUES ONLY    %
     MODEL     --CLOSE-    --FAR-- FORWARD REVERSE MIXED FOUND  CA                  SEQ
               RMSD   N      N       N       N      N          SCORE  SEQ MATCH(%)  SCORE  MEAN LENGTH  FRAGMENTS BAD CONNECTIONS

 Unique_target 1.67   58     64     15      29      14   42.6   0.25     8.6        0.04    4.5         16               8"""
  found_text="\n".join(f.getvalue().strip().splitlines()[-17:-6])
  if remove_blank(found_text)!=remove_blank(expected_text):
    print("\n\nExpected tst_06: \n%s \n\nFound: \n%s" %(expected_text,found_text))
    raise AssertionError("FAILED")

  print("OK")

if __name__=="__main__":
  tst_01()
  tst_02()
  tst_03()
  tst_04()
  tst_05()
  tst_06()
