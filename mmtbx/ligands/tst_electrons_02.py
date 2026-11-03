from __future__ import absolute_import, division, print_function
'''
GeoStd ligands whose restraint files put -1 on a carbon with four bonds next to a
sulfonium S+ (DSK CAL, EU9 C7, KTL CAW, SSD C7; SSD also has its sulfate as in the
CCD, S9 with three O-). The formal charges cannot be built as given; without the
carbon's charge the file's other charges give the total: DSK +1, EU9 0, KTL 0,
SSD -2. AP5 (Ap5A, one O- on each of five phosphates; the model's H24 is five fewer
than the CCD's neutral H29) for contrast: the file's charges as given, -5, no carbon
charge ignored. Models from the Ligands_Bench_Charge set (3l4u, 6ltv, 3l4v, 3l4z,
5ycb), restraints from GeoStd: if those files are corrected the notes go and this
test needs updating.
'''
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.ligands import electrons
from mmtbx.ligands import rdkit_utils

pdbs = {}
pdbs['AP5'] = '''
CRYST1   36.977   43.089   38.640  90.00  90.00  90.00 P 1
SCALE1      0.027044  0.000000  0.000000        0.00000
SCALE2      0.000000  0.023208  0.000000        0.00000
SCALE3      0.000000  0.000000  0.025880        0.00000
HETATM    1  C1F AP5 A 201      22.555  15.208  19.826  1.00 39.69           C
HETATM    2  C1J AP5 A 201      15.472  29.187  10.985  1.00 39.82           C
HETATM    3  C2A AP5 A 201      21.874  14.252  24.238  1.00 26.18           C
HETATM    4  C2B AP5 A 201      11.878  31.726  10.757  1.00 53.70           C
HETATM    5  C2F AP5 A 201      21.354  14.764  19.513  1.00 38.76           C
HETATM    6  C2J AP5 A 201      16.757  29.219  11.047  1.00 45.37           C
HETATM    7  C3F AP5 A 201      20.737  15.853  18.600  1.00 37.05           C
HETATM    8  C3J AP5 A 201      17.207  27.816  10.587  1.00 47.32           C
HETATM    9  C4A AP5 A 201      22.276  15.593  22.326  1.00 25.71           C
HETATM   10  C4B AP5 A 201      13.734  30.717  11.807  1.00 50.48           C
HETATM   11  C4F AP5 A 201      21.946  16.224  17.778  1.00 33.23           C
HETATM   12  C4J AP5 A 201      16.126  26.903  11.168  1.00 40.36           C
HETATM   13  C5A AP5 A 201      22.283  16.703  23.095  1.00 31.92           C
HETATM   14  C5B AP5 A 201      13.401  31.337  12.973  1.00 52.94           C
HETATM   15  C5F AP5 A 201      21.715  17.660  17.258  1.00 32.25           C
HETATM   16  C5J AP5 A 201      16.556  26.381  12.579  1.00 38.34           C
HETATM   17  C6A AP5 A 201      22.070  16.576  24.458  1.00 30.14           C
HETATM   18  C6B AP5 A 201      12.285  32.172  13.032  1.00 51.90           C
HETATM   19  C8A AP5 A 201      22.637  17.322  21.037  1.00 32.00           C
HETATM   20  C8B AP5 A 201      15.211  30.137  13.322  1.00 49.10           C
HETATM   21  N1A AP5 A 201      21.874  15.365  25.016  1.00 28.86           N
HETATM   22  N1B AP5 A 201      11.544  32.353  11.928  1.00 54.73           N
HETATM   23  N3A AP5 A 201      22.078  14.373  22.895  1.00 30.97           N
HETATM   24  N3B AP5 A 201      12.973  30.910  10.700  1.00 53.71           N
HETATM   25  N6A AP5 A 201      22.022  17.577  25.490  1.00 30.16           N
HETATM   26  N6B AP5 A 201      11.731  32.908  14.061  1.00 41.74           N
HETATM   27  N7A AP5 A 201      22.506  17.785  22.294  1.00 31.14           N
HETATM   28  N7B AP5 A 201      14.323  30.975  13.913  1.00 47.14           N
HETATM   29  N9A AP5 A 201      22.495  15.969  21.046  1.00 33.81           N
HETATM   30  N9B AP5 A 201      14.860  29.974  12.019  1.00 44.60           N
HETATM   31  O1A AP5 A 201      22.690  20.636  17.432  1.00 36.12           O
HETATM   32  O1B AP5 A 201      18.443  20.447  18.794  1.00 49.92           O
HETATM   33  O1D AP5 A 201      16.759  24.738  17.421  1.00 53.59           O
HETATM   34  O1E AP5 A 201      18.496  28.434  15.177  1.00 39.30           O
HETATM   35  O1G AP5 A 201      20.271  23.182  15.573  1.00 43.38           O
HETATM   36  O2A AP5 A 201      21.271  20.765  19.498  1.00 34.42           O
HETATM   37  O2B AP5 A 201      18.119  19.460  16.551  1.00 40.61           O
HETATM   38  O2D AP5 A 201      19.034  25.543  17.042  1.00 29.67           O
HETATM   39  O2E AP5 A 201      19.437  26.446  14.082  1.00 35.33           O
HETATM   40  O2F AP5 A 201      21.409  13.489  18.743  1.00 37.14           O
HETATM   41  O2G AP5 A 201      18.608  21.693  14.406  1.00 52.01           O
HETATM   42  O2J AP5 A 201      17.424  30.238  10.176  1.00 52.43           O
HETATM   43  O3A AP5 A 201      20.365  20.339  17.102  1.00 61.00           O
HETATM   44  O3B AP5 A 201      18.423  21.997  16.864  1.00 53.52           O
HETATM   45  O3D AP5 A 201      17.249  26.279  15.473  1.00 53.27           O
HETATM   46  O3F AP5 A 201      19.768  15.407  17.753  1.00 30.23           O
HETATM   47  O3G AP5 A 201      17.949  23.929  15.363  1.00 58.60           O
HETATM   48  O3J AP5 A 201      17.220  27.749   9.235  1.00 51.95           O
HETATM   49  O4F AP5 A 201      22.945  16.198  18.624  1.00 36.14           O
HETATM   50  O4J AP5 A 201      15.059  27.638  11.247  1.00 37.04           O
HETATM   51  O5F AP5 A 201      21.545  18.470  18.427  1.00 40.38           O
HETATM   52  O5J AP5 A 201      17.304  27.415  13.208  1.00 46.96           O
HETATM   53  PA  AP5 A 201      21.474  20.086  18.164  1.00 36.63           P
HETATM   54  PB  AP5 A 201      18.834  20.549  17.330  1.00 45.69           P
HETATM   55  PD  AP5 A 201      17.781  25.093  16.368  1.00 36.67           P
HETATM   56  PE  AP5 A 201      18.166  27.134  14.487  1.00 40.14           P
HETATM   57  PG  AP5 A 201      18.854  22.670  15.526  1.00 42.78           P
HETATM   58  H1F AP5 A 201      23.291  14.413  19.944  1.00 39.69           H
HETATM   59  H1J AP5 A 201      15.106  29.570  10.032  1.00 39.82           H
HETATM   60  H2A AP5 A 201      21.713  13.281  24.684  1.00 26.18           H
HETATM   61  H2B AP5 A 201      11.270  31.881   9.878  1.00 53.70           H
HETATM   62  H2F AP5 A 201      20.767  14.735  20.431  1.00 38.76           H
HETATM   63  H2J AP5 A 201      16.987  29.408  12.095  1.00 45.37           H
HETATM   64  H3F AP5 A 201      20.335  16.674  19.194  1.00 37.05           H
HETATM   65  H3J AP5 A 201      18.148  27.579  11.084  1.00 47.32           H
HETATM   66  H4F AP5 A 201      22.192  15.536  16.969  1.00 33.23           H
HETATM   67  H4J AP5 A 201      15.887  26.071  10.506  1.00 40.36           H
HETATM   68  H8A AP5 A 201      22.824  17.905  20.147  1.00 32.00           H
HETATM   69  H8B AP5 A 201      16.064  29.663  13.785  1.00 49.10           H
HETATM   70 H51A AP5 A 201      20.845  17.664  16.602  1.00 32.25           H
HETATM   71 H51B AP5 A 201      17.129  25.464  12.444  1.00 38.34           H
HETATM   72 H52A AP5 A 201      22.568  17.958  16.648  1.00 32.25           H
HETATM   73 H52B AP5 A 201      15.657  26.110  13.132  1.00 38.34           H
HETATM   74 H61A AP5 A 201      21.852  17.310  26.460  1.00 30.16           H
HETATM   75 H61B AP5 A 201      10.894  33.468  13.900  1.00 41.74           H
HETATM   76 H62A AP5 A 201      22.157  18.561  25.260  1.00 30.16           H
HETATM   77 H62B AP5 A 201      12.155  32.896  14.989  1.00 41.74           H
HETATM   78 HO2A AP5 A 201      20.560  13.362  18.271  1.00 37.14           H
HETATM   79 HO2B AP5 A 201      18.395  30.131  10.259  1.00 52.43           H
HETATM   80 HO3A AP5 A 201      19.451  16.166  17.221  1.00 30.23           H
HETATM   81 HO3B AP5 A 201      17.904  27.092   8.987  1.00 51.95           H
'''

pdbs['DSK'] = '''
CRYST1   31.038   30.965   29.294  90.00  90.00  90.00 P 1
SCALE1      0.032219  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032295  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034137        0.00000
HETATM    1  CAJ DSK A4001      13.572   5.562  12.786  1.00 25.34           C
HETATM    2  CAK DSK A4001      14.373  15.403  17.938  1.00 19.12           C
HETATM    3  CAL DSK A4001      15.426  12.386  15.680  1.00 19.18           C
HETATM    4  CAM DSK A4001      16.472  14.966  14.520  1.00 22.15           C
HETATM    5  CAN DSK A4001      13.836   7.038  12.462  1.00 25.00           C
HETATM    6  CAO DSK A4001      14.869  11.561  14.478  1.00 21.53           C
HETATM    7  CAP DSK A4001      17.508  15.213  15.636  1.00 21.80           C
HETATM    8  CAQ DSK A4001      14.257   7.786  13.746  1.00 21.54           C
HETATM    9  CAR DSK A4001      14.911  10.054  14.717  1.00 22.14           C
HETATM   10  CAS DSK A4001      14.361   9.288  13.491  1.00 21.62           C
HETATM   11  CAT DSK A4001      16.672  15.881  16.762  1.00 20.06           C
HETATM   12  CAU DSK A4001      15.564  14.846  17.110  1.00 20.32           C
HETATM   13  OAA DSK A4001      12.831   5.000  11.716  1.00 30.82           O
HETATM   14  OAB DSK A4001      13.850  16.639  17.441  1.00 20.68           O
HETATM   15  OAC DSK A4001      14.888   7.149  11.492  1.00 24.89           O
HETATM   16  OAD DSK A4001      15.573  11.860  13.281  1.00 22.81           O
HETATM   17  OAE DSK A4001      18.572  16.123  15.222  1.00 20.83           O
HETATM   18  OAF DSK A4001      15.551   7.314  14.138  1.00 24.40           O
HETATM   19  OAG DSK A4001      13.990   9.804  15.793  1.00 22.00           O
HETATM   20  OAH DSK A4001      13.072   9.786  13.136  1.00 22.77           O
HETATM   21  OAI DSK A4001      17.552  16.066  17.878  1.00 18.81           O
HETATM   22  SAV DSK A4001      15.119  14.158  15.469  1.00 21.46           S
HETATM   23  HAJ DSK A4001      13.037   5.503  13.734  1.00 25.34           H
HETATM   24  HAK DSK A4001      13.607  14.627  17.941  1.00 19.12           H
HETATM   25  HAL DSK A4001      16.507  12.294  15.788  1.00 19.18           H
HETATM   26  HAM DSK A4001      16.813  14.282  13.743  1.00 22.15           H
HETATM   27  HAN DSK A4001      12.925   7.487  12.065  1.00 25.00           H
HETATM   28  HAO DSK A4001      13.835  11.889  14.374  1.00 21.53           H
HETATM   29  HAP DSK A4001      17.937  14.268  15.968  1.00 21.80           H
HETATM   30  HAQ DSK A4001      13.538   7.591  14.542  1.00 21.54           H
HETATM   31  HAR DSK A4001      15.932   9.751  14.949  1.00 22.14           H
HETATM   32  HAS DSK A4001      15.046   9.471  12.663  1.00 21.62           H
HETATM   33  HAT DSK A4001      16.180  16.815  16.491  1.00 20.06           H
HETATM   34  HAU DSK A4001      15.905  13.966  17.655  1.00 20.32           H
HETATM   35 HAJA DSK A4001      14.528   5.060  12.934  1.00 25.34           H
HETATM   36 HAKA DSK A4001      14.731  15.514  18.962  1.00 19.12           H
HETATM   37 HALA DSK A4001      14.951  12.124  16.625  1.00 19.18           H
HETATM   38 HAMA DSK A4001      16.086  15.878  14.065  1.00 22.15           H
HETATM   39 HOAA DSK A4001      11.939   5.405  11.719  1.00 30.82           H
HETATM   40 HOAB DSK A4001      13.649  16.527  16.488  1.00 20.68           H
HETATM   41 HOAC DSK A4001      15.364   7.996  11.621  1.00 24.89           H
HETATM   42 HOAD DSK A4001      15.688  12.830  13.200  1.00 22.81           H
HETATM   43 HOAE DSK A4001      19.110  15.697  14.522  1.00 20.83           H
HETATM   44 HOAF DSK A4001      15.792   7.737  14.989  1.00 24.40           H
HETATM   45 HOAG DSK A4001      14.281  10.300  16.587  1.00 22.00           H
HETATM   46 HOAH DSK A4001      13.066  10.763  13.206  1.00 22.77           H
HETATM   47 HOAI DSK A4001      17.944  15.202  18.123  1.00 18.81           H
'''

pdbs['EU9'] = '''
CRYST1   26.461   29.778   34.340  90.00  90.00  90.00 P 1
SCALE1      0.037791  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033582  0.000000        0.00000
SCALE3      0.000000  0.000000  0.029121        0.00000
HETATM    1  C1  EU9 A 503      11.117  16.309  16.456  1.00 32.42           C
HETATM    2  C12 EU9 A 503      12.854  17.112  22.681  1.00 24.30           C
HETATM    3  C13 EU9 A 503      13.221  18.585  22.356  1.00 24.76           C
HETATM    4  C14 EU9 A 503      14.383  19.126  23.161  1.00 29.93           C
HETATM    5  C15 EU9 A 503      15.113  20.261  22.440  1.00 30.36           C
HETATM    6  C19 EU9 A 503      11.968  14.846  21.960  1.00 23.95           C
HETATM    7  C2  EU9 A 503      11.455  13.426  14.739  1.00 25.97           C
HETATM    8  C20 EU9 A 503      11.970  13.565  21.143  1.00 20.71           C
HETATM    9  C21 EU9 A 503      11.106  16.388  17.856  1.00 29.75           C
HETATM   10  C22 EU9 A 503      13.034  11.862  20.004  1.00 20.20           C
HETATM   11  C23 EU9 A 503      12.224  12.709  19.031  1.00 20.59           C
HETATM   12  C27 EU9 A 503      13.271  12.828  21.168  1.00 19.56           C
HETATM   13  C3  EU9 A 503      12.263  16.779  18.582  1.00 33.42           C
HETATM   14  C4  EU9 A 503      12.724  13.332  16.694  1.00 22.51           C
HETATM   15  C41 EU9 A 503      13.445  17.133  17.864  1.00 35.01           C
HETATM   16  C5  EU9 A 503      13.891  13.778  15.926  1.00 25.28           C
HETATM   17  C51 EU9 A 503      13.447  17.005  16.466  1.00 31.34           C
HETATM   18  C6  EU9 A 503      13.714  14.001  14.565  1.00 25.43           C
HETATM   19  C61 EU9 A 503      12.291  16.634  15.760  1.00 31.67           C
HETATM   20  C7  EU9 A 503      12.097  16.894  20.099  1.00 27.66           C
HETATM   21  C8  EU9 A 503      14.439  13.494  18.027  1.00 23.75           C
HETATM   22  N1  EU9 A 503      12.513  13.794  14.015  1.00 22.83           N
HETATM   23  N18 EU9 A 503      15.272  18.254  23.886  1.00 25.02           N
HETATM   24  N3  EU9 A 503      11.535  13.193  16.074  1.00 21.86           N
HETATM   25  N6  EU9 A 503      14.766  14.373  13.779  1.00 29.19           N
HETATM   26  N7  EU9 A 503      14.921  13.825  16.793  1.00 25.09           N
HETATM   27  N8  EU9 A 503      14.648  17.504  18.505  1.00 43.75           N
HETATM   28  N9  EU9 A 503      13.127  13.185  17.946  1.00 23.16           N
HETATM   29  O10 EU9 A 503      14.660  17.975  19.804  1.00 47.86           O
HETATM   30  O16 EU9 A 503      16.213  20.544  22.849  1.00 27.88           O
HETATM   31  O17 EU9 A 503      14.482  20.957  21.631  1.00 30.93           O
HETATM   32  O24 EU9 A 503      11.705  13.813  19.795  1.00 19.52           O
HETATM   33  O25 EU9 A 503      13.502  12.082  22.345  1.00 20.55           O
HETATM   34  O26 EU9 A 503      12.269  10.718  20.379  1.00 19.54           O
HETATM   35  O9  EU9 A 503      15.883  17.368  17.855  1.00 58.14           O
HETATM   36  S11 EU9 A 503      12.948  16.120  21.333  1.00 23.91           S
HETATM   37  H2  EU9 A 503      13.713  19.439  23.962  1.00 29.93           H
HETATM   38  H3  EU9 A 503      16.177  18.750  24.013  1.00 25.02           H
HETATM   39  H10 EU9 A 503      12.243  17.930  20.406  1.00 27.66           H
HETATM   40  H11 EU9 A 503      11.082  16.610  20.379  1.00 27.66           H
HETATM   41  H12 EU9 A 503      10.203  16.148  18.398  1.00 29.75           H
HETATM   42  H13 EU9 A 503      10.229  16.001  15.924  1.00 32.42           H
HETATM   43  H14 EU9 A 503      12.316  16.602  14.681  1.00 31.67           H
HETATM   44  H15 EU9 A 503      14.357  17.195  15.916  1.00 31.34           H
HETATM   45  H16 EU9 A 503      10.957  15.248  22.031  1.00 23.95           H
HETATM   46  H17 EU9 A 503      12.332  14.652  22.969  1.00 23.95           H
HETATM   47  H18 EU9 A 503      11.166  12.936  21.525  1.00 20.71           H
HETATM   48  H19 EU9 A 503      11.377  12.196  18.575  1.00 20.59           H
HETATM   49  H20 EU9 A 503      13.966  11.547  19.535  1.00 20.20           H
HETATM   50  H21 EU9 A 503      11.369  11.013  20.630  1.00 19.54           H
HETATM   51  H22 EU9 A 503      14.089  13.516  20.956  1.00 19.56           H
HETATM   52  H23 EU9 A 503      13.887  12.665  23.033  1.00 20.55           H
HETATM   53  H24 EU9 A 503      15.009  13.477  18.944  1.00 23.75           H
HETATM   54  H25 EU9 A 503      14.632  14.536  12.781  1.00 29.19           H
HETATM   55  H26 EU9 A 503      15.695  14.492  14.183  1.00 29.19           H
HETATM   56  H27 EU9 A 503      10.504  13.311  14.241  1.00 25.97           H
HETATM   57  H4  EU9 A 503      14.764  17.368  24.083  1.00 25.02           H
HETATM   58  H6  EU9 A 503      13.447  18.656  21.292  1.00 24.76           H
HETATM   59  H7  EU9 A 503      12.339  19.204  22.523  1.00 24.76           H
HETATM   60  H8  EU9 A 503      13.525  16.681  23.424  1.00 24.30           H
HETATM   61  H9  EU9 A 503      11.835  17.021  23.056  1.00 24.30           H
'''

pdbs['KTL'] = '''
CRYST1   30.975   30.826   29.259  90.00  90.00  90.00 P 1
SCALE1      0.032284  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032440  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034178        0.00000
HETATM    1  CAL KTL A1001      13.850   5.379  12.809  1.00 32.77           C
HETATM    2  CAM KTL A1001      14.459  15.186  17.852  1.00 21.88           C
HETATM    3  CAN KTL A1001      16.501  14.594  14.560  1.00 23.48           C
HETATM    4  CAO KTL A1001      15.350  12.411  15.600  1.00 24.34           C
HETATM    5  CAQ KTL A1001      14.141   6.860  12.495  1.00 29.37           C
HETATM    6  CAR KTL A1001      17.501  15.128  15.605  1.00 22.77           C
HETATM    7  CAS KTL A1001      14.899  11.490  14.449  1.00 25.78           C
HETATM    8  CAT KTL A1001      14.596   7.614  13.746  1.00 27.37           C
HETATM    9  CAU KTL A1001      16.676  15.705  16.747  1.00 22.52           C
HETATM   10  CAV KTL A1001      14.531   9.142  13.567  1.00 25.81           C
HETATM   11  CAW KTL A1001      15.584  14.654  16.972  1.00 22.73           C
HETATM   12  CAX KTL A1001      15.117   9.987  14.713  1.00 26.12           C
HETATM   13  OAA KTL A1001      12.700   5.000  12.018  1.00 36.11           O
HETATM   14  OAB KTL A1001      13.849  16.399  17.379  1.00 23.87           O
HETATM   15  OAC KTL A1001      15.113   6.977  11.449  1.00 28.82           O
HETATM   16  OAD KTL A1001      18.485  16.104  15.176  1.00 25.95           O
HETATM   17  OAE KTL A1001      15.525  11.825  13.208  1.00 24.30           O
HETATM   18  OAF KTL A1001      15.918   7.177  14.030  1.00 24.99           O
HETATM   19  OAG KTL A1001      17.473  15.943  17.908  1.00 20.52           O
HETATM   20  OAH KTL A1001      13.158   9.526  13.385  1.00 22.75           O
HETATM   21  OAI KTL A1001      16.035   7.955  16.926  1.00 31.15           O
HETATM   22  OAJ KTL A1001      16.001  10.139  17.822  1.00 29.78           O
HETATM   23  OAK KTL A1001      14.030   8.695  18.107  1.00 29.74           O
HETATM   24  OAP KTL A1001      14.387   9.629  15.875  1.00 29.80           O
HETATM   25  SAY KTL A1001      15.160  14.084  15.439  1.00 24.99           S
HETATM   26  SAZ KTL A1001      15.136   9.095  17.217  1.00 29.17           S
HETATM   27  H22 KTL A1001      16.871  13.728  14.012  1.00 23.48           H
HETATM   28  H23 KTL A1001      16.418  12.282  15.778  1.00 24.34           H
HETATM   29  H24 KTL A1001      15.959  13.720  17.391  1.00 22.73           H
HETATM   30  HAL KTL A1001      13.650   5.263  13.874  1.00 32.77           H
HETATM   31  HAM KTL A1001      13.724  14.384  17.920  1.00 21.88           H
HETATM   32  HAN KTL A1001      16.154  15.358  13.864  1.00 23.48           H
HETATM   33  HAO KTL A1001      14.795  12.168  16.506  1.00 24.34           H
HETATM   34  HAQ KTL A1001      13.196   7.274  12.144  1.00 29.37           H
HETATM   35  HAR KTL A1001      18.047  14.219  15.857  1.00 22.77           H
HETATM   36  HAS KTL A1001      13.833  11.705  14.371  1.00 25.78           H
HETATM   37  HAT KTL A1001      13.913   7.385  14.564  1.00 27.37           H
HETATM   38  HAU KTL A1001      16.253  16.664  16.448  1.00 22.52           H
HETATM   39  HAV KTL A1001      15.155   9.382  12.706  1.00 25.81           H
HETATM   40  HAX KTL A1001      16.189   9.826  14.827  1.00 26.12           H
HETATM   41 HALA KTL A1001      14.719   4.773  12.554  1.00 32.77           H
HETATM   42 HAMA KTL A1001      14.893  15.331  18.841  1.00 21.88           H
HETATM   43 HOAA KTL A1001      11.925   5.501  12.347  1.00 36.11           H
HETATM   44 HOAB KTL A1001      13.609  16.284  16.436  1.00 23.87           H
HETATM   45 HOAC KTL A1001      16.005   7.044  11.849  1.00 28.82           H
HETATM   46 HOAD KTL A1001      19.005  15.706  14.447  1.00 25.95           H
HETATM   47 HOAE KTL A1001      16.348  11.303  13.109  1.00 24.30           H
HETATM   48 HOAF KTL A1001      16.517   7.795  13.562  1.00 24.99           H
HETATM   49 HOAG KTL A1001      17.507  15.127  18.450  1.00 20.52           H
HETATM   50 HOAH KTL A1001      12.992  10.365  13.864  1.00 22.75           H
'''

pdbs['SSD'] = '''
CRYST1   31.000   29.415   29.144  90.00  90.00  90.00 P 1
SCALE1      0.032258  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033996  0.000000        0.00000
SCALE3      0.000000  0.000000  0.034312        0.00000
HETATM    1  C1  SSD A1001      16.759  14.353  16.620  1.00 35.40           C
HETATM    2  C10 SSD A1001      14.742   7.700  13.589  1.00 42.72           C
HETATM    3  C2  SSD A1001      17.605  13.688  15.539  1.00 36.32           C
HETATM    4  C3  SSD A1001      16.608  13.435  14.389  1.00 37.91           C
HETATM    5  C5  SSD A1001      15.628  13.334  16.936  1.00 35.16           C
HETATM    6  C6  SSD A1001      14.489  13.972  17.754  1.00 33.93           C
HETATM    7  C7  SSD A1001      15.494  10.987  15.447  1.00 40.05           C
HETATM    8  C8  SSD A1001      14.914  10.089  14.362  1.00 39.81           C
HETATM    9  C9  SSD A1001      15.223   8.625  14.706  1.00 42.04           C
HETATM   10  O1  SSD A1001      17.585  14.545  17.768  1.00 33.65           O
HETATM   11  O10 SSD A1001      13.442   8.164  13.251  1.00 43.53           O
HETATM   12  O11 SSD A1001      16.184   8.571  17.710  1.00 43.22           O
HETATM   13  O12 SSD A1001      14.228   7.211  18.068  1.00 44.46           O
HETATM   14  O13 SSD A1001      16.013   6.460  16.506  1.00 45.97           O
HETATM   15  O2  SSD A1001      18.622  14.654  15.166  1.00 34.50           O
HETATM   16  O6  SSD A1001      14.014  15.154  17.129  1.00 35.63           O
HETATM   17  O8  SSD A1001      15.479  10.403  13.098  1.00 39.78           O
HETATM   18  O9  SSD A1001      14.535   8.270  15.919  1.00 42.78           O
HETATM   19  S4  SSD A1001      15.157  12.771  15.286  1.00 40.56           S
HETATM   20  S9  SSD A1001      15.242   7.638  17.021  1.00 43.36           S
HETATM   21  H1  SSD A1001      16.440  15.354  16.329  1.00 35.40           H
HETATM   22  H2  SSD A1001      18.165  12.822  15.892  1.00 36.32           H
HETATM   23  H10 SSD A1001      13.160   8.872  13.867  1.00 43.53           H
HETATM   24  H31 SSD A1001      16.985  12.719  13.658  1.00 37.91           H
HETATM   25  H32 SSD A1001      16.356  14.346  13.846  1.00 37.91           H
HETATM   26  H5  SSD A1001      15.891  12.435  17.493  1.00 35.16           H
HETATM   27  H61 SSD A1001      13.601  13.342  17.815  1.00 33.93           H
HETATM   28  H62 SSD A1001      14.809  14.315  18.738  1.00 33.93           H
HETATM   29  H71 SSD A1001      16.575  10.846  15.468  1.00 40.05           H
HETATM   30  H72 SSD A1001      15.108  10.655  16.411  1.00 40.05           H
HETATM   31  H8  SSD A1001      13.869  10.323  14.160  1.00 39.81           H
HETATM   32  H9  SSD A1001      16.290   8.455  14.848  1.00 42.04           H
HETATM   33  HO1 SSD A1001      17.802  13.686  18.188  1.00 33.65           H
HETATM   34  HO2 SSD A1001      19.307  14.233  14.605  1.00 34.50           H
HETATM   35  HO6 SSD A1001      13.884  15.010  16.168  1.00 35.63           H
HETATM   36  HO8 SSD A1001      15.390  11.368  12.950  1.00 39.78           H
HETATM   37 H101 SSD A1001      15.352   7.775  12.689  1.00 42.72           H
HETATM   38 H102 SSD A1001      14.638   6.665  13.916  1.00 42.72           H
'''

# resname: (total charge, carbon whose -1 in the restraint file is ignored, if any)
expected = {
  'AP5': (-5, None),
  'DSK': (1, 'CAL'),
  'EU9': (0, 'C7'),
  'KTL': (0, 'CAW'),
  'SSD': (-2, 'C7'),
}

def exercise_residue_molecule(code, pdb_str):
  total, carbon = expected[code]
  inp = iotbx.pdb.input(lines=pdb_str.split('\n'), source_info=None)
  model = mmtbx.model.manager(model_input=inp, log=null_out())
  model.process(make_restraints=True)
  rg = model.get_hierarchy().only_residue_group()
  r = rdkit_utils.residue_rigid_components(model=model, residue_group=rg).molecule
  assert r.ok, r.reason
  assert (r.total_charge, r.charge_certain) == (total, True), (code, r.total_charge,
    r.charge_certain)
  assert r.total_charge_source in ('restraint file', 'formal charges'), \
    r.total_charge_source
  ignored = [n for n in r.charge_notes if n.startswith('formal charge ignored on ')]
  if carbon is None:
    assert ignored == [] and r.differences['charges'] == [], (code, r.charge_notes,
      r.differences)
  else:
    assert ignored == ['formal charge ignored on %s (-1) (a carbon with four bonds): '
      'total %d, not %d' % (carbon, total, total - 1)], (code, r.charge_notes)

def exercise_program(code, pdb_str):
  model_filename = 'tst_electrons_02_%s.pdb' % code
  with open(model_filename, 'w') as f:
    f.write(pdb_str)
  rc = run_program(program_class=electrons.Program,
                   args=[model_filename],
                   logger=null_out(),
                  )
  assert rc.total_charge == expected[code][0], (code, rc.total_charge)

def main():
  for code in sorted(pdbs):
    exercise_residue_molecule(code, pdbs[code])
    exercise_program(code, pdbs[code])
  print('OK')

if __name__ == '__main__':
  main()
