from __future__ import absolute_import, division, print_function
import os
from libtbx import easy_run
from six.moves import range

pdbs = {"linking_test_CYS_CYS_alt_loc.pdb" : """
ATOM    274  N   CYS A  27      17.541   4.439  12.897  1.00 13.99           N
ATOM    275  CA  CYS A  27      16.566   5.527  12.862  1.00 14.57           C
ATOM    276  C   CYS A  27      16.236   6.026  11.467  1.00 14.53           C
ATOM    277  O   CYS A  27      15.254   6.760  11.351  1.00 16.95           O
ATOM    278  CB  CYS A  27      17.114   6.698  13.662  1.00 15.77           C
ATOM    279  SG  CYS A  27      17.230   6.332  15.443  1.00 17.57           S
ATOM   1225  CB  CYS A 123      14.607   7.591  16.260  1.00 24.16           C
ATOM   1226  SG  CYS A 123      15.316   5.939  15.946  1.00 20.05           S
ATOM   1217  N  ACYS A 123      15.023   7.279  18.624  0.58 26.40           N
ATOM   1219  CA ACYS A 123      15.266   8.190  17.491  0.58 25.69           C
ATOM   1221  C  ACYS A 123      14.764   9.599  17.776  0.58 26.33           C
ATOM   1223  O  ACYS A 123      14.197  10.238  16.886  0.58 28.70           O
ATOM   1227  OXTACYS A 123      14.975  10.081  18.878  0.58 28.31           O
ATOM   1218  N  BCYS A 123      15.023   7.288  18.685  0.42 25.68           N
ATOM   1220  CA BCYS A 123      15.108   8.205  17.548  0.42 25.86           C
ATOM   1222  C  BCYS A 123      14.270   9.460  17.813  0.42 26.42           C
ATOM   1224  O  BCYS A 123      13.915  10.125  16.837  0.42 27.75           O
ATOM   1228  OXTBCYS A 123      13.981   9.728  18.968  0.42 28.04           O
""",
        "linking_test_NAG-NAG.pdb" : """
ATOM    982  N   ASN A 131      59.894  31.406  61.164  1.00  8.27           N
ATOM    983  CA  ASN A 131      60.727  32.209  60.265  1.00  8.51           C
ATOM    984  C   ASN A 131      62.218  31.862  60.359  1.00  9.14           C
ATOM    985  O   ASN A 131      63.057  32.700  59.954  1.00 11.54           O
ATOM    986  CB  ASN A 131      60.527  33.684  60.538  1.00  8.67           C
ATOM    987  CG  ASN A 131      60.938  34.568  59.363  1.00  9.10           C
ATOM    988  OD1 ASN A 131      60.572  34.301  58.221  1.00  8.69           O
ATOM    989  ND2 ASN A 131      61.660  35.631  59.674  1.00 10.39           N
HETATM 2693  C1  NAG A 361      61.950  36.682  58.740  1.00 10.88           C
HETATM 2694  C2  NAG A 361      61.484  38.019  59.269  1.00 11.23           C
HETATM 2695  C3  NAG A 361      61.917  39.166  58.345  1.00 12.36           C
HETATM 2696  C4  NAG A 361      63.416  39.059  58.041  1.00 12.10           C
HETATM 2697  C5  NAG A 361      63.760  37.632  57.548  1.00 12.32           C
HETATM 2698  C6  NAG A 361      65.221  37.430  57.213  1.00 16.14           C
HETATM 2699  C7  NAG A 361      59.438  37.677  60.568  1.00 11.18           C
HETATM 2700  C8  NAG A 361      57.935  37.538  60.529  1.00 12.29           C
HETATM 2701  N2  NAG A 361      60.056  37.983  59.429  1.00 11.34           N
HETATM 2702  O3  NAG A 361      61.602  40.407  58.980  1.00 14.07           O
HETATM 2703  O4  NAG A 361      63.792  39.965  56.971  1.00 13.40           O
HETATM 2704  O5  NAG A 361      63.357  36.678  58.496  1.00 11.54           O
HETATM 2705  O6  NAG A 361      65.973  37.639  58.364  1.00 17.57           O
HETATM 2706  O7  NAG A 361      60.051  37.500  61.614  1.00 12.67           O
HETATM 2707  C1  NAG A 362      64.268  41.214  57.404  1.00 15.39           C
HETATM 2708  C2  NAG A 362      65.186  41.796  56.314  1.00 16.15           C
HETATM 2709  C3  NAG A 362      65.630  43.192  56.758  1.00 19.28           C
HETATM 2710  C4  NAG A 362      64.431  44.087  57.131  1.00 19.23           C
HETATM 2711  C5  NAG A 362      63.526  43.362  58.126  1.00 18.32           C
HETATM 2712  C6  NAG A 362      62.219  44.114  58.366  1.00 18.44           C
HETATM 2713  C7  NAG A 362      66.649  40.270  55.063  1.00 15.73           C
HETATM 2714  C8  NAG A 362      67.859  39.361  55.043  1.00 17.34           C
HETATM 2715  N2  NAG A 362      66.303  40.872  56.196  1.00 16.94           N
HETATM 2716  O3  NAG A 362      66.321  43.771  55.701  1.00 23.91           O
HETATM 2717  O4  NAG A 362      64.893  45.314  57.706  1.00 22.28           O
HETATM 2718  O5  NAG A 362      63.154  42.081  57.584  1.00 16.19           O
HETATM 2719  O6  NAG A 362      61.496  43.532  59.419  1.00 20.18           O
HETATM 2720  O7  NAG A 362      66.009  40.389  54.022  1.00 14.38           O
""",
  "linking_test_MAN-SER.pdb" : """
HETATM 2721  C1  MAN A 364      43.704  27.933  15.645  1.00 16.07           C
HETATM 2722  C2  MAN A 364      42.296  27.563  15.194  1.00 20.45           C
HETATM 2723  C3  MAN A 364      42.280  26.343  14.240  1.00 21.29           C
HETATM 2724  C4  MAN A 364      43.262  26.552  13.084  1.00 19.50           C
HETATM 2725  C5  MAN A 364      44.616  27.026  13.597  1.00 15.55           C
HETATM 2726  C6  MAN A 364      45.515  27.538  12.507  1.00 19.27           C
HETATM 2727  O2  MAN A 364      41.697  28.638  14.482  1.00 21.51           O
HETATM 2728  O3  MAN A 364      40.973  26.218  13.718  1.00 26.03           O
HETATM 2729  O4  MAN A 364      43.352  25.365  12.310  1.00 23.56           O
HETATM 2730  O5  MAN A 364      44.525  28.138  14.481  1.00 15.77           O
HETATM 2731  O6  MAN A 364      46.856  27.797  12.809  1.00 32.01           O
ATOM   2544  N   SER A 336      46.708  25.446  16.017  1.00 10.08           N
ATOM   2545  CA  SER A 336      46.632  26.822  16.443  1.00 11.00           C
ATOM   2546  C   SER A 336      47.835  27.172  17.291  1.00 10.18           C
ATOM   2547  O   SER A 336      48.284  26.389  18.124  1.00 11.56           O
ATOM   2548  CB  SER A 336      45.335  27.061  17.256  1.00 12.96           C
ATOM   2549  OG  SER A 336      44.189  26.866  16.408  1.00 15.06           O
""",
    "linking_test_ASN-NAG.pdb" : """
ATOM   1989  N   ASN A 270      20.738  61.827 110.156  1.00 11.35           N
ATOM   1990  CA  ASN A 270      20.866  63.240 110.482  1.00 11.45           C
ATOM   1991  C   ASN A 270      20.238  64.075 109.377  1.00 11.25           C
ATOM   1992  O   ASN A 270      20.902  64.453 108.408  1.00 10.95           O
ATOM   1993  CB  ASN A 270      22.330  63.643 110.648  1.00 12.36           C
ATOM   1994  CG  ASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1 ASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2 ASN A 270      23.440  65.083 112.236  1.00 13.89           N
HETATM 3221  C1  NAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2  NAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3  NAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4  NAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5  NAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6  NAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7  NAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8  NAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2  NAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3  NAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4  NAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5  NAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6  NAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7  NAG A 435      22.047  63.959 114.195  1.00 19.58           O
""",
        "linking_test_CYS-VSP.pdb" : """
ATOM    631  N   CYS A 338     -12.642  13.597  13.517  1.00  7.89           N
ATOM    632  CA  CYS A 338     -12.754  14.691  12.554  1.00  6.21           C
ATOM    633  C   CYS A 338     -14.168  14.844  12.007  1.00  8.51           C
ATOM    634  O   CYS A 338     -15.020  13.986  12.195  1.00 13.75           O
ATOM    635  CB  CYS A 338     -11.756  14.486  11.426  1.00  5.09           C
ATOM    636  SG  CYS A 338      -9.991  14.515  11.993  1.00 18.62           S
ATOM   2670  N   CYS B 338      -8.615  45.421 -38.743  1.00 19.16           N
ATOM   2671  CA  CYS B 338      -8.357  46.490 -37.804  1.00  9.38           C
ATOM   2672  C   CYS B 338      -6.940  46.477 -37.283  1.00 13.38           C
ATOM   2673  O   CYS B 338      -6.183  45.531 -37.509  1.00 12.78           O
ATOM   2674  CB  CYS B 338      -9.340  46.408 -36.640  1.00  6.81           C
ATOM   2675  SG  CYS B 338     -11.059  46.627 -37.137  1.00 18.08           S
HETATM 4105  C10 VSP A 600     -12.073  10.268   5.878  1.00 15.21           C
HETATM 4106  C13 VSP A 600      -9.415  12.009   8.138  1.00 17.34           C
HETATM 4107  C17 VSP A 600      -9.173  15.500   9.538  1.00 18.06           C
HETATM 4108  C21 VSP A 600      -9.724  11.447  10.488  1.00 13.59           C
HETATM 4109  C22 VSP A 600     -10.796  10.611  10.141  1.00 13.09           C
HETATM 4110  O20 VSP A 600      -6.512  15.097   9.546  1.00 36.73           O
HETATM 4111  N01 VSP A 600     -13.584  12.637   7.376  1.00  2.48           N
HETATM 4112  C11 VSP A 600     -10.820  11.061   6.336  1.00 16.83           C
HETATM 4113  C02 VSP A 600     -14.234  11.590   6.685  1.00  0.00           C
HETATM 4114  C12 VSP A 600     -10.480  11.165   7.808  1.00 19.86           C
HETATM 4115  C03 VSP A 600     -13.499  10.594   6.048  1.00  7.99           C
HETATM 4116  C23 VSP A 600     -11.178  10.454   8.807  1.00 13.40           C
HETATM 4117  N24 VSP A 600     -15.578   9.539   5.351  1.00  5.37           N
HETATM 4118  C04 VSP A 600     -14.214   9.541   5.368  1.00  6.59           C
HETATM 4119  C14 VSP A 600      -9.040  12.146   9.490  1.00 20.32           C
HETATM 4120  N05 VSP A 600     -13.277   8.720   4.845  1.00 14.36           N
HETATM 4121  N15 VSP A 600      -7.941  13.011   9.844  1.00 30.28           N
HETATM 4122  C25 VSP A 600     -16.266  10.507   5.972  1.00  8.74           C
HETATM 4123  N26 VSP A 600     -15.595  11.531   6.617  1.00  7.31           N
HETATM 4124  C06 VSP A 600     -13.564   7.524   4.064  1.00 21.61           C
HETATM 4125  S16 VSP A 600      -7.758  14.494   9.137  1.00 29.97           S
HETATM 4126  C07 VSP A 600     -12.688   6.369   4.532  1.00 18.74           C
HETATM 4127  C08 VSP A 600     -13.266   7.825   2.614  1.00 11.95           C
HETATM 4128  C18 VSP A 600      -9.254  15.824  11.014  1.00 27.59           C
HETATM 4129  N09 VSP A 600     -12.038   9.143   5.150  1.00 15.70           N
HETATM 4130  O19 VSP A 600      -7.560  14.438   7.701  1.00 31.30           O
HETATM 4131  C10 VSP B 600      -9.100  42.047 -31.076  1.00 15.42           C
HETATM 4132  C13 VSP B 600     -11.662  43.913 -33.270  1.00 14.18           C
HETATM 4133  C17 VSP B 600     -11.854  47.483 -34.510  1.00 33.18           C
HETATM 4134  C21 VSP B 600     -11.316  43.454 -35.629  1.00 18.55           C
HETATM 4135  C22 VSP B 600     -10.247  42.605 -35.296  1.00 16.17           C
HETATM 4136  O20 VSP B 600     -14.506  47.137 -34.713  1.00 28.59           O
HETATM 4137  N01 VSP B 600      -7.431  44.349 -32.620  1.00  3.96           N
HETATM 4138  C11 VSP B 600     -10.309  42.931 -31.468  1.00 13.44           C
HETATM 4139  C02 VSP B 600      -6.859  43.222 -31.970  1.00  9.31           C
HETATM 4140  C12 VSP B 600     -10.611  43.066 -32.949  1.00 16.06           C
HETATM 4141  C03 VSP B 600      -7.659  42.259 -31.325  1.00 10.60           C
HETATM 4142  C23 VSP B 600      -9.886  42.396 -33.954  1.00 14.26           C
HETATM 4143  N24 VSP B 600      -5.646  41.044 -30.732  1.00 11.16           N
HETATM 4144  C04 VSP B 600      -7.014  41.143 -30.686  1.00 10.43           C
HETATM 4145  C14 VSP B 600     -12.022  44.125 -34.611  1.00 21.96           C
HETATM 4146  N05 VSP B 600      -8.006  40.380 -30.149  1.00 12.14           N
HETATM 4147  N15 VSP B 600     -13.131  45.011 -34.903  1.00 21.91           N
HETATM 4148  C25 VSP B 600      -4.907  41.964 -31.365  1.00 14.73           C
HETATM 4149  N26 VSP B 600      -5.504  43.046 -31.973  1.00  6.61           N
HETATM 4150  C06 VSP B 600      -7.852  39.153 -29.374  1.00 16.29           C
HETATM 4151  S16 VSP B 600     -13.310  46.515 -34.202  1.00 24.65           S
HETATM 4152  C07 VSP B 600      -8.917  38.154 -29.784  1.00 12.46           C
HETATM 4153  C08 VSP B 600      -8.172  39.532 -27.956  1.00 15.34           C
HETATM 4154  C18 VSP B 600     -11.699  47.833 -35.975  1.00 29.79           C
HETATM 4155  N09 VSP B 600      -9.203  40.922 -30.367  1.00 10.22           N
HETATM 4156  O19 VSP B 600     -13.564  46.514 -32.785  1.00 21.41           O
""",
        "vsp.cif" : """data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
VSP        VSP 'N-(3-{[4-amino-1-(propan-2-yl)-1H-pyrazolo[3,4-d]pyrimidin-3-yl]methyl}phenyl)ethanesulfonamide' ligand 48 26 .
#
data_comp_VSP
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
VSP         C10    C   CR5   .         -1.0980    0.5259   -1.4219
VSP         C13    C   CR16  .          1.7614   -0.3172    0.5726
VSP         C17    C   CH2   .          5.6014   -2.4120    3.1666
VSP         C21    C   CR16  .          1.7453    0.7693    2.7632
VSP         C22    C   CR16  .          0.8735    1.7367    2.2769
VSP         O20    O   OS    .          5.2933   -0.3503    1.4516
VSP         N01    N   NH2   .         -0.3458    2.0703   -4.1108
VSP         C11    C   CH2   .          0.4050    0.6157   -1.3625
VSP         C02    C   CR6   .         -1.6892    1.7101   -3.6579
VSP         C12    C   CR6   .          0.8759    0.6638    0.0868
VSP         C03    C   CR56  .         -1.9069    1.0097   -2.4216
VSP         C23    C   CR16  .          0.4365    1.6918    0.9409
VSP         N24    N   N     .         -4.2441    1.0766   -2.8909
VSP         C04    C   CR56  .         -3.2296    0.7025   -2.0603
VSP         C14    C   CR6   .          2.1976   -0.2670    1.9141
VSP         N05    N   NR5   .         -3.1716    0.0434   -0.8597
VSP         N15    N   NC1   .          3.1146   -1.2864    2.4232
VSP         C25    C   CR16  .         -3.9926    1.7210   -4.0202
VSP         N26    N   N     .         -2.7539    2.0345   -4.4038
VSP         C06    C   CH1   .         -4.2964   -0.4622   -0.1100
VSP         S16    S   S     .          4.7064   -0.9300    2.6495
VSP         C07    C   CH3   .         -5.2645    0.6852    0.2164
VSP         C08    C   CH3   .         -5.0317   -1.5364   -0.9389
VSP         C18    C   CH3   .          6.7580   -2.6668    2.2013
VSP         N09    N   N     .         -1.8758   -0.0508   -0.4964
VSP         O19    O   OS    .          4.8764    0.1446    3.5931
VSP         H13    H   HCR6  .          2.0684   -1.0255   -0.0103
VSP         H17    H   HCH2  .          5.9543   -2.2830    4.0830
VSP        H17A    H   HCH2  .          4.9925   -3.1829    3.1559
VSP         H21    H   HCR6  .          2.0424    0.8054    3.6816
VSP         H22    H   HCR6  .          0.5656    2.4466    2.8623
VSP        HN01    H   HNH2  .         -0.0669    2.9442   -4.0707
VSP        HN0A    H   HNH2  .          0.2203    1.4223   -4.4341
VSP         H11    H   HCH2  .          0.7930   -0.1642   -1.8001
VSP        H11A    H   HCH2  .          0.7032    1.4281   -1.8250
VSP         H23    H   HCR6  .         -0.1657    2.3626    0.6121
VSP        HN15    H   HNC1  .          2.8065   -2.1149    2.6093
VSP         H25    H   HCR6  .         -4.7360    1.9744   -4.5899
VSP         H06    H   HCH1  .         -3.9722   -0.8645    0.7323
VSP         H07    H   HCH3  .         -5.8747    0.4033    0.9157
VSP        H07A    H   HCH3  .         -4.7525    1.4713    0.5237
VSP        H07B    H   HCH3  .         -5.7736    0.9202   -0.5815
VSP         H08    H   HCH3  .         -5.5202   -1.1073   -1.6737
VSP        H08A    H   HCH3  .         -4.3705   -2.1772   -1.3107
VSP        H08B    H   HCH3  .         -5.6638   -2.0152   -0.3647
VSP         H18    H   HCH3  .          6.4134   -2.7373    1.2907
VSP        H18A    H   HCH3  .          7.3897   -1.9326    2.2516
VSP        H18B    H   HCH3  .          7.2035   -3.4950    2.4420
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
VSP   C10     C11   single        1.507 0.020
VSP   C10     C03   aromatic      1.374 0.020
VSP   C10     N09   aromatic      1.339 0.020
VSP   C13     C12   aromatic      1.408 0.020
VSP   C13     C14   aromatic      1.412 0.020
VSP   C17     S16   single        1.807 0.020
VSP   C17     C18   single        1.528 0.020
VSP   C21     C22   aromatic      1.390 0.020
VSP   C21     C14   aromatic      1.414 0.020
VSP   C22     C23   aromatic      1.406 0.020
VSP   O20     S16   double        1.454 0.020
VSP   N01     C02   single        1.463 0.020
VSP   C11     C12   single        1.525 0.020
VSP   C02     C03   aromatic      1.437 0.020
VSP   C02     N26   aromatic      1.340 0.020
VSP   C12     C23   aromatic      1.407 0.020
VSP   C03     C04   aromatic      1.405 0.020
VSP   N24     C04   aromatic      1.363 0.020
VSP   N24     C25   aromatic      1.324 0.020
VSP   C04     N05   aromatic      1.371 0.020
VSP   C14     N15   single        1.463 0.020
VSP   N05     C06   single        1.443 0.020
VSP   N05     N09   aromatic      1.349 0.020
VSP   N15     S16   single        1.647 0.020
VSP   C25     N26   aromatic      1.334 0.020
VSP   C06     C07   single        1.536 0.020
VSP   C06     C08   single        1.543 0.020
VSP   S16     O19   double        1.440 0.020
VSP   C13     H13   single        0.967 0.020
VSP   C17     H17   single        0.990 0.020
VSP   C17    H17A   single        0.982 0.020
VSP   C21     H21   single        0.966 0.020
VSP   C22     H22   single        0.970 0.020
VSP   N01    HN01   single        0.918 0.020
VSP   N01    HN0A   single        0.919 0.020
VSP   C11     H11   single        0.975 0.020
VSP   C11    H11A   single        0.981 0.020
VSP   C23     H23   single        0.959 0.020
VSP   N15    HN15   single        0.903 0.020
VSP   C25     H25   single        0.970 0.020
VSP   C06     H06   single        0.988 0.020
VSP   C07     H07   single        0.970 0.020
VSP   C07    H07A   single        0.987 0.020
VSP   C07    H07B   single        0.975 0.020
VSP   C08     H08   single        0.981 0.020
VSP   C08    H08A   single        0.993 0.020
VSP   C08    H08B   single        0.979 0.020
VSP   C18     H18   single        0.976 0.020
VSP   C18    H18A   single        0.970 0.020
VSP   C18    H18B   single        0.971 0.020
#
""",
        "linking_test_LYS-ABA-GLU.pdb" : """
ATOM     41  N   LYS S   7     -13.818  -0.993  15.491  1.00 17.94           N
ATOM     42  CA  LYS S   7     -13.936  -0.779  14.054  1.00 17.73           C
ATOM     43  C   LYS S   7     -13.883  -2.110  13.316  1.00 17.51           C
ATOM     44  O   LYS S   7     -14.643  -2.338  12.374  1.00 14.65           O
ATOM     45  CB  LYS S   7     -12.808   0.129  13.558  1.00 19.47           C
ATOM     46  CG  LYS S   7     -12.855   0.415  12.068  1.00 24.22           C
ATOM     47  CD  LYS S   7     -11.746   1.372  11.655  1.00 27.47           C
ATOM     48  CE  LYS S   7     -11.804   1.661  10.162  1.00 31.97           C
ATOM     49  NZ  LYS S   7     -10.718   2.586   9.730  1.00 34.51           N
HETATM   50  N   ABA S   8     -12.984  -2.991  13.745  1.00 15.85           N
HETATM   51  CA  ABA S   8     -12.846  -4.298  13.109  1.00 15.74           C
HETATM   52  C   ABA S   8     -14.142  -5.102  13.246  1.00 14.66           C
HETATM   53  O   ABA S   8     -14.565  -5.778  12.308  1.00 16.07           O
HETATM   54  CB  ABA S   8     -11.676  -5.068  13.736  1.00 17.05           C
HETATM   55  CG  ABA S   8     -11.423  -6.426  13.100  1.00 20.97           C
ATOM     56  N   GLU S   9     -14.778  -5.025  14.413  1.00 14.08           N
ATOM     57  CA  GLU S   9     -16.026  -5.746  14.648  1.00 15.05           C
ATOM     58  C   GLU S   9     -17.147  -5.199  13.765  1.00 15.55           C
ATOM     59  O   GLU S   9     -17.893  -5.956  13.136  1.00 15.55           O
ATOM     60  CB  GLU S   9     -16.441  -5.626  16.114  1.00 16.16           C
ATOM     61  CG  GLU S   9     -15.413  -6.165  17.096  1.00 20.84           C
ATOM     62  CD  GLU S   9     -15.822  -5.942  18.542  1.00 26.26           C
ATOM     63  OE1 GLU S   9     -16.348  -4.851  18.850  1.00 29.12           O
ATOM     64  OE2 GLU S   9     -15.608  -6.850  19.370  1.00 30.05           O
""",
        "linking_test_ASN-NAG-altloc1.pdb" : """
ATOM   1989  N   ASN A 270      20.738  61.827 110.156  1.00 11.35           N
ATOM   1990  CA  ASN A 270      20.866  63.240 110.482  1.00 11.45           C
ATOM   1991  C   ASN A 270      20.238  64.075 109.377  1.00 11.25           C
ATOM   1992  O   ASN A 270      20.902  64.453 108.408  1.00 10.95           O
ATOM   1993  CB  ASN A 270      22.330  63.643 110.648  1.00 12.36           C
ATOM   1994  CG  ASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1 ASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2 ASN A 270      23.440  65.083 112.236  1.00 13.89           N
HETATM 3221  C1 ANAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 ANAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 ANAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 ANAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 ANAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 ANAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 ANAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 ANAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 ANAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 ANAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 ANAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 ANAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 ANAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 ANAG A 435      22.047  63.959 114.195  1.00 19.58           O
HETATM 3221  C1 BNAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 BNAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 BNAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 BNAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 BNAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 BNAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 BNAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 BNAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 BNAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 BNAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 BNAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 BNAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 BNAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 BNAG A 435      22.047  63.959 114.195  1.00 19.58           O
""",
        "linking_test_ASN-NAG-altloc2.pdb" : """
ATOM   1989  N  AASN A 270      20.738  61.827 110.156  1.00 11.35           N
ATOM   1990  CA AASN A 270      20.866  63.240 110.482  1.00 11.45           C
ATOM   1991  C  AASN A 270      20.238  64.075 109.377  1.00 11.25           C
ATOM   1992  O  AASN A 270      20.902  64.453 108.408  1.00 10.95           O
ATOM   1993  CB AASN A 270      22.330  63.643 110.648  1.00 12.36           C
ATOM   1994  CG AASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1AASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2AASN A 270      23.440  65.083 112.236  1.00 13.89           N
ATOM   1989  N  BASN A 270      20.738  61.827 110.156  1.00 11.35           N
ATOM   1990  CA BASN A 270      20.866  63.240 110.482  1.00 11.45           C
ATOM   1991  C  BASN A 270      20.238  64.075 109.377  1.00 11.25           C
ATOM   1992  O  BASN A 270      20.902  64.453 108.408  1.00 10.95           O
ATOM   1993  CB BASN A 270      22.330  63.643 110.648  1.00 12.36           C
ATOM   1994  CG BASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1BASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2BASN A 270      23.440  65.083 112.236  1.00 13.89           N
HETATM 3221  C1 ANAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 ANAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 ANAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 ANAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 ANAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 ANAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 ANAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 ANAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 ANAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 ANAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 ANAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 ANAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 ANAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 ANAG A 435      22.047  63.959 114.195  1.00 19.58           O
HETATM 3221  C1 BNAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 BNAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 BNAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 BNAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 BNAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 BNAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 BNAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 BNAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 BNAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 BNAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 BNAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 BNAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 BNAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 BNAG A 435      22.047  63.959 114.195  1.00 19.58           O
""",
        "linking_test_ASN-NAG-altloc3.pdb" : """
ATOM   1989  N   ASN A 270      20.738  61.827 110.156  1.00 11.35           N
ATOM   1990  CA  ASN A 270      20.866  63.240 110.482  1.00 11.45           C
ATOM   1991  C   ASN A 270      20.238  64.075 109.377  1.00 11.25           C
ATOM   1992  O   ASN A 270      20.902  64.453 108.408  1.00 10.95           O
ATOM   1993  CB  ASN A 270      22.330  63.643 110.648  1.00 12.36           C
ATOM   1994  CG AASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1AASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2AASN A 270      23.440  65.083 112.236  1.00 13.89           N
ATOM   1994  CG BASN A 270      22.475  64.985 111.329  1.00 13.25           C
ATOM   1995  OD1BASN A 270      21.731  65.920 111.033  1.00 13.36           O
ATOM   1996  ND2BASN A 270      23.440  65.083 112.236  1.00 13.89           N
HETATM 3221  C1 ANAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 ANAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 ANAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 ANAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 ANAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 ANAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 ANAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 ANAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 ANAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 ANAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 ANAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 ANAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 ANAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 ANAG A 435      22.047  63.959 114.195  1.00 19.58           O
HETATM 3221  C1 BNAG A 435      23.764  66.364 112.833  1.00 15.11           C
HETATM 3222  C2 BNAG A 435      23.381  66.347 114.325  1.00 16.09           C
HETATM 3223  C3 BNAG A 435      23.832  67.642 115.011  1.00 16.01           C
HETATM 3224  C4 BNAG A 435      25.317  67.871 114.758  1.00 15.48           C
HETATM 3225  C5 BNAG A 435      25.583  67.877 113.254  1.00 14.91           C
HETATM 3226  C6 BNAG A 435      27.059  68.022 112.960  1.00 14.68           C
HETATM 3227  C7 BNAG A 435      21.392  64.986 114.378  1.00 18.41           C
HETATM 3228  C8 BNAG A 435      19.882  64.903 114.521  1.00 19.12           C
HETATM 3229  N2 BNAG A 435      21.944  66.193 114.461  1.00 17.50           N
HETATM 3230  O3 BNAG A 435      23.601  67.552 116.408  1.00 16.86           O
HETATM 3231  O4 BNAG A 435      25.725  69.107 115.324  1.00 15.49           O
HETATM 3232  O5 BNAG A 435      25.166  66.624 112.676  1.00 14.58           O
HETATM 3233  O6 BNAG A 435      27.784  66.911 113.470  1.00 14.35           O
HETATM 3234  O7 BNAG A 435      22.047  63.959 114.195  1.00 19.58           O
""",
        "linking_test_DT-4MF-DA.pdb" : """
ATOM   2093  P    DT C  26       0.288   1.691 -25.901  1.00 51.23           P
ATOM   2094  OP1  DT C  26       1.525   1.058 -26.431  1.00 51.10           O
ATOM   2095  OP2  DT C  26      -0.538   2.550 -26.800  1.00 50.36           O
ATOM   2096  O5'  DT C  26      -0.678   0.604 -25.234  1.00 47.26           O
ATOM   2097  C5'  DT C  26      -0.147  -0.278 -24.258  1.00 44.30           C
ATOM   2098  C4'  DT C  26      -1.108  -0.531 -23.117  1.00 40.69           C
ATOM   2099  O4'  DT C  26      -1.860   0.647 -22.788  1.00 39.38           O
ATOM   2100  C3'  DT C  26      -2.247  -1.484 -23.405  1.00 40.36           C
ATOM   2101  O3'  DT C  26      -1.754  -2.822 -23.434  1.00 39.62           O
ATOM   2102  C2'  DT C  26      -3.214  -1.187 -22.243  1.00 39.03           C
ATOM   2103  C1'  DT C  26      -2.726   0.186 -21.778  1.00 36.96           C
ATOM   2104  N1   DT C  26      -3.767   1.212 -21.523  1.00 36.37           N
ATOM   2105  C2   DT C  26      -4.324   1.300 -20.260  1.00 35.70           C
ATOM   2106  O2   DT C  26      -4.018   0.575 -19.325  1.00 35.41           O
ATOM   2107  N3   DT C  26      -5.271   2.290 -20.125  1.00 34.72           N
ATOM   2108  C4   DT C  26      -5.714   3.176 -21.081  1.00 34.88           C
ATOM   2109  O4   DT C  26      -6.571   4.022 -20.842  1.00 35.55           O
ATOM   2110  C5   DT C  26      -5.100   3.038 -22.379  1.00 35.49           C
ATOM   2111  C7   DT C  26      -5.637   3.830 -23.539  1.00 34.73           C
ATOM   2112  C6   DT C  26      -4.166   2.075 -22.536  1.00 36.24           C
HETATM 2113  O2P 4MF C  27      -3.441  -4.407 -23.957  1.00 41.52           O
HETATM 2114  P   4MF C  27      -2.133  -4.347 -23.209  1.00 40.92           P
HETATM 2115  O1P 4MF C  27      -0.830  -4.701 -23.874  1.00 40.35           O
HETATM 2116  O5' 4MF C  27      -2.273  -5.251 -21.880  1.00 43.61           O
HETATM 2117  C5' 4MF C  27      -2.804  -4.709 -20.672  1.00 45.90           C
HETATM 2118  C4' 4MF C  27      -4.324  -4.877 -20.528  1.00 48.09           C
HETATM 2119  C3' 4MF C  27      -5.014  -5.783 -21.538  1.00 48.57           C
HETATM 2120  O3' 4MF C  27      -5.183  -7.086 -20.990  1.00 49.62           O
HETATM 2121  C2' 4MF C  27      -6.376  -5.166 -21.819  1.00 49.28           C
HETATM 2122  C1' 4MF C  27      -6.327  -3.804 -21.133  1.00 49.96           C
HETATM 2123  O4' 4MF C  27      -4.988  -3.610 -20.660  1.00 50.00           O
HETATM 2124  N1  4MF C  27      -6.624  -2.783 -22.129  1.00 50.63           N
HETATM 2125  C7A 4MF C  27      -7.134  -1.592 -21.893  1.00 50.94           C
HETATM 2126  C7  4MF C  27      -7.511  -0.950 -20.876  1.00 50.46           C
HETATM 2127  C6  4MF C  27      -7.823   0.199 -21.023  1.00 50.24           C
HETATM 2128  C5  4MF C  27      -8.031   0.699 -22.101  1.00 50.59           C
HETATM 2129  C4  4MF C  27      -7.681   0.140 -23.189  1.00 50.69           C
HETATM 2130  C4M 4MF C  27      -7.799   0.755 -24.356  1.00 50.25           C
HETATM 2131  C3A 4MF C  27      -7.199  -1.092 -23.071  1.00 50.73           C
HETATM 2132  C3  4MF C  27      -6.765  -1.912 -23.925  1.00 51.26           C
HETATM 2133  C2  4MF C  27      -6.425  -2.920 -23.332  1.00 50.64           C
ATOM   2134  P    DA C  28      -5.472  -8.466 -21.727  1.00 49.25           P
ATOM   2135  OP1  DA C  28      -4.838  -9.205 -20.608  1.00 49.49           O
ATOM   2136  OP2  DA C  28      -4.975  -8.640 -23.106  1.00 48.81           O
ATOM   2137  O5'  DA C  28      -7.046  -8.727 -21.708  1.00 49.80           O
ATOM   2138  C5'  DA C  28      -7.693  -9.165 -20.528  1.00 51.19           C
ATOM   2139  C4'  DA C  28      -9.182  -8.886 -20.616  1.00 52.05           C
ATOM   2140  O4'  DA C  28      -9.387  -7.491 -20.963  1.00 51.89           O
ATOM   2141  C3'  DA C  28      -9.926  -9.711 -21.664  1.00 52.66           C
ATOM   2142  O3'  DA C  28     -11.105 -10.227 -21.073  1.00 54.53           O
ATOM   2143  C2'  DA C  28     -10.233  -8.695 -22.768  1.00 51.97           C
ATOM   2144  C1'  DA C  28     -10.388  -7.412 -21.957  1.00 50.53           C
ATOM   2145  N9   DA C  28     -10.194  -6.169 -22.690  1.00 49.74           N
ATOM   2146  C8   DA C  28      -9.269  -5.919 -23.662  1.00 49.76           C
ATOM   2147  N7   DA C  28      -9.328  -4.705 -24.145  1.00 49.73           N
ATOM   2148  C5   DA C  28     -10.361  -4.112 -23.430  1.00 49.84           C
ATOM   2149  C6   DA C  28     -10.936  -2.820 -23.456  1.00 49.81           C
ATOM   2150  N6   DA C  28     -10.525  -1.842 -24.264  1.00 49.75           N
ATOM   2151  N1   DA C  28     -11.956  -2.559 -22.614  1.00 49.98           N
ATOM   2152  C2   DA C  28     -12.373  -3.529 -21.793  1.00 49.56           C
ATOM   2153  N3   DA C  28     -11.917  -4.774 -21.688  1.00 49.25           N
ATOM   2154  C4   DA C  28     -10.901  -5.004 -22.529  1.00 48.98           C
""",
    "4mf.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
4MF        4MF '1-(2-deoxy-5-O-phosphono-beta-D-erythro-pentofuranosyl)-4-methyl-1H-indole' ligand 38 22 .
#
data_comp_4MF
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
4MF         O3P    O   OP    .          0.3880   -6.5014   -0.0977
4MF         O2P    O   OP    .          2.8128   -6.3045    0.2744
4MF         P      P   P     .          1.4780   -5.6201    0.4475
4MF         O1P    O   O     .          1.2312   -5.3632    1.9128
4MF        O5'     O   O2    .          1.4842   -4.2849   -0.3118
4MF        C5'     C   CH2   .          0.4523   -3.3972   -0.1241
4MF        C4'     C   CH1   .          0.4545   -2.3849   -1.2304
4MF        C3'     C   CH1   .          1.6882   -1.6217   -1.1949
4MF        O3'     O   OH1   .          2.4943   -1.9836   -2.2794
4MF        C2'     C   CH2   .          1.3646   -0.1729   -1.2778
4MF        C1'     C   CH1   .         -0.1502   -0.1860   -1.2349
4MF        O4'     O   O2    .         -0.6137   -1.5216   -1.0689
4MF         N1     N   NR5   .         -0.5939    0.6028   -0.1623
4MF         C7A    C   CR56  .         -0.7732    1.9661   -0.2075
4MF         C7     C   CR16  .         -0.5799    2.8750   -1.2494
4MF         C6     C   CR16  .         -0.8345    4.2251   -1.0458
4MF         C5     C   CR16  .         -1.2864    4.6714    0.2113
4MF         C4     C   CR6   .         -1.4808    3.7274    1.2736
4MF         C4M    C   CH3   .         -1.9588    4.1914    2.6065
4MF         C3A    C   CR56  .         -1.2204    2.3671    1.0551
4MF         C3     C   CR15  .         -1.3024    1.2223    1.8455
4MF         C2     C   CR15  .         -0.9093    0.1452    1.0684
4MF        H5'1    H   HCH2  .         -0.3941   -3.8790   -0.1265
4MF        H5'2    H   HCH2  .          0.5662   -2.9370    0.7439
4MF        H4'     H   HCH1  .          0.3777   -2.8403   -2.0880
4MF        H3'     H   HCH1  .          2.1650   -1.8026   -0.3580
4MF         H3T    H   HOH1  .          3.2930   -2.2169   -1.9870
4MF        H2'1    H   HCH2  .          1.7298    0.3100   -0.5098
4MF        H2'2    H   HCH2  .          1.6862    0.2071   -2.1162
4MF        H1'     H   HCH1  .         -0.5001    0.1731   -2.0687
4MF         H7     H   HCR6  .         -0.2675    2.5671   -2.1190
4MF         H6     H   HCR6  .         -0.6998    4.8655   -1.7753
4MF         H5     H   HCR6  .         -1.4667    5.6202    0.3592
4MF        H4M1    H   HCH3  .         -1.2789    3.9974    3.2823
4MF        H4M2    H   HCH3  .         -2.1218    5.1482    2.5767
4MF        H4M3    H   HCH3  .         -2.7822    3.7299    2.8339
4MF         H3     H   HCR5  .         -1.5857    1.1856    2.7788
4MF         H2     H   HCR5  .         -0.8680   -0.7821    1.3604
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
4MF   O3P     P     deloc         1.504 0.020
4MF   O2P     P     deloc         1.510 0.020
4MF   P       O1P   deloc         1.508 0.020
4MF   P      O5'    single        1.536 0.020
4MF  O5'     C5'    single        1.374 0.020
4MF  C5'     C4'    single        1.500 0.020
4MF  C5'     H5'1   single        0.974 0.020
4MF  C5'     H5'2   single        0.989 0.020
4MF  C4'     C3'    single        1.451 0.020
4MF  C4'     O4'    single        1.383 0.020
4MF  C4'     H4'    single        0.974 0.020
4MF  C3'     C2'    single        1.487 0.020
4MF  C3'     O3'    single        1.399 0.020
4MF  C3'     H3'    single        0.980 0.020
4MF  O3'      H3T   single        0.882 0.020
4MF  C2'     C1'    single        1.516 0.020
4MF  C2'     H2'1   single        0.978 0.020
4MF  C2'     H2'2   single        0.975 0.020
4MF  C1'      N1    single        1.403 0.020
4MF  C1'     O4'    single        1.423 0.020
4MF  C1'     H1'    single        0.973 0.020
4MF   N1      C2    aromatic      1.350 0.020
4MF   N1      C7A   aromatic      1.376 0.020
4MF   C7A     C3A   aromatic      1.398 0.020
4MF   C7A     C7    aromatic      1.396 0.020
4MF   C7      C6    aromatic      1.389 0.020
4MF   C7      H7    single        0.974 0.020
4MF   C6      C5    aromatic      1.408 0.020
4MF   C6      H6    single        0.980 0.020
4MF   C5      C4    aromatic      1.434 0.020
4MF   C5      H5    single        0.977 0.020
4MF   C4      C4M   single        1.490 0.020
4MF   C4      C3A   aromatic      1.402 0.020
4MF   C4M    H4M1   single        0.978 0.020
4MF   C4M    H4M2   single        0.971 0.020
4MF   C4M    H4M3   single        0.971 0.020
4MF   C3A     C3    aromatic      1.394 0.020
4MF   C3      C2    aromatic      1.385 0.020
4MF   C3      H3    single        0.976 0.020
4MF   C2      H2    single        0.973 0.020
""",
    "linking_test_NAG-FU4.pdb" : """
HETATM 9423  C1  NAG D1584       6.798 -15.752  22.420  1.00 22.72           C
HETATM 9424  C2  NAG D1584       5.572 -16.434  23.033  1.00 24.18           C
HETATM 9425  C3  NAG D1584       4.596 -15.354  23.498  1.00 23.13           C
HETATM 9426  C4  NAG D1584       5.251 -14.345  24.447  1.00 23.39           C
HETATM 9427  C5  NAG D1584       6.580 -13.815  23.895  1.00 24.93           C
HETATM 9428  C6  NAG D1584       7.379 -13.088  24.991  1.00 27.57           C
HETATM 9429  C7  NAG D1584       5.155 -18.617  22.021  1.00 29.09           C
HETATM 9430  C8  NAG D1584       4.408 -19.417  20.984  1.00 29.54           C
HETATM 9431  N2  NAG D1584       4.916 -17.305  22.067  1.00 25.34           N
HETATM 9432  O3  NAG D1584       3.454 -15.926  24.089  1.00 23.03           O
HETATM 9433  O4  NAG D1584       4.387 -13.247  24.573  1.00 22.66           O
HETATM 9434  O5  NAG D1584       7.367 -14.889  23.394  1.00 23.22           O
HETATM 9435  O7  NAG D1584       5.951 -19.182  22.779  1.00 31.04           O
HETATM 9436  O6  NAG D1584       7.429 -13.988  26.081  1.00 31.80           O
HETATM 9462  C1  FU4 D1588       8.042 -13.514  27.304  1.00 37.58           C
HETATM 9463  C2  FU4 D1588       7.470 -14.288  28.512  1.00 38.91           C
HETATM 9464  C3  FU4 D1588       7.920 -15.750  28.448  1.00 41.29           C
HETATM 9465  C4  FU4 D1588       9.448 -15.793  28.444  1.00 41.76           C
HETATM 9466  C5  FU4 D1588       9.980 -14.949  27.275  1.00 40.83           C
HETATM 9467  C6  FU4 D1588      11.507 -14.975  27.258  1.00 40.52           C
HETATM 9468  O2  FU4 D1588       6.059 -14.235  28.600  1.00 39.96           O
HETATM 9469  O3  FU4 D1588       7.398 -16.498  29.524  1.00 40.18           O
HETATM 9470  O4  FU4 D1588       9.951 -15.298  29.675  1.00 43.17           O
HETATM 9471  O5  FU4 D1588       9.470 -13.611  27.319  1.00 39.54           O
""",
    "fu4.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
FU4        FU4 '2,6-anhydro-1-deoxy-D-galactitol' ligand 22 10 .
#
data_comp_FU4
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
FU4         C1     C   CH2   .         -0.9215    0.7090   -1.3477
FU4         C2     C   CH1   .          0.5905    0.7090   -1.3477
FU4         C3     C   CH1   .          1.0944    0.7090    0.0773
FU4         C4     C   CH1   .          0.4272   -0.4074    0.8476
FU4         C5     C   CH1   .         -1.0718   -0.3124    0.6757
FU4         C6     C   CH3   .         -1.7403   -1.4199    1.4597
FU4         O2     O   OH1   .          1.0598    1.8584   -2.0114
FU4         O3     O   OH1   .          2.4889    0.5143    0.0821
FU4         O4     O   OH1   .          0.8778   -1.6488    0.3597
FU4         O5     O   O2    .         -1.3966   -0.4405   -0.6867
FU4        H1C1    H   HCH2  .         -1.2461    0.7114   -2.2683
FU4        H1C2    H   HCH2  .         -1.2462    1.5059   -0.8870
FU4         H2     H   HCH1  .          0.9159   -0.0880   -1.8079
FU4         H3     H   HCH1  .          0.8839    1.5648    0.4968
FU4         H4     H   HCH1  .          0.6518   -0.3289    1.7943
FU4         H5     H   HCH1  .         -1.3829    0.5524    1.0045
FU4        H6C1    H   HCH3  .         -2.7057   -1.2758    1.4678
FU4        H6C2    H   HCH3  .         -1.5431   -2.2795    1.0413
FU4        H6C3    H   HCH3  .         -1.4036   -1.4191    2.3759
FU4         HA     H   HOH1  .          1.3561    1.6345   -2.8176
FU4         HB     H   HOH1  .          2.8770    1.1651    0.5444
FU4         HC     H   HOH1  .          1.4344   -2.0133    0.9472
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
FU4   C1      C2    single        1.512 0.020
FU4   C1      O5    single        1.408 0.020
FU4   C1     H1C1   single        0.976 0.020
FU4   C1     H1C2   single        0.976 0.020
FU4   C2      C3    single        1.511 0.020
FU4   C2      O2    single        1.408 0.020
FU4   C2      H2    single        0.976 0.020
FU4   C3      C4    single        1.511 0.020
FU4   C3      O3    single        1.408 0.020
FU4   C3      H3    single        0.976 0.020
FU4   C4      C5    single        1.512 0.020
FU4   C4      O4    single        1.408 0.020
FU4   C4      H4    single        0.976 0.020
FU4   C5      C6    single        1.513 0.020
FU4   C5      O5    single        1.406 0.020
FU4   C5      H5    single        0.976 0.020
FU4   C6     H6C1   single        0.976 0.020
FU4   C6     H6C2   single        0.976 0.020
FU4   C6     H6C3   single        0.976 0.020
FU4   O2      HA    single        0.888 0.020
FU4   O3      HB    single        0.888 0.020
FU4   O4      HC    single        0.888 0.020
""",
    "linking_test_BGC-BGC.pdb" : """
HETATM 2374  C2  BGC B 401       5.858  17.444  53.947  1.00 10.73           C
HETATM 2375  C3  BGC B 401       6.199  18.564  52.974  1.00 10.02           C
HETATM 2376  C4  BGC B 401       4.883  19.259  52.650  1.00  9.66           C
HETATM 2377  C5  BGC B 401       4.352  19.862  53.947  1.00  9.19           C
HETATM 2378  C6  BGC B 401       3.073  20.685  53.802  1.00  9.97           C
HETATM 2379  C1  BGC B 401       5.210  18.029  55.193  1.00 12.05           C
HETATM 2380  O1  BGC B 401       4.688  16.987  55.955  1.00 13.71           O
HETATM 2381  O2  BGC B 401       7.016  16.779  54.404  1.00 12.15           O
HETATM 2382  O3  BGC B 401       6.709  18.077  51.750  1.00 10.15           O
HETATM 2383  O4  BGC B 401       5.090  20.298  51.721  1.00  9.48           O
HETATM 2384  O5  BGC B 401       4.087  18.799  54.833  1.00  8.55           O
HETATM 2385  O6  BGC B 401       2.004  19.920  53.229  1.00  8.55           O
HETATM 2386  C2  BGC B 402       8.524  17.329  50.431  1.00 12.41           C
HETATM 2387  C3  BGC B 402      10.050  17.335  50.257  1.00 12.69           C
HETATM 2388  C4  BGC B 402      10.659  18.712  50.487  1.00 12.30           C
HETATM 2389  C5  BGC B 402      10.171  19.272  51.825  1.00 10.88           C
HETATM 2390  C6  BGC B 402      10.610  20.705  52.119  1.00 10.57           C
HETATM 2391  C1  BGC B 402       8.143  18.001  51.734  1.00 11.41           C
HETATM 2392  O2  BGC B 402       7.987  16.015  50.403  1.00 13.00           O
HETATM 2393  O3  BGC B 402      10.347  17.051  48.908  1.00 15.37           O
HETATM 2394  O4  BGC B 402      12.073  18.535  50.546  1.00 13.31           O
HETATM 2395  O5  BGC B 402       8.743  19.283  51.808  1.00  9.88           O
HETATM 2396  O6  BGC B 402      10.164  21.605  51.121  1.00 11.09           O
""",
    "bgc.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
BGC        BGC 'beta-D-glucopyranose     ' ligand 24 12 .
#
data_comp_BGC
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
BGC         C2     C   CH1   .          1.0181    0.7466   -1.2014
BGC         C3     C   CH1   .         -0.4936    0.7466   -1.2014
BGC         C4     C   CH1   .         -0.9974   -0.4285   -0.3952
BGC         C5     C   CH1   .         -0.4042   -0.3777    0.9944
BGC         C6     C   CH2   .         -0.8493   -1.5929    1.7774
BGC         C1     C   CH1   .          1.5218    0.7466    0.2233
BGC         O1     O   OH1   .          2.9275    0.6794    0.2239
BGC         O2     O   OH1   .          1.4875    1.8963   -1.8652
BGC         O3     O   OH1   .         -0.9629    0.6435   -2.5248
BGC         O4     O   OH1   .         -2.4015   -0.3707   -0.3086
BGC         O5     O   O2    .          0.9993   -0.3659    0.9049
BGC         O6     O   OH1   .         -0.2097   -1.6024    3.0344
BGC         H2     H   HCH1  .          1.3435   -0.0504   -1.6615
BGC         H3     H   HCH1  .         -0.8190    1.5769   -0.8045
BGC         H4     H   HCH1  .         -0.7324   -1.2603   -0.8319
BGC         H5     H   HCH1  .         -0.7081    0.4312    1.4486
BGC        H6C1    H   HCH2  .         -1.8166   -1.5616    1.9048
BGC        H6C2    H   HCH2  .         -0.6121   -2.4019    1.2854
BGC         H1     H   HCH1  .          1.2351    1.5659    0.6699
BGC         HA     H   HOH1  .          3.2455    1.1381    0.9141
BGC         HB     H   HOH1  .          2.0042    1.6551   -2.5454
BGC         HC     H   HOH1  .         -1.4335    1.3679   -2.7292
BGC         HD     H   HOH1  .         -2.7570   -0.9641   -0.8647
BGC         H6     H   HOH1  .         -0.5848   -2.2177    3.5527
#
""",
  "linking_test_LEU-CSY-VAL.pdb" : """
ATOM    470  N   LEU A  64      11.530  10.554   5.290  1.00 15.95           N
ATOM    471  CA  LEU A  64      12.272   9.726   6.192  1.00 21.72           C
ATOM    472  CB  LEU A  64      13.571  10.414   6.526  1.00 18.65           C
ATOM    473  CG  LEU A  64      14.569  10.566   5.400  1.00 16.06           C
ATOM    474  CD1 LEU A  64      15.684  11.496   5.782  1.00 15.79           C
ATOM    475  CD2 LEU A  64      15.142   9.234   5.000  1.00 19.26           C
ATOM    476  C   LEU A  64      11.524   9.463   7.499  1.00 47.31           C
ATOM    477  O   LEU A  64      11.704   8.447   8.057  1.00 19.74           O
ATOM    478  O   CSY A  66       7.079   6.443   8.014  1.00 20.23           O
ATOM    479  C   CSY A  66       7.307   7.455   8.592  1.00 42.84           C
ATOM    480  CA3 CSY A  66       7.949   8.535   7.772  1.00 22.48           C
ATOM    481  N3  CSY A  66       7.755   9.898   8.205  1.00 20.36           N
ATOM    482  C2  CSY A  66       6.699  10.604   7.953  1.00 18.89           C
ATOM    483  O2  CSY A  66       5.661  10.200   7.389  1.00 18.87           O
ATOM    484  CA2 CSY A  66       6.889  11.895   8.529  1.00 17.00           C
ATOM    485  N2  CSY A  66       8.098  11.934   9.059  1.00 17.35           N
ATOM    486  C1  CSY A  66       8.613  10.698   8.850  1.00 23.65           C
ATOM    487  CA1 CSY A  66       9.984  10.257   9.208  1.00 15.87           C
ATOM    488  N   CSY A  66      10.796  10.290   7.982  1.00 16.92           N
ATOM    489  CB1 CSY A  66      10.657  11.187  10.190  1.00 22.76           C
ATOM    490  OG1 CSY A  66      10.041  11.081  11.425  1.00 26.93           O
ATOM    491  CB2 CSY A  66       5.901  12.951   8.489  1.00 21.22           C
ATOM    492  CG2 CSY A  66       6.045  14.317   8.943  1.00 17.77           C
ATOM    493  CD1 CSY A  66       7.186  14.821   9.489  1.00 16.75           C
ATOM    494  CE1 CSY A  66       7.230  16.153   9.841  1.00 17.97           C
ATOM    495  CZ  CSY A  66       6.167  16.999   9.597  1.00 17.84           C
ATOM    496  OH  CSY A  66       6.180  18.295   9.906  1.00 18.10           O
ATOM    497  CE2 CSY A  66       5.053  16.481   9.007  1.00 17.93           C
ATOM    498  CD2 CSY A  66       5.000  15.167   8.689  1.00 16.52           C
ATOM    499  N   VAL A  68       7.337   7.591   9.827  1.00 22.68           N
ATOM    500  CA  VAL A  68       6.874   6.579  10.743  1.00 17.94           C
ATOM    501  CB  VAL A  68       8.031   5.907  11.484  1.00 18.83           C
ATOM    502  CG1 VAL A  68       8.808   5.000  10.546  1.00 18.49           C
ATOM    503  CG2 VAL A  68       8.940   6.946  12.060  1.00 21.81           C
ATOM    504  C   VAL A  68       5.837   7.169  11.692  1.00 18.71           C
ATOM    505  O   VAL A  68       6.046   7.339  12.840  1.00 18.66           O
""",
  "csy.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
CSY        CSY 'Unknown                  ' ligand 36 21 .
#
data_comp_CSY
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
CSY         O      O   O     .          6.5225    8.8391   10.1578
CSY         C      C   C1    .          6.8040    8.0495    9.2786
CSY         CA3    C   CH2   .          7.7300    8.4648    8.1571
CSY         N3     N   NR5   .          7.7692    9.9112    8.0435
CSY         C2     C   CR5   .          6.8095   10.6833    7.5020
CSY         O2     O   OH1   .          5.5652   10.2401    6.9578
CSY         CA2    C   CR5   .          7.2462   12.0320    7.6413
CSY         N2     N   N     .          8.4469   11.9889    8.2606
CSY         C1     C   CR5   .          8.7507   10.7094    8.4972
CSY         CA1    C   CH1   .         10.0456   10.2262    9.1064
CSY         N      N   NH2   .         11.0914   10.3565    8.1586
CSY         CB1    C   C1    .         10.3683   11.0559   10.3320
CSY         OG1    O   O     .          9.6909   10.9426   11.3400
CSY         CB2    C   CH2   .          6.5545   13.2812    7.0704
CSY         CG2    C   CR6   .          6.3534   14.3333    8.2072
CSY         CD1    C   CR16  .          7.3740   15.2590    8.4985
CSY         CE1    C   CR16  .          7.1903   16.2102    9.5124
CSY         CZ     C   CR6   .          5.9969   16.2401   10.2329
CSY         OH     O   OH1   .          5.8390   17.1584   11.3068
CSY         CE2    C   CR16  .          4.9783   15.3171    9.9430
CSY         CD2    C   CR16  .          5.1631   14.3623    8.9250
CSY        HC1     H   H     .          6.4050    7.1514    9.2787
CSY        HA31    H   HCH2  .          7.4089    8.0824    7.3139
CSY        HA32    H   HCH2  .          8.6295    8.1304    8.3411
CSY        H21     H   HOH1  .          5.6829    9.9885    6.1182
CSY        HA11    H   HCH1  .          9.9540    9.2881    9.3635
CSY        HN2     H   HNH2  .         11.6123    9.6000    8.1723
CSY        HN3     H   HNH2  .         11.6050   11.0843    8.3689
CSY        HB11    H   H     .         11.0974   11.7100   10.3027
CSY        HB21    H   HCH2  .          5.6862   13.0307    6.7011
CSY        HB22    H   HCH2  .          7.1096   13.6673    6.3644
CSY        HD11    H   HCR6  .          8.1999   15.2679    7.9588
CSY        HE11    H   HCR6  .          7.9314   16.8265    9.7542
CSY        HH1     H   HOH1  .          5.2232   17.7727   11.0875
CSY        HE21    H   HCR6  .          4.1503   15.3164   10.4733
CSY        HD21    H   HCR6  .          4.4341   13.7482    8.6900
#
""",
  "linking_test_ALA-ALA-ALA.pdb" : """
HETATM    1  N   ALA A   1      -0.424   1.960   3.877  1.00 20.00      A    N+1
HETATM    2  H   ALA A   1       0.452   1.694   3.861  1.00 20.00      A    H
HETATM    3  H2  ALA A   1      -0.472   2.843   4.121  1.00 20.00      A    H
HETATM    4  H3  ALA A   1      -0.888   1.448   4.484  1.00 20.00      A    H
HETATM    5  CA  ALA A   1      -0.994   1.794   2.582  1.00 20.00      A    C
HETATM    6  HA  ALA A   1      -0.584   2.442   1.977  1.00 20.00      A    H
HETATM    7  CB  ALA A   1      -2.476   2.066   2.644  1.00 20.00      A    C
HETATM    8  HB1 ALA A   1      -2.627   2.987   2.945  1.00 20.00      A    H
HETATM    9  HB2 ALA A   1      -2.867   1.945   1.755  1.00 20.00      A    H
HETATM   10  HB3 ALA A   1      -2.896   1.444   3.272  1.00 20.00      A    H
HETATM   11  C   ALA A   1      -0.695   0.397   2.070  1.00 20.00      A    C
HETATM   12  O   ALA A   1      -1.218  -0.566   2.593  1.00 20.00      A    O
HETATM   13  N   ALA A   2       0.165   0.219   0.945  1.00 20.00      A    N
HETATM   14  H   ALA A   2       0.578   0.941   0.568  1.00 20.00      A    H
HETATM   15  CA  ALA A   2       0.460  -1.094   0.477  1.00 20.00      A    C
HETATM   16  HA  ALA A   2      -0.095  -1.682   1.027  1.00 20.00      A    H
HETATM   17  CB  ALA A   2       1.882  -1.540   0.755  1.00 20.00      A    C
HETATM   18  HB1 ALA A   2       2.509  -0.921   0.323  1.00 20.00      A    H
HETATM   19  HB2 ALA A   2       2.039  -1.544   1.721  1.00 20.00      A    H
HETATM   20  HB3 ALA A   2       2.016  -2.441   0.399  1.00 20.00      A    H
HETATM   21  C   ALA A   2      -0.001  -1.357  -0.949  1.00 20.00      A    C
HETATM   22  O   ALA A   2      -1.001  -2.020  -1.152  1.00 20.00      A    O
HETATM   23  N   ALA A   3       0.750  -0.851  -2.042  1.00 20.00      A    N
HETATM   24  H   ALA A   3       1.494  -0.350  -1.890  1.00 20.00      A    H
HETATM   25  CA  ALA A   3       0.361  -1.157  -3.372  1.00 20.00      A    C
HETATM   26  HA  ALA A   3      -0.113  -2.013  -3.369  1.00 20.00      A    H
HETATM   27  CB  ALA A   3       1.587  -1.293  -4.240  1.00 20.00      A    C
HETATM   28  HB1 ALA A   3       2.124  -0.477  -4.175  1.00 20.00      A    H
HETATM   29  HB2 ALA A   3       1.315  -1.431  -5.171  1.00 20.00      A    H
HETATM   30  HB3 ALA A   3       2.118  -2.058  -3.937  1.00 20.00      A    H
HETATM   31  C   ALA A   3      -0.576  -0.091  -3.884  1.00 20.00      A    C
HETATM   32  O   ALA A   3      -1.633  -0.418  -4.488  1.00 20.00      A    O
HETATM   33  OXT ALA A   3      -0.293   1.123  -3.725  1.00 20.00      A    O-1
""",
  "linking_test_3g2j-LLP.pdb" : """
ATOM   5294  N   MET A 679      23.170  19.173  20.583  1.00 27.07           N
ATOM   5295  CA  MET A 679      22.542  17.860  20.608  1.00 27.59           C
ATOM   5296  C   MET A 679      21.147  17.931  21.225  1.00 27.78           C
ATOM   5297  O   MET A 679      20.213  17.315  20.716  1.00 28.06           O
ATOM   5298  CB  MET A 679      23.413  16.867  21.380  1.00 27.20           C
ATOM   5299  CG  MET A 679      24.806  16.705  20.814  1.00 27.16           C
ATOM   5300  SD  MET A 679      25.889  15.923  22.006  1.00 27.84           S
ATOM   5301  CE  MET A 679      27.462  16.037  21.161  1.00 27.18           C
HETATM 5302  N1  LLP A 680      24.297  25.452  30.133  1.00 26.73           N
HETATM 5303  C2  LLP A 680      23.189  24.669  30.383  1.00 27.24           C
HETATM 5304  C2' LLP A 680      22.258  25.076  31.486  1.00 27.21           C
HETATM 5305  C3  LLP A 680      22.947  23.518  29.621  1.00 28.75           C
HETATM 5306  O3  LLP A 680      21.948  22.817  29.843  1.00 28.91           O
HETATM 5307  C4  LLP A 680      23.842  23.169  28.598  1.00 28.00           C
HETATM 5308  C4' LLP A 680      23.628  21.946  27.737  1.00 28.32           C
HETATM 5309  C5  LLP A 680      24.962  23.975  28.379  1.00 27.80           C
HETATM 5310  C6  LLP A 680      25.178  25.121  29.134  1.00 27.34           C
HETATM 5311  C5' LLP A 680      25.930  23.666  27.283  1.00 27.48           C
HETATM 5312  OP4 LLP A 680      25.584  24.381  26.078  1.00 28.12           O
HETATM 5313  P   LLP A 680      26.626  24.603  24.884  1.00 30.04           P
HETATM 5314  OP1 LLP A 680      25.728  25.192  23.877  1.00 30.41           O
HETATM 5315  OP2 LLP A 680      27.647  25.508  25.455  1.00 29.33           O
HETATM 5316  OP3 LLP A 680      27.132  23.253  24.540  1.00 28.52           O
HETATM 5317  N   LLP A 680      21.015  18.702  22.304  1.00 28.08           N
HETATM 5318  CA  LLP A 680      19.732  18.900  22.970  1.00 28.30           C
HETATM 5319  CB  LLP A 680      19.877  19.833  24.180  1.00 28.35           C
HETATM 5320  CG  LLP A 680      20.859  19.371  25.247  1.00 28.13           C
HETATM 5321  CD  LLP A 680      21.068  20.451  26.293  1.00 28.09           C
HETATM 5322  CE  LLP A 680      22.289  20.159  27.151  1.00 29.38           C
HETATM 5323  NZ  LLP A 680      22.551  21.255  28.114  1.00 29.67           N
HETATM 5324  C   LLP A 680      18.699  19.460  21.992  1.00 28.76           C
HETATM 5325  O   LLP A 680      17.553  18.996  21.953  1.00 28.82           O
ATOM   5326  N   PHE A 681      19.111  20.458  21.211  1.00 28.66           N
ATOM   5327  CA  PHE A 681      18.259  21.049  20.185  1.00 29.08           C
ATOM   5328  C   PHE A 681      17.927  20.071  19.069  1.00 29.48           C
ATOM   5329  O   PHE A 681      16.805  20.055  18.579  1.00 29.87           O
ATOM   5330  CB  PHE A 681      18.920  22.265  19.557  1.00 28.88           C
ATOM   5331  CG  PHE A 681      18.889  23.490  20.417  1.00 28.62           C
ATOM   5332  CD1 PHE A 681      17.751  24.287  20.474  1.00 28.27           C
ATOM   5333  CD2 PHE A 681      20.013  23.867  21.142  1.00 28.47           C
ATOM   5334  CE1 PHE A 681      17.732  25.427  21.255  1.00 28.69           C
ATOM   5335  CE2 PHE A 681      20.003  25.007  21.930  1.00 28.03           C
ATOM   5336  CZ  PHE A 681      18.871  25.787  21.985  1.00 28.57           C
""",
  "linking_test_ASN-NAG-altloc4.pdb" : """
HETATM  125  C1 ANAG A   1      30.978  40.626 -25.446  0.50 98.96      E    C
HETATM  126  C2 ANAG A   1      30.428  41.897 -26.138  0.50100.81      E    C
HETATM  127  N2 ANAG A   1      29.783  41.751 -27.447  0.50 94.26      E    N
HETATM  128  C7 ANAG A   1      28.807  42.561 -27.924  0.50 89.45      E    C
HETATM  129  O7 ANAG A   1      28.325  42.377 -29.035  0.50 83.16      E    O
HETATM  130  C8 ANAG A   1      28.252  43.720 -27.128  0.50 87.63      E    C
HETATM  131  C3 ANAG A   1      31.618  42.847 -26.221  0.50102.77      E    C
HETATM  132  O3 ANAG A   1      31.420  43.888 -27.159  0.50102.25      E    O
HETATM  133  C4 ANAG A   1      31.840  43.369 -24.802  0.50103.37      E    C
HETATM  134  O4 ANAG A   1      32.969  44.217 -24.782  0.50104.91      E    O
HETATM  135  C5 ANAG A   1      32.040  42.223 -23.785  0.50 94.99      E    C
HETATM  136  C6 ANAG A   1      31.564  42.633 -22.381  0.50 86.60      E    C
HETATM  137  O6 ANAG A   1      32.632  42.591 -21.462  0.50 77.76      E    O
HETATM  138  O5 ANAG A   1      31.458  40.954 -24.130  0.50 99.19      E    O
HETATM  139  C1 BNAG A   1      30.271  40.925 -24.108  0.50 98.96      E    C
HETATM  140  C2 BNAG A   1      31.415  41.878 -23.684  0.50100.81      E    C
HETATM  141  N2 BNAG A   1      32.048  41.664 -22.378  0.50 94.26      E    N
HETATM  142  C7 BNAG A   1      33.331  41.981 -22.079  0.50 89.45      E    C
HETATM  143  O7 BNAG A   1      33.784  41.778 -20.959  0.50 83.16      E    O
HETATM  144  C8 BNAG A   1      34.276  42.586 -23.092  0.50 87.63      E    C
HETATM  145  C3 BNAG A   1      30.816  43.277 -23.775  0.50102.77      E    C
HETATM  146  O3 BNAG A   1      31.569  44.240 -23.061  0.50102.25      E    O
HETATM  147  C4 BNAG A   1      30.718  43.597 -25.266  0.50103.37      E    C
HETATM  148  O4 BNAG A   1      30.115  44.862 -25.440  0.50104.91      E    O
HETATM  149  C5 BNAG A   1      29.907  42.532 -26.038  0.50 94.99      E    C
HETATM  150  C6 BNAG A   1      30.373  42.426 -27.500  0.50 86.60      E    C
HETATM  151  O6 BNAG A   1      29.320  42.743 -28.381  0.50 77.76      E    O
HETATM  152  O5 BNAG A   1      29.866  41.216 -25.458  0.50 99.19      E    O
ATOM    745  N   ASN A 116      27.207  36.475 -25.453  1.00 62.15      A    N
ATOM    746  CA AASN A 116      28.475  37.144 -25.135  0.50 61.77      A    C
ATOM    747  CB AASN A 116      28.182  38.336 -24.223  0.50 69.02      A    C
ATOM    748  CG AASN A 116      29.037  39.555 -24.495  0.50 79.45      A    C
ATOM    749  OD1AASN A 116      28.656  40.644 -24.047  0.50 83.18      A    O
ATOM    750  ND2AASN A 116      30.178  39.419 -25.218  0.50 87.72      A    N
ATOM    751  C   ASN A 116      29.353  36.143 -24.389  1.00 53.28      A    C
ATOM    752  O   ASN A 116      29.778  36.380 -23.266  1.00 48.83      A    O
ATOM    753  CA BASN A 116      28.575  37.144 -25.135  0.50 61.77      A    C
ATOM    754  CB BASN A 116      28.282  38.336 -24.223  0.50 69.02      A    C
ATOM    755  CG BASN A 116      29.487  39.201 -23.920  0.50 79.45      A    C
ATOM    756  OD1BASN A 116      30.611  38.710 -24.084  0.50 83.18      A    O
ATOM    757  ND2BASN A 116      29.299  40.471 -23.478  0.50 87.72      A    N
""",
  "linking_test_ASN_A-NAG_B.pdb" : """
ATOM   2279  N   ASN A 303       7.500  12.479   6.505  1.00 50.55           N
ATOM   2280  CA  ASN A 303       8.366  11.834   7.487  1.00 48.33           C
ATOM   2281  C   ASN A 303       7.525  11.114   8.549  1.00 49.96           C
ATOM   2282  O   ASN A 303       6.418  11.545   8.868  1.00 49.99           O
ATOM   2283  CB  ASN A 303       9.328  12.845   8.125  1.00 49.14           C
ATOM   2284  CG  ASN A 303      10.450  13.261   7.184  1.00 51.22           C
ATOM   2285  OD1 ASN A 303      11.016  12.431   6.474  1.00 50.95           O
ATOM   2286  ND2 ASN A 303      10.787  14.549   7.186  1.00 55.18           N
TER
HETATM11738  C1  NAG B1349      11.771  15.301   6.383  1.00 64.93           C
HETATM11739  C2  NAG B1349      11.260  16.717   6.129  1.00 71.42           C
HETATM11740  C3  NAG B1349      12.340  17.588   5.490  1.00 76.60           C
HETATM11741  C4  NAG B1349      13.598  17.552   6.354  1.00 85.65           C
HETATM11742  C5  NAG B1349      14.040  16.100   6.541  1.00 77.74           C
HETATM11743  C6  NAG B1349      15.243  16.002   7.471  1.00 78.29           C
HETATM11744  C7  NAG B1349       8.880  17.098   5.780  1.00 66.42           C
HETATM11745  C8  NAG B1349       8.026  16.100   6.508  1.00 39.77           C
HETATM11746  N2  NAG B1349      10.056  16.673   5.317  1.00 65.92           N
HETATM11747  O3  NAG B1349      11.875  18.912   5.361  1.00 77.58           O
HETATM11748  O4  NAG B1349      14.644  18.346   5.809  1.00 96.99           O
HETATM11749  O5  NAG B1349      12.999  15.308   7.087  1.00 74.70           O
HETATM11750  O6  NAG B1349      16.436  16.098   6.729  1.00 77.54           O
HETATM11751  O7  NAG B1349       8.483  18.252   5.630  1.00 80.43           O
""",
  "linking_test_nstd_rna_dna_h_bond.pdb" : """
HETATM  100  P   A44 A   5      15.790  -6.178  18.620  1.00 33.68           P
HETATM  101  OP2 A44 A   5      16.475  -4.926  18.273  1.00 37.12           O
HETATM  102  OP1 A44 A   5      16.553  -7.268  19.263  1.00 37.20           O
HETATM  103  O5' A44 A   5      15.148  -6.673  17.251  1.00 35.19           O
HETATM  104  C5' A44 A   5      14.132  -7.683  17.238  1.00 36.66           C
HETATM  105  C4' A44 A   5      13.552  -7.828  15.849  1.00 38.49           C
HETATM  106  O4' A44 A   5      12.637  -6.734  15.559  1.00 39.35           O
HETATM  107  C3' A44 A   5      14.540  -7.792  14.691  1.00 39.06           C
HETATM  108  O3' A44 A   5      15.263  -8.992  14.487  1.00 39.12           O
HETATM  109  C2' A44 A   5      13.644  -7.425  13.514  1.00 40.78           C
HETATM  110  O2' A44 A   5      12.914  -8.512  12.973  1.00 41.20           O
HETATM  111  C1' A44 A   5      12.679  -6.437  14.168  1.00 38.59           C
HETATM  112  N9  A44 A   5      13.167  -5.067  13.981  1.00 37.63           N
HETATM  113  C8  A44 A   5      13.886  -4.289  14.859  1.00 38.01           C
HETATM  114  N7  A44 A   5      14.200  -3.106  14.380  1.00 38.45           N
HETATM  115  C5  A44 A   5      13.643  -3.103  13.107  1.00 36.06           C
HETATM  116  C6  A44 A   5      13.607  -2.138  12.099  1.00 34.61           C
HETATM  117  N6  A44 A   5      14.159  -0.933  12.209  1.00 36.45           N
HETATM  118  N1  A44 A   5      12.967  -2.450  10.958  1.00 33.41           N
HETATM  119  C2  A44 A   5      12.399  -3.647  10.845  1.00 32.30           C
HETATM  120  N3  A44 A   5      12.361  -4.638  11.717  1.00 34.10           N
HETATM  121  C4  A44 A   5      13.006  -4.301  12.846  1.00 36.44           C
HETATM  122  CA' A44 A   5      12.063  -8.031  11.930  1.00 46.36           C
HETATM  123  CD' A44 A   5      11.616 -11.212  12.261  1.00 58.12           C
HETATM  124  OC' A44 A   5      10.850 -10.006  12.555  1.00 53.39           O
HETATM  125  CB' A44 A   5      11.090  -9.064  11.485  1.00 49.86           C
HETATM  478  P   U36 B  20      11.477   3.696   3.004  1.00 42.30           P
HETATM  479  OP1 U36 B  20      10.335   4.004   2.145  1.00 44.51           O
HETATM  480  OP2 U36 B  20      11.842   4.626   4.093  1.00 45.81           O
HETATM  481  O5' U36 B  20      11.173   2.296   3.679  1.00 42.73           O
HETATM  482  C5' U36 B  20      10.757   1.210   2.877  1.00 42.45           C
HETATM  483  C4' U36 B  20      10.319   0.065   3.737  1.00 42.68           C
HETATM  484  O4' U36 B  20      11.457  -0.572   4.347  1.00 42.31           O
HETATM  485  C3' U36 B  20       9.406   0.368   4.909  1.00 41.73           C
HETATM  486  O3' U36 B  20       8.059   0.493   4.530  1.00 44.16           O
HETATM  487  C2' U36 B  20       9.567  -0.885   5.739  1.00 42.76           C
HETATM  488  O2' U36 B  20       8.836  -1.946   5.173  1.00 46.09           O
HETATM  489  C1' U36 B  20      11.046  -1.178   5.554  1.00 41.52           C
HETATM  490  N1  U36 B  20      11.789  -0.588   6.668  1.00 40.08           N
HETATM  491  C2  U36 B  20      11.863  -1.336   7.813  1.00 41.61           C
HETATM  492  O2  U36 B  20      11.363  -2.445   7.904  1.00 42.96           O
HETATM  493  N3  U36 B  20      12.545  -0.747   8.850  1.00 43.77           N
HETATM  494  C4  U36 B  20      13.154   0.495   8.855  1.00 42.08           C
HETATM  495  O4  U36 B  20      13.738   0.884   9.875  1.00 37.89           O
HETATM  496  C5  U36 B  20      13.040   1.208   7.617  1.00 41.51           C
HETATM  497  C6  U36 B  20      12.374   0.651   6.587  1.00 39.50           C
HETATM  498  CA' U36 B  20       9.181  -3.182   5.789  1.00 49.64           C
HETATM  499  CB' U36 B  20       9.334  -4.295   4.806  1.00 50.46           C
HETATM  500  CD' U36 B  20      10.943  -5.943   4.104  1.00 51.43           C
HETATM  501  OC' U36 B  20      10.740  -4.603   4.641  1.00 52.78           O
""",
  "a44.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
A44        A44 '2'-O-(2-methoxyethyl)adenosine 5'-(dihydrogen phosphate)' ligand 45 27 .
#
data_comp_A44
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
A44         P      P   P     .          2.5224   -5.5637    0.6446
A44         OP2    O   OP    .          4.0244   -5.2909    0.6504
A44         OP1    O   O     .          2.0353   -5.7363    2.0663
A44        O5'     O   O2    .          1.7520   -4.3205   -0.0401
A44        C5'     C   CH2   .          1.7304   -3.1036    0.6336
A44        C4'     C   CR15  .          1.2062   -2.0228   -0.2918
A44        O4'     O   O     .         -0.0783   -1.9855   -0.2471
A44        C3'     C   CR15  .          1.6764   -0.6041    0.2038
A44        O3'     O   OH1   .          2.5794   -0.0871   -0.6532
A44        C2'     C   CR15  .          0.4371    0.2496    0.2212
A44        O2'     O   O2    .          0.6982    1.5749   -0.5050
A44        C1'     C   CR15  .         -0.4813   -0.4608   -0.4288
A44         N9     N   NR5   .         -1.7808   -0.2191    0.1349
A44         C8     C   CR15  .         -2.0670   -0.3505    1.4454
A44         N7     N   N     .         -3.4269   -0.2922    1.5891
A44         C5     C   CR56  .         -3.9797   -0.1222    0.3440
A44         C6     C   CR6   .         -5.3029   -0.0049   -0.0829
A44         N6     N   NH2   .         -6.5639    0.0390    0.5147
A44         N1     N   N     .         -5.5783    0.1542   -1.3948
A44         C2     C   CR16  .         -4.5689    0.2002   -2.3007
A44         N3     N   N     .         -3.2570    0.0844   -1.8806
A44         C4     C   CR56  .         -2.9767   -0.0766   -0.5532
A44        CA'     C   CH2   .          1.2611    2.6105    0.3430
A44        CD'     C   CH3   .          3.8911    5.0297   -0.3369
A44        OC'     O   O2    .          2.8946    4.3578    0.3441
A44        CB'     C   CH2   .          2.0410    3.5848   -0.4972
A44         OP3    O   OP    .          2.2422   -6.8280   -0.1418
""",
  "u36.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
U36        U36 '2'-O-(2-methoxyethyl)uridine 5'-(dihydrogen phosphate)' ligand 42 25 .
#
data_comp_U36
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
U36         P      P   P     .         -0.0923    3.8242    3.7169
U36         OP1    O   O     .         -0.7464    3.7757    5.0772
U36         OP2    O   OP    .         -1.0067    4.5301    2.7434
U36         OP3    O   OP    .          1.2153    4.5737    3.8136
U36        O5'     O   O2    .          0.1807    2.3510    3.2104
U36        C5'     C   CH2   .         -0.1828    2.0165    1.8914
U36        C4'     C   CH1   .          0.1665    0.7096    1.6366
U36        O4'     O   O2    .         -1.0654   -0.1293    1.4258
U36        C3'     C   CH1   .          1.0034    0.6117    0.2648
U36        O3'     O   OH1   .          2.4358    0.2844    0.5536
U36        C2'     C   CH1   .          0.4751   -0.3054   -0.3636
U36        O2'     O   O2    .          1.1615   -1.5865   -0.0670
U36        C1'     C   CH1   .         -1.0908   -0.3701    0.1523
U36         N1     N   NR6   .         -1.8533    0.5989   -0.5044
U36         C2     C   CR6   .         -2.9446    1.2709    0.2003
U36         O2     O   O     .         -3.2291    0.9400    1.3340
U36         N3     N   NR16  .         -3.6934    2.3410   -0.4578
U36         C4     C   CR6   .         -3.4325    2.6580   -1.8160
U36         O4     O   O     .         -4.1270    3.4693   -2.3943
U36         C5     C   CR16  .         -2.3699    1.9506   -2.5425
U36         C6     C   CR16  .         -1.5867    0.9104   -1.8626
U36        CA'     C   CH2   .          1.3872   -2.3894   -1.2028
U36        CB'     C   CH2   .          2.2425   -3.5784   -0.8184
U36        CD'     C   CH3   .          2.9132   -5.6983   -1.6208
U36        OC'     O   O2    .          1.9630   -4.6568   -1.6813
""",
  "linking_test_nstd_rna_dna.pdb" : """
HETATM  478  P   U36 B  20      11.477   3.696   3.004  1.00 42.30           P
HETATM  479  OP1 U36 B  20      10.335   4.004   2.145  1.00 44.51           O
HETATM  480  OP2 U36 B  20      11.842   4.626   4.093  1.00 45.81           O
HETATM  481  O5' U36 B  20      11.173   2.296   3.679  1.00 42.73           O
HETATM  482  C5' U36 B  20      10.757   1.210   2.877  1.00 42.45           C
HETATM  483  C4' U36 B  20      10.319   0.065   3.737  1.00 42.68           C
HETATM  484  O4' U36 B  20      11.457  -0.572   4.347  1.00 42.31           O
HETATM  485  C3' U36 B  20       9.406   0.368   4.909  1.00 41.73           C
HETATM  486  O3' U36 B  20       8.059   0.493   4.530  1.00 44.16           O
HETATM  487  C2' U36 B  20       9.567  -0.885   5.739  1.00 42.76           C
HETATM  488  O2' U36 B  20       8.836  -1.946   5.173  1.00 46.09           O
HETATM  489  C1' U36 B  20      11.046  -1.178   5.554  1.00 41.52           C
HETATM  490  N1  U36 B  20      11.789  -0.588   6.668  1.00 40.08           N
HETATM  491  C2  U36 B  20      11.863  -1.336   7.813  1.00 41.61           C
HETATM  492  O2  U36 B  20      11.363  -2.445   7.904  1.00 42.96           O
HETATM  493  N3  U36 B  20      12.545  -0.747   8.850  1.00 43.77           N
HETATM  494  C4  U36 B  20      13.154   0.495   8.855  1.00 42.08           C
HETATM  495  O4  U36 B  20      13.738   0.884   9.875  1.00 37.89           O
HETATM  496  C5  U36 B  20      13.040   1.208   7.617  1.00 41.51           C
HETATM  497  C6  U36 B  20      12.374   0.651   6.587  1.00 39.50           C
HETATM  498  CA' U36 B  20       9.181  -3.182   5.789  1.00 49.64           C
HETATM  499  CB' U36 B  20       9.334  -4.295   4.806  1.00 50.46           C
HETATM  500  CD' U36 B  20      10.943  -5.943   4.104  1.00 51.43           C
HETATM  501  OC' U36 B  20      10.740  -4.603   4.641  1.00 52.78           O
HET    C43  B  21      24
LINK         O3' U36 B  20                 P   C43 B  21     1555   1555  1.60
LINK         O3' C43 B  21                 P   G48 B  22     1555   1555  1.61
HETATM  502  P   C43 B  21       7.042   1.245   5.508  1.00 45.71           P
HETATM  503  OP2 C43 B  21       7.731   2.473   5.993  1.00 47.37           O
HETATM  504  OP1 C43 B  21       5.771   1.363   4.736  1.00 43.78           O
HETATM  505  O5' C43 B  21       6.853   0.247   6.738  1.00 39.90           O
HETATM  506  C5' C43 B  21       6.411  -1.086   6.501  1.00 38.58           C
HETATM  507  C4' C43 B  21       6.481  -1.910   7.754  1.00 37.96           C
HETATM  508  O4' C43 B  21       7.858  -2.119   8.141  1.00 37.13           O
HETATM  509  C3' C43 B  21       5.834  -1.356   9.001  1.00 39.94           C
HETATM  510  O3' C43 B  21       4.419  -1.513   9.009  1.00 38.74           O
HETATM  511  C2' C43 B  21       6.517  -2.174  10.095  1.00 42.40           C
HETATM  512  O2' C43 B  21       5.975  -3.463  10.290  1.00 51.49           O
HETATM  513  C1' C43 B  21       7.925  -2.350   9.532  1.00 38.54           C
HETATM  514  N1  C43 B  21       8.857  -1.402  10.159  1.00 36.39           N
HETATM  515  C2  C43 B  21       9.294  -1.667  11.451  1.00 33.63           C
HETATM  516  O2  C43 B  21       8.947  -2.703  11.998  1.00 28.79           O
HETATM  517  N3  C43 B  21      10.082  -0.780  12.075  1.00 34.10           N
HETATM  518  C4  C43 B  21      10.449   0.335  11.459  1.00 36.14           C
HETATM  519  N4  C43 B  21      11.208   1.196  12.143  1.00 39.09           N
HETATM  520  C5  C43 B  21      10.051   0.621  10.126  1.00 34.73           C
HETATM  521  C6  C43 B  21       9.265  -0.271   9.515  1.00 36.40           C
HETATM  522  CA' C43 B  21       6.672  -4.189  11.334  1.00 58.64           C
HETATM  523  CB' C43 B  21       7.917  -4.946  10.870  1.00 64.99           C
HETATM  524  CD' C43 B  21       8.788  -5.769  13.046  1.00 68.99           C
HETATM  525  OC' C43 B  21       9.059  -5.026  11.816  1.00 67.88           O
""",
  "c43.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
C43        C43 '2'-O-(2-methoxyethyl)cytidine 5'-(dihydrogen phosphate)' ligand 43 25 .
#
data_comp_C43
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
C43         P      P   P     .          3.7565    3.2915    1.5876
C43         OP3    O   OP    .          4.2694    2.5857    2.8218
C43         OP2    O   OP    .          2.7755    4.3524    1.9881
C43         OP1    O   O     .          4.9068    3.9231    0.8595
C43        O5'     O   O2    .          3.0486    2.2550    0.6486
C43        C5'     C   CH2   .          1.8203    1.6916    1.0686
C43        C4'     C   CR15  .          1.6741    0.4618    0.5078
C43        O4'     O   O     .          0.9711   -0.4760    1.4642
C43        C3'     C   CR15  .          0.7330    0.5479   -0.7879
C43        O3'     O   OH1   .          1.5421    0.3327   -2.0324
C43        C2'     C   CR15  .         -0.0850   -0.3670   -0.6482
C43        O2'     O   O2    .          0.4137   -1.6016   -1.3197
C43        C1'     C   CR15  .         -0.2173   -0.6200    0.9744
C43         N1     N   NR6   .         -1.0961    0.3267    1.5414
C43         C2     C   CR6   .         -2.2763    0.7649    0.7937
C43         O2     O   O     .         -2.5555    0.2318   -0.2825
C43         N3     N   N     .         -3.1751    1.7982    1.3634
C43         C4     C   CR6   .         -2.9944    2.2380    2.7088
C43         N4     N   NH2   .         -3.9728    3.1027    3.3229
C43         C5     C   CR16  .         -1.8514    1.7296    3.5031
C43         C6     C   CR16  .         -0.9017    0.7650    2.8947
C43        CA'     C   CH2   .         -0.4824   -2.1632   -2.1831
C43        CB'     C   CH2   .          0.1430   -3.3861   -2.8405
C43        CD'     C   CH3   .         -0.1994   -4.9988   -4.6237
C43        OC'     O   O2    .         -0.7247   -3.8624   -3.9069
""",
  "linking_test_Mg_HOH.pdb" : """
HETATM    1 MG   MG  A   1      14.481  28.862  20.807  1.00 20.00      A   MG+2
HETATM    2  O   HOH A   2      16.505  29.075  21.158  1.00 20.00      A    O
HETATM    3  O   HOH A   3      12.607  28.705  20.389  1.00 20.00      A    O
HETATM    4  O   HOH A   4      15.029  27.124  19.772  1.00 20.00      A    O
HETATM    5  O   HOH A   5      14.172  30.496  21.689  1.00 20.00      A    O
HETATM    6  O   HOH A   6      14.474  27.674  22.398  1.00 20.00      A    O
HETATM    7  O   HOH A   7      14.683  29.906  19.056  1.00 20.00      A    O
""",
  "linking_test_Mg_HOH_CRYST1.pdb" : """
CRYST1   13.898   13.372   13.342  90.00  90.00  90.00 P 1
SCALE1      0.071953  0.000000  0.000000        0.00000
SCALE2      0.000000  0.074783  0.000000        0.00000
SCALE3      0.000000  0.000000  0.074951        0.00000
HETATM    1 MG   MG  A   1       6.874   6.738   6.751  1.00 20.00      A   MG+2
HETATM    2  O   HOH A   2       8.898   6.951   7.102  1.00 20.00      A    O
HETATM    3  O   HOH A   3       5.000   6.581   6.333  1.00 20.00      A    O
HETATM    4  O   HOH A   4       7.422   5.000   5.716  1.00 20.00      A    O
HETATM    5  O   HOH A   5       6.565   8.372   7.633  1.00 20.00      A    O
HETATM    6  O   HOH A   6       6.867   5.550   8.342  1.00 20.00      A    O
HETATM    7  O   HOH A   7       7.076   7.782   5.000  1.00 20.00      A    O
TER
""",
  "linking_test_Mg_EDT.pdb" : """
HETATM    1 MG    MG A1501      -7.869  -0.167  32.075  0.76 19.71          Mg
HETATM    2  O01 EDT A   1      -8.014   0.505  29.989  0.49 19.08           O
HETATM    3  C02 EDT A   1      -8.838  -0.073  29.218  0.49 21.60           C
HETATM    4  O03 EDT A   1      -9.035   0.272  28.025  0.49 21.74           O
HETATM    5  C04 EDT A   1      -9.648  -1.250  29.726  0.49 20.37           C
HETATM    6  N05 EDT A   1      -9.634  -1.304  31.200  0.49 22.50           N
HETATM    7  C06 EDT A   1      -9.423  -2.672  31.724  0.49 20.66           C
HETATM    8  C07 EDT A   1      -7.955  -3.048  31.586  0.49 23.87           C
HETATM    9  O08 EDT A   1      -7.636  -4.245  31.436  0.49 23.73           O
HETATM   10  O09 EDT A   1      -7.092  -2.152  31.618  0.49 25.09           O
HETATM   11  C10 EDT A   1     -10.892  -0.735  31.713  0.49 21.64           C
HETATM   12  C11 EDT A   1     -10.806   0.778  31.984  0.49 21.39           C
HETATM   13  N12 EDT A   1      -9.522   1.212  32.554  0.49 22.50           N
HETATM   14  C13 EDT A   1      -9.643   1.195  34.028  0.49 20.27           C
HETATM   15  C14 EDT A   1      -8.884   0.036  34.628  0.49 21.97           C
HETATM   16  O15 EDT A   1      -9.056  -0.252  35.836  0.49 20.76           O
HETATM   17  O16 EDT A   1      -8.088  -0.608  33.907  0.49 19.08           O
HETATM   18  C17 EDT A   1      -9.242   2.584  32.090  0.49 20.69           C
HETATM   19  C18 EDT A   1      -7.769   2.967  32.220  0.49 23.96           C
HETATM   20  O19 EDT A   1      -6.885   2.088  32.322  0.49 25.12           O
HETATM   21  O20 EDT A   1      -7.479   4.176  32.238  0.49 23.80           O
""",
  "edt.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
EDT        EDT 'Unknown                  ' ligand 32 20 .
#
data_comp_EDT
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
EDT         O01    O   O     .         -0.9770   -5.1640  -35.4940
EDT         C02    C   C     .         -0.8990   -6.3850  -35.1670
EDT         O03    O   OC    .         -0.1220   -7.2090  -35.7080
EDT         C04    C   CH2   .         -1.7640   -6.9100  -34.0410
EDT         N05    N   NT    .         -1.9780   -5.8560  -33.0320
EDT         C06    C   CH2   .         -0.7560   -5.5740  -32.2430
EDT         C07    C   C     .         -0.5840   -4.0730  -32.0610
EDT         O08    O   O     .          0.3270   -3.6420  -31.3290
EDT         O09    O   OC    .         -1.3550   -3.2970  -32.6510
EDT         C10    C   CH2   .         -3.0920   -6.2330  -32.1430
EDT         C11    C   CH2   .         -4.0630   -5.0930  -31.8120
EDT         N12    N   NT    .         -4.3530   -4.1550  -32.9060
EDT         C13    C   CH2   .         -4.7520   -2.8590  -32.3240
EDT         C14    C   C     .         -4.2850   -1.7280  -33.2010
EDT         O15    O   O     .         -4.6870   -0.5650  -32.9710
EDT         O16    O   OC    .         -3.5120   -1.9950  -34.1490
EDT         C17    C   CH2   .         -5.4790   -4.6970  -33.6920
EDT         C18    C   C     .         -5.0620   -5.2790  -35.0410
EDT         O19    O   O     .         -3.9340   -5.0370  -35.5190
EDT         O20    O   OC    .         -5.8920   -5.9720  -35.6540
EDT        H041    H   HCH2  .         -1.3187   -7.6788  -33.6195
EDT        H042    H   HCH2  .         -2.6312   -7.1948  -34.4028
EDT        H061    H   HCH2  .          0.0263   -5.9362  -32.7143
EDT        H062    H   HCH2  .         -0.8318   -6.0020  -31.3616
EDT        H101    H   HCH2  .         -2.7149   -6.5735  -31.3006
EDT        H102    H   HCH2  .         -3.5994   -6.9574  -32.5721
EDT        H111    H   HCH2  .         -3.6877   -4.5805  -31.0618
EDT        H112    H   HCH2  .         -4.9130   -5.4915  -31.5177
EDT        H131    H   HCH2  .         -4.3504   -2.7658  -31.4322
EDT        H132    H   HCH2  .         -5.7312   -2.8259  -32.2447
EDT        H171    H   HCH2  .         -6.1277   -3.9755  -33.8490
EDT        H172    H   HCH2  .         -5.9159   -5.4037  -33.1675
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
EDT   C02     O01   deloc         1.266 0.020
EDT   O03     C02   deloc         1.255 0.020
EDT   C04     C02   single        1.514 0.020
EDT   N05     C04   single        1.475 0.020
EDT   C07     C06   single        1.522 0.020
EDT   O08     C07   deloc         1.246 0.020
EDT   O09     C07   deloc         1.243 0.020
EDT   C06     N05   single        1.482 0.020
EDT   C10     N05   single        1.474 0.020
EDT   C11     C10   single        1.534 0.020
EDT   N12     C11   single        1.470 0.020
EDT   C14     C13   single        1.505 0.020
EDT   O15     C14   deloc         1.252 0.020
EDT   O16     C14   deloc         1.252 0.020
EDT   C13     N12   single        1.476 0.020
EDT   C17     N12   single        1.476 0.020
EDT   C18     C17   single        1.527 0.020
EDT   O19     C18   deloc         1.249 0.020
EDT   O20     C18   deloc         1.243 0.020
EDT  H041     C04   single        0.983 0.020
EDT  H042     C04   single        0.982 0.020
EDT  H061     C06   single        0.982 0.020
EDT  H062     C06   single        0.983 0.020
EDT  H101     C10   single        0.984 0.020
EDT  H102     C10   single        0.983 0.020
EDT  H111     C11   single        0.983 0.020
EDT  H112     C11   single        0.984 0.020
EDT  H131     C13   single        0.983 0.020
EDT  H132     C13   single        0.983 0.020
EDT  H171     C17   single        0.983 0.020
EDT  H172     C17   single        0.983 0.020
#
""",
  "linking_test_1jbe_ALA-SNN-ACY-ALA.pdb" : """
ATOM    600  N   ALA A  74      28.223  10.182  26.006  1.00  9.40           N
ATOM    601  CA  ALA A  74      28.012  10.148  27.455  1.00 11.26           C
ATOM    602  C   ALA A  74      29.306  10.199  28.232  1.00 13.02           C
ATOM    603  O   ALA A  74      29.290  10.257  29.458  1.00 16.02           O
ATOM    604  CB  ALA A  74      27.262   8.877  27.829  1.00 15.69           C
HETATM  605  N1  SNN A  75      33.758  11.091  28.225  1.00 13.76           N
HETATM  606  C2  SNN A  75      32.448  11.310  28.512  1.00 12.79           C
HETATM  607  C3  SNN A  75      31.700  10.019  28.212  1.00 11.71           C
HETATM  608  N3  SNN A  75      30.432  10.162  27.536  1.00 11.34           N
HETATM  609  C4  SNN A  75      32.736   9.250  27.383  1.00 12.39           C
HETATM  610  C5  SNN A  75      34.123   9.826  27.531  1.00 12.20           C
HETATM  611  O2  SNN A  75      31.994  12.355  28.940  1.00 17.38           O
HETATM  612  O5  SNN A  75      35.186   9.317  27.190  1.00 17.87           O
HETATM  613  C   ACY A  76      35.504  12.683  27.565  1.00 16.21           C
HETATM  614  O   ACY A  76      34.812  13.387  26.829  1.00 17.55           O
HETATM  615  CH3 ACY A  76      34.852  11.912  28.684  1.00 16.85           C
ATOM    616  N   ALA A  77      36.810  12.555  27.393  1.00 17.53           N
ATOM    617  CA  ALA A  77      37.522  13.396  26.441  1.00 17.22           C
ATOM    618  C   ALA A  77      36.978  13.243  25.024  1.00 15.73           C
ATOM    619  O   ALA A  77      36.998  14.214  24.277  1.00 16.94           O
ATOM    620  CB  ALA A  77      39.018  13.091  26.481  1.00 27.39           C
""",
  "snn.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
SNN        SNN '(3S)-3-aminopyrrolidine-2,5-dione' ligand 14 8 .
#
data_comp_SNN
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
SNN         N1     N   NH1   .         -1.1318   -0.2275   -1.2446
SNN         C2     C   C     .          0.2247   -0.2275   -1.2446
SNN         C3     C   CH1   .          0.6936   -0.2275    0.2180
SNN         N3     N   NH2   .          1.6615    0.7434    0.3846
SNN         C4     C   CH2   .         -0.6212   -0.2275    1.0200
SNN         C5     C   C     .         -1.5673    0.3826   -0.0137
SNN         O2     O   O     .          0.9386   -0.3193   -2.2232
SNN         O5     O   O     .         -2.4854    1.1684    0.1080
SNN         HN     H   HNH1  .         -1.6677   -0.4697   -1.9417
SNN         H3     H   HCH1  .          1.2207   -0.9991    0.5003
SNN        HN31    H   HNH2  .          1.7785    0.9065    1.2742
SNN        HN32    H   HNH2  .          2.4472    0.4548    0.0225
SNN         H41    H   HCH2  .         -0.5948    0.1814    1.9060
SNN         H42    H   HCH2  .         -0.8965   -1.1391    1.2345
#
""",
  "linking_test_3gmq_NAG-FUC.pdb" : """
HETATM 3108  C1  NAG A 294       5.929 -23.329  20.979  1.00 36.07           C
HETATM 3109  C2  NAG A 294       6.923 -23.190  22.146  1.00 36.71           C
HETATM 3110  C3  NAG A 294       8.360 -23.145  21.617  1.00 36.96           C
HETATM 3111  C4  NAG A 294       8.692 -24.336  20.749  1.00 36.82           C
HETATM 3112  C5  NAG A 294       7.591 -24.460  19.696  1.00 38.98           C
HETATM 3113  C6  NAG A 294       7.798 -25.785  18.952  1.00 38.89           C
HETATM 3114  C7  NAG A 294       5.795 -21.923  23.855  1.00 34.86           C
HETATM 3115  C8  NAG A 294       5.669 -20.618  24.563  1.00 36.07           C
HETATM 3116  N2  NAG A 294       6.712 -21.966  22.871  1.00 35.60           N
HETATM 3117  O3  NAG A 294       9.221 -23.104  22.736  1.00 39.69           O
HETATM 3118  O4  NAG A 294       9.926 -24.109  20.077  1.00 36.71           O
HETATM 3119  O5  NAG A 294       6.288 -24.471  20.285  1.00 34.73           O
HETATM 3120  O6  NAG A 294       6.889 -25.880  17.874  1.00 42.57           O
HETATM 3121  O7  NAG A 294       5.008 -22.854  24.057  1.00 40.29           O
HETATM 3169  C1  FUC A 299       7.292 -25.423  16.517  1.00 46.98           C
HETATM 3170  C2  FUC A 299       6.227 -26.090  15.616  1.00 49.47           C
HETATM 3171  C3  FUC A 299       4.835 -25.506  15.941  1.00 49.89           C
HETATM 3172  C4  FUC A 299       4.818 -23.991  15.655  1.00 53.30           C
HETATM 3173  C5  FUC A 299       5.984 -23.378  16.459  1.00 50.98           C
HETATM 3174  C6  FUC A 299       6.038 -21.848  16.375  1.00 45.12           C
HETATM 3175  O2  FUC A 299       6.203 -27.506  15.815  1.00 49.49           O
HETATM 3176  O3  FUC A 299       3.810 -26.238  15.266  1.00 47.77           O
HETATM 3177  O4  FUC A 299       4.964 -23.690  14.277  1.00 51.21           O
HETATM 3178  O5  FUC A 299       7.240 -24.036  16.163  1.00 46.78           O
""",
  "linking_test_CD_GHE_A_B.pdb" : """
CRYST1  154.379  108.579  120.263  90.00 138.53  90.00 C 1 2 1      20
HETATM 8577 CD    CD C1205     -85.285  28.175  80.928  1.00 13.07          CD
HETATM 8578 CD    CD C1206     -86.132  30.445  78.070  1.00 13.45          CD
HETATM 8579  N1 AGHE C1207     -87.632  33.273  82.068  0.25 15.93           N
HETATM 8580  C2 AGHE C1207     -87.139  31.954  82.363  0.25 13.77           C
HETATM 8581  C3 AGHE C1207     -85.848  32.181  82.970  0.25 13.59           C
HETATM 8582  C4 AGHE C1207     -85.379  33.497  82.427  0.25 15.16           C
HETATM 8583  C5 AGHE C1207     -86.627  34.267  82.225  0.25 18.36           C
HETATM 8584  C1 AGHE C1207     -86.965  31.040  81.200  0.25 12.73           C
HETATM 8585  O2 AGHE C1207     -86.877  31.463  79.997  0.25 13.10           O
HETATM 8586  O1 AGHE C1207     -86.879  29.793  81.461  0.25 11.38           O
HETATM 8587  C6 AGHE C1207     -88.988  33.572  81.719  0.25 20.72           C
HETATM 8588  O6 AGHE C1207     -89.317  34.719  81.513  0.25 24.50           O
HETATM 8589  C7 AGHE C1207     -89.928  32.396  81.589  0.25 13.71           C
HETATM 8590  C8 AGHE C1207     -90.599  32.349  80.226  0.25 23.85           C
HETATM 8591  C9 AGHE C1207     -89.628  32.441  79.013  0.25 23.81           C
HETATM 8592  C10AGHE C1207     -90.298  32.455  77.696  0.25 13.75           C
HETATM 8593  C11AGHE C1207     -91.284  33.597  77.576  0.25 20.78           C
HETATM 8594  O11AGHE C1207     -90.998  34.811  77.764  0.25 24.27           O
HETATM 8595  N2 AGHE C1207     -92.616  33.250  77.274  0.25 15.91           N
HETATM 8596  C15AGHE C1207     -93.060  31.918  77.024  0.25 13.70           C
HETATM 8597  C14AGHE C1207     -94.339  32.116  76.287  0.25 13.63           C
HETATM 8598  C13AGHE C1207     -94.859  33.436  76.736  0.25 15.20           C
HETATM 8599  C12AGHE C1207     -93.641  34.215  77.139  0.25 18.35           C
HETATM 8600  C16AGHE C1207     -93.276  31.124  78.251  0.25 12.59           C
HETATM 8601  O16AGHE C1207     -93.343  29.864  78.193  0.25 10.81           O
HETATM 8602  O17AGHE C1207     -93.406  31.667  79.357  0.25 13.15           O
HETATM 8579  N1 BGHE C1207     -87.762  33.127  82.051  0.25 15.93           N
HETATM 8580  C2 BGHE C1207     -87.157  31.870  82.281  0.25 13.77           C
HETATM 8581  C3 BGHE C1207     -85.824  32.206  82.976  0.25 13.59           C
HETATM 8582  C4 BGHE C1207     -85.481  33.523  82.395  0.25 15.16           C
HETATM 8583  C5 BGHE C1207     -86.817  34.170  82.230  0.25 18.36           C
HETATM 8584  C1 BGHE C1207     -86.997  31.025  81.014  0.25 12.73           C
HETATM 8585  O2 BGHE C1207     -86.879  31.590  79.969  0.25 13.10           O
HETATM 8586  O1 BGHE C1207     -86.924  29.833  81.336  0.25 11.38           O
HETATM 8587  C6 BGHE C1207     -89.137  33.243  81.727  0.25 20.72           C
HETATM 8588  O6 BGHE C1207     -89.853  32.297  81.625  0.25 24.50           O
HETATM 8589  C7 BGHE C1207     -89.559  34.637  81.532  0.25 13.71           C
HETATM 8590  C8 BGHE C1207     -90.563  34.607  80.287  0.25 23.85           C
HETATM 8591  C9 BGHE C1207     -89.667  34.602  79.002  0.25 23.81           C
HETATM 8592  C10BGHE C1207     -90.628  34.667  77.787  0.25 13.75           C
HETATM 8593  C11BGHE C1207     -91.092  33.238  77.563  0.25 20.78           C
HETATM 8594  O11BGHE C1207     -90.376  32.288  77.659  0.25 24.27           O
HETATM 8595  N2 BGHE C1207     -92.458  33.129  77.232  0.25 15.91           N
HETATM 8596  C15BGHE C1207     -93.058  31.883  76.995  0.25 13.70           C
HETATM 8597  C14BGHE C1207     -94.408  32.218  76.304  0.25 13.63           C
HETATM 8598  C13BGHE C1207     -94.743  33.524  76.887  0.25 15.20           C
HETATM 8599  C12BGHE C1207     -93.408  34.169  77.056  0.25 18.35           C
HETATM 8600  C16BGHE C1207     -93.228  31.034  78.271  0.25 12.59           C
HETATM 8601  O16BGHE C1207     -93.298  29.819  77.942  0.25 10.81           O
HETATM 8602  O17BGHE C1207     -93.349  31.594  79.315  0.25 13.15           O
""",
  "ghe.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
GHE        GHE 'Unknown                  ' ligand 46 24 .
#
data_comp_GHE
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
GHE         N1     N   NT    .        -79.2922   87.5583   63.7529
GHE         C2     C   CH1   .        -78.2204   88.2930   63.9823
GHE         C3     C   CH2   .        -77.1377   87.8430   62.8003
GHE         C4     C   CH2   .        -77.3209   86.5460   62.6252
GHE         C5     C   CH2   .        -78.8231   86.2738   63.0451
GHE         C1     C   C     .        -77.6682   88.1432   65.4087
GHE         O2     O   O     .        -77.7245   87.0287   65.9917
GHE         O1     O   OC    .        -77.1497   89.1322   65.9901
GHE         C6     C   C     .        -80.6286   87.8903   64.0949
GHE         O6     O   O     .        -80.8690   88.9411   64.6552
GHE         C7     C   CH2   .        -81.7587   86.9547   63.7250
GHE         C8     C   CH2   .        -83.1158   87.5330   64.0610
GHE         C9     C   CH2   .        -84.1621   87.1799   63.0266
GHE         C10    C   CH2   .        -85.4675   87.9082   63.2593
GHE         C11    C   C     .        -86.3747   87.8603   62.0495
GHE         O11    O   O     .        -86.0710   87.1812   61.0892
GHE         N2     N   NT    .        -87.5618   88.6367   62.0221
GHE         C15    C   CH1   .        -88.4134   88.7124   61.0336
GHE         C14    C   CH2   .        -89.0906   90.1723   61.1953
GHE         C13    C   CH2   .        -88.6543   90.6746   62.5810
GHE         C12    C   CH2   .        -88.1072   89.5871   63.2049
GHE         C16    C   C     .        -89.4518   87.5401   61.0112
GHE         O16    O   O     .        -89.6653   86.9087   59.9431
GHE         O17    O   OC    .        -90.0884   87.2431   62.0560
GHE         H21    H   HCH1  .        -78.4208   89.2507   63.9474
GHE         H31    H   HCH2  .        -77.3239   88.3266   61.9687
GHE         H32    H   HCH2  .        -76.2204   88.0184   63.0967
GHE         H41    H   HCH2  .        -77.1829   86.3077   61.6850
GHE         H42    H   HCH2  .        -76.7108   86.0367   63.1983
GHE         H51    H   HCH2  .        -78.8705   85.5130   63.6581
GHE         H52    H   HCH2  .        -79.3707   86.1051   62.2524
GHE         H71    H   HCH2  .        -81.7199   86.7754   62.7642
GHE         H72    H   HCH2  .        -81.6444   86.1127   64.2096
GHE         H81    H   HCH2  .        -83.0404   88.5067   64.1157
GHE         H82    H   HCH2  .        -83.4008   87.1876   64.9306
GHE         H91    H   HCH2  .        -84.3269   86.2163   63.0589
GHE         H92    H   HCH2  .        -83.8224   87.4141   62.1397
GHE        H101    H   HCH2  .        -85.2749   88.8429   63.4740
GHE        H102    H   HCH2  .        -85.9288   87.4971   64.0175
GHE        H151    H   HCH1  .        -87.9759   88.6068   60.1655
GHE        H141    H   HCH2  .        -88.7644   90.7760   60.4982
GHE        H142    H   HCH2  .        -90.0651   90.1000   61.1528
GHE        H131    H   HCH2  .        -87.9860   91.3852   62.4895
GHE        H132    H   HCH2  .        -89.4297   90.9997   63.0837
GHE        H121    H   HCH2  .        -88.7897   89.1182   63.7284
GHE        H122    H   HCH2  .        -87.3675   89.8655   63.7838
""",
  "linking_test_ALY_MCM.pdb" : """
HETATM 2888  OH  ALY C   5     -15.286  24.739   8.076  1.00 45.73      C    O
HETATM 2889  CH  ALY C   5     -16.455  25.116   7.877  1.00 41.87      C    C
HETATM 2890  CH3 ALY C   5     -16.828  25.980   6.695  1.00 36.75      C    C
HETATM 2891  NZ  ALY C   5     -17.470  24.749   8.740  1.00 48.84      C    N
HETATM 2892  CE  ALY C   5     -17.263  23.930   9.901  1.00 39.21      C    C
HETATM 2893  CD  ALY C   5     -18.562  23.629  10.640  1.00 41.37      C    C
HETATM 2894  CG  ALY C   5     -18.329  22.814  11.892  1.00 41.99      C    C
HETATM 2895  CB  ALY C   5     -19.662  22.656  12.590  1.00 45.24      C    C
HETATM 2896  CA  ALY C   5     -19.524  22.011  13.974  1.00 51.20      C    C
HETATM 2897  N   ALY C   5     -20.755  21.670  14.665  1.00 49.87      C    N
HETATM 2898  C   ALY C   5     -18.619  22.936  14.781  1.00 47.11      C    C
HETATM 2899  O   ALY C   5     -17.427  22.728  15.024  1.00 49.30      C    O
HETATM 2900  N   MCM C   6     -19.324  23.360  15.960  1.00 49.19      C    N
HETATM 2901  CA  MCM C   6     -18.698  24.321  16.695  1.00 44.46      C    C
HETATM 2902  C2  MCM C   6     -19.284  25.579  16.752  1.00 49.13      C    C
HETATM 2903  C3  MCM C   6     -18.708  26.580  17.525  1.00 48.36      C    C
HETATM 2904  C4  MCM C   6     -18.700  28.830  18.338  1.00 45.36      C    C
HETATM 2905  C5  MCM C   6     -17.544  28.627  19.099  1.00 42.82      C    C
HETATM 2906  C6  MCM C   6     -16.902  27.390  19.114  1.00 45.52      C    C
HETATM 2907  C7  MCM C   6     -17.484  26.294  18.302  1.00 46.42      C    C
HETATM 2908  C8  MCM C   6     -16.929  25.030  18.246  1.00 40.92      C    C
HETATM 2909  C9  MCM C   6     -17.552  24.060  17.457  1.00 43.70      C    C
HETATM 2910  C10 MCM C   6     -15.654  27.180  19.937  1.00 42.32      C    C
HETATM 2911  O1  MCM C   6     -19.260  27.832  17.578  1.00 46.61      C    O
HETATM 2912  O2  MCM C   6     -19.248  29.950  18.344  1.00 44.55      C    O
""",
  "mcm.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
MCM        MCM '7-amino-4-methyl-2H-chromen-2-one' ligand 22 13 .
#
data_comp_MCM
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
MCM         N      N   NH2   .         -3.8313    0.0126   -0.1203
MCM         CA     C   CR6   .         -2.4124    0.0126   -0.1203
MCM         C2     C   CR16  .         -1.7183    0.0126    1.0827
MCM         C3     C   CR66  .         -0.3102    0.0126    1.0833
MCM         C4     C   CR6   .          1.7583    0.0126    2.2891
MCM         C5     C   CR16  .          2.4557    0.0126    1.0847
MCM         C6     C   CR6   .          1.7703    0.0126   -0.1062
MCM         C7     C   CR66  .          0.3770    0.0126   -0.1073
MCM         C8     C   CR16  .         -0.3279    0.0126   -1.3263
MCM         C9     C   CR16  .         -1.7167    0.0126   -1.3257
MCM         C10    C   CH3   .          2.5277   -0.0477   -1.4146
MCM         O1     O   O     .          0.3860    0.0126    2.2882
MCM         O2     O   OC    .          2.3647    0.0146    3.3417
MCM         HN1    H   HNH2  .         -4.2873    0.8024   -0.1159
MCM         HN2    H   HNH2  .         -4.2873   -0.7772   -0.1160
MCM         H2     H   HCR6  .         -2.2066    0.0126    1.9279
MCM         H5     H   HCR6  .          3.4319    0.0126    1.0860
MCM         H8     H   HCR6  .          0.1597    0.0126   -2.1720
MCM         H9     H   HCR6  .         -2.2050    0.0131   -2.1709
MCM        H101    H   HCH3  .          2.3139    0.7394   -1.9509
MCM        H102    H   HCH3  .          3.4869   -0.0680   -1.2348
MCM        H103    H   HCH3  .          2.2709   -0.8533   -1.9024
""",
  "linking_test_XYP_XYP.pdb" : """
HETATM 2826  O3B XYP E   1     -19.496  30.622 -25.625  1.00 40.53           O
HETATM 2827  C3B XYP E   1     -20.078  29.454 -25.182  1.00 38.50           C
HETATM 2828  C4B XYP E   1     -19.553  28.253 -25.964  1.00 42.57           C
HETATM 2829  O4B XYP E   1     -20.091  28.171 -27.284  1.00 43.82           O
HETATM 2830  C5B XYP E   1     -19.957  27.023 -25.162  1.00 43.75           C
HETATM 2831  O5B XYP E   1     -19.354  27.144 -23.900  1.00 41.03           O
HETATM 2832  C1B XYP E   1     -19.814  28.095 -22.997  1.00 48.42           C
HETATM 2833  C2B XYP E   1     -19.795  29.449 -23.683  1.00 38.61           C
HETATM 2834  O2B XYP E   1     -20.419  30.487 -22.970  1.00 36.77           O
HETATM 2835  O3B XYP E   2     -17.702  26.306 -20.601  1.00 21.95           O
HETATM 2836  C3B XYP E   2     -18.690  27.191 -20.166  1.00 22.83           C
HETATM 2837  C4B XYP E   2     -19.848  27.282 -21.144  1.00 29.33           C
HETATM 2838  O4B XYP E   2     -19.535  27.935 -22.359  1.00 39.74           O
HETATM 2839  C5B XYP E   2     -21.016  27.978 -20.518  1.00 32.56           C
HETATM 2840  O5B XYP E   2     -21.364  27.415 -19.297  1.00 36.29           O
HETATM 2841  C1B XYP E   2     -20.357  27.458 -18.342  1.00 36.98           C
HETATM 2842  C2B XYP E   2     -19.138  26.714 -18.808  1.00 27.91           C
HETATM 2843  O2B XYP E   2     -18.170  26.657 -17.835  1.00 29.93           O
""",
  "xyp.cif" : """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
XYP        XYP 'beta-D-xylopyranose      ' ligand 20 10 .
#
data_comp_XYP
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
XYP         O4A    O   OH1   .         -2.6653    0.5410   -0.7788
XYP         C1B    C   CH1   .         -1.2854    0.5410   -0.7788
XYP         C2B    C   CH1   .         -0.7949    0.5410    0.6557
XYP         C3B    C   CH1   .          0.6771    0.3610    0.7492
XYP         C4B    C   CH1   .          1.1508   -0.8250   -0.0209
XYP         C5B    C   CH2   .          0.6485   -0.8204   -1.4579
XYP         O2B    O   OH1   .         -1.1499    1.7930    1.2682
XYP         O3B    O   OH1   .          1.0411    0.1907    2.1434
XYP         O4B    O   OH1   .          2.5828   -0.8283   -0.0313
XYP         O5B    O   O2    .         -0.8048   -0.6378   -1.5443
XYP        HO4A    H   HOH1  .         -2.9613    0.7962   -1.5757
XYP         H1B    H   HCH1  .         -0.9700    1.3528   -1.2196
XYP         H2B    H   HCH1  .         -1.2337   -0.1841    1.1400
XYP         H3B    H   HCH1  .          1.1168    1.1604    0.4021
XYP         H4B    H   HCH1  .          0.8356   -1.6360    0.4216
XYP        H5B1    H   HCH2  .          0.8851   -1.6694   -1.8775
XYP        H5B2    H   HCH2  .          1.0868   -0.0955   -1.9429
XYP        HO2B    H   HOH1  .         -1.7421    1.6494    1.9136
XYP        HO3B    H   HOH1  .          0.7049   -0.5731    2.4459
XYP        HO4B    H   HOH1  .          2.8778   -1.6570    0.0880
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
XYP   O4A     C1B   single        1.359 0.020     1.359
XYP   O4A    HO4A   single        0.850 0.020     0.980
XYP   C1B     C2B   single        1.528 0.020     1.528
XYP   C1B     O5B   single        1.500 0.020     1.500
XYP   C1B     H1B   single        0.970 0.020     1.090
XYP   C2B     C3B   single        1.489 0.020     1.489
XYP   C2B     O2B   single        1.435 0.020     1.435
XYP   C2B     H2B   single        0.970 0.020     1.090
XYP   C3B     C4B   single        1.495 0.020     1.495
XYP   C3B     O3B   single        1.453 0.020     1.453
XYP   C3B     H3B   single        0.970 0.020     1.090
XYP   C4B     C5B   single        1.536 0.020     1.536
XYP   C4B     O4B   single        1.428 0.020     1.428
XYP   C4B     H4B   single        0.970 0.020     1.090
XYP   C5B     O5B   single        1.479 0.020     1.479
XYP   C5B    H5B1   single        0.970 0.020     1.090
XYP   C5B    H5B2   single        0.970 0.020     1.090
XYP   O2B    HO2B   single        0.850 0.020     0.980
XYP   O3B    HO3B   single        0.850 0.020     0.980
XYP   O4B    HO4B   single        0.850 0.020     0.980
""",
# CRYST1   16.302   21.238   70.013  90.07  90.04  67.31 P 1
# SCALE1      0.061342 -0.025647  0.000016        0.00000
# SCALE2      0.000000  0.051035  0.000053        0.00000
# SCALE3      0.000000  0.000000  0.014283        0.00000
  'linking_test_cyclic_main_chain.pdb' : '''
ATOM      1  N   SER A   1      -4.265  -4.954  -1.607  1.00  8.10           N
ATOM      2  CA  SER A   1      -4.144  -5.756  -2.800  1.00 12.96           C
ATOM      3  C   SER A   1      -4.846  -5.023  -3.931  1.00 15.20           C
ATOM      4  O   SER A   1      -5.852  -4.345  -3.752  1.00 13.88           O
ATOM      5  CB  SER A   1      -4.783  -7.118  -2.597  1.00 11.98           C
ATOM      6  OG  SER A   1      -6.157  -6.982  -2.308  1.00 20.63           O
ATOM      7  H1  SER A   1      -3.830  -5.404  -0.877  1.00  9.73           H
ATOM      8  H3  SER A   1      -5.194  -4.834  -1.375  1.00  9.73           H
ATOM      9  HA  SER A   1      -3.212  -5.906  -3.021  1.00 15.57           H
ATOM     10  HB2 SER A   1      -4.679  -7.641  -3.407  1.00 14.39           H
ATOM     11  HB3 SER A   1      -4.346  -7.566  -1.855  1.00 14.39           H
ATOM     12  HG  SER A   1      -6.551  -6.619  -2.955  1.00 24.78           H
ATOM     13  N   GLY A   2      -4.267  -5.151  -5.101  1.00  9.33           N
ATOM     14  CA  GLY A   2      -4.825  -4.626  -6.311  1.00 12.49           C
ATOM     15  C   GLY A   2      -5.008  -5.783  -7.263  1.00  9.40           C
ATOM     16  O   GLY A   2      -4.714  -6.935  -6.932  1.00 11.00           O
ATOM     17  H   GLY A   2      -3.518  -5.556  -5.222  1.00 11.21           H
ATOM     18  HA2 GLY A   2      -5.683  -4.208  -6.136  1.00 15.01           H
ATOM     19  HA3 GLY A   2      -4.229  -3.970  -6.704  1.00 15.01           H
ATOM     20  N   CYS A  33      -2.546   1.149   1.766  1.00  4.81           N
ATOM     21  CA  CYS A  33      -2.165   0.125   0.812  1.00  7.98           C
ATOM     22  C   CYS A  33      -3.153  -1.036   0.819  1.00  7.86           C
ATOM     23  O   CYS A  33      -3.760  -1.377   1.846  1.00  9.25           O
ATOM     24  CB  CYS A  33      -0.764  -0.397   1.128  1.00  8.49           C
ATOM     25  SG  CYS A  33       0.515   0.890   1.039  1.00  6.44           S
ATOM     26  H   CYS A  33      -2.063   1.185   2.476  1.00  5.78           H
ATOM     27  HA  CYS A  33      -2.161   0.511  -0.078  1.00  9.59           H
ATOM     28  HB2 CYS A  33      -0.760  -0.760   2.027  1.00 10.20           H
ATOM     29  HB3 CYS A  33      -0.536  -1.090   0.489  1.00 10.20           H
ATOM     30  N   GLY A  34      -3.323  -1.639  -0.350  1.00  6.42           N
ATOM     31  CA  GLY A  34      -4.109  -2.847  -0.461  1.00  8.51           C
ATOM     32  C   GLY A  34      -3.699  -3.646  -1.668  1.00  5.95           C
ATOM     33  O   GLY A  34      -2.828  -3.223  -2.415  1.00  5.28           O
ATOM     34  H   GLY A  34      -2.991  -1.363  -1.093  1.00  7.72           H
ATOM     35  HA2 GLY A  34      -3.984  -3.392   0.331  1.00 10.22           H
ATOM     36  HA3 GLY A  34      -5.048  -2.619  -0.543  1.00 10.22           H
''',
  "linking_test_cyclic_side_chain.pdb" : """
ATOM   3269  N   SER P   1      18.059  15.260 -24.071  1.00 30.30           N
ATOM   3270  CA  SER P   1      17.744  15.799 -25.407  1.00 28.03           C
ATOM   3271  C   SER P   1      17.804  14.711 -26.452  1.00 26.79           C
ATOM   3272  O   SER P   1      17.637  14.921 -27.647  1.00 26.96           O
ATOM   3273  CB  SER P   1      18.762  16.955 -25.664  1.00 24.81           C
ATOM   3274  OG  SER P   1      20.066  16.481 -25.290  1.00 25.91           O
ATOM   3275  N   HIS P   2      18.186  13.483 -26.051  1.00 27.55           N
ATOM   3276  CA  HIS P   2      18.293  12.344 -26.970  1.00 26.35           C
ATOM   3277  C   HIS P   2      17.629  11.119 -26.311  1.00 28.18           C
ATOM   3278  O   HIS P   2      17.571  11.054 -25.068  1.00 26.85           O
ATOM   3279  CB  HIS P   2      19.701  11.995 -27.457  1.00 25.52           C
ATOM   3280  CG  HIS P   2      20.361  13.249 -28.028  1.00 28.42           C
ATOM   3281  ND1 HIS P   2      20.863  14.172 -27.147  1.00 24.54           N
ATOM   3282  CD2 HIS P   2      20.448  13.734 -29.305  1.00 26.62           C
ATOM   3283  CE1 HIS P   2      21.328  15.194 -27.878  1.00 26.65           C
ATOM   3284  NE2 HIS P   2      21.104  14.954 -29.164  1.00 26.82           N
ATOM   3285  N   PHE P   3      17.144  10.201 -27.143  1.00 29.14           N
ATOM   3286  CA  PHE P   3      16.344   9.075 -26.655  1.00 30.32           C
ATOM   3287  C   PHE P   3      17.184   8.115 -25.808  1.00 30.69           C
ATOM   3288  O   PHE P   3      16.528   7.489 -24.951  1.00 28.96           O
ATOM   3289  CB  PHE P   3      15.661   8.283 -27.795  1.00 30.32           C
ATOM   3290  CG  PHE P   3      16.571   7.318 -28.492  1.00 31.30           C
ATOM   3291  CD1 PHE P   3      17.474   7.781 -29.449  1.00 31.02           C
ATOM   3292  CD2 PHE P   3      16.540   5.956 -28.197  1.00 32.14           C
ATOM   3293  CE1 PHE P   3      18.367   6.926 -30.034  1.00 32.25           C
ATOM   3294  CE2 PHE P   3      17.438   5.082 -28.804  1.00 31.93           C
ATOM   3295  CZ  PHE P   3      18.354   5.566 -29.735  1.00 32.42           C
ATOM   3296  N   ASN P   4      18.483   8.081 -26.027  1.00 29.01           N
ATOM   3297  CA  ASN P   4      19.372   7.142 -25.362  1.00 30.69           C
ATOM   3298  C   ASN P   4      20.157   7.719 -24.189  1.00 28.07           C
ATOM   3299  O   ASN P   4      21.300   7.375 -23.910  1.00 30.39           O
ATOM   3300  CB  ASN P   4      20.314   6.536 -26.417  1.00 30.85           C
ATOM   3301  CG  ASN P   4      21.272   7.509 -27.112  1.00 30.81           C
ATOM   3302  OD1 ASN P   4      20.921   8.699 -27.117  1.00 29.68           O
ATOM   3303  ND2 ASN P   4      22.387   7.033 -27.684  1.00 24.16           N
ATOM   3304  N   GLU P   5      19.570   8.666 -23.470  1.00 29.50           N
ATOM   3305  CA  GLU P   5      20.265   9.227 -22.291  1.00 29.41           C
ATOM   3306  C   GLU P   5      19.129   9.685 -21.395  1.00 29.39           C
ATOM   3307  O   GLU P   5      18.001   9.896 -21.859  1.00 29.14           O
ATOM   3308  CB  GLU P   5      21.196  10.384 -22.704  1.00 26.92           C
ATOM   3309  CG  GLU P   5      20.427  11.613 -23.230  1.00 25.95           C
ATOM   3310  CD  GLU P   5      21.249  12.501 -24.174  1.00 24.67           C
ATOM   3311  OE1 GLU P   5      22.366  12.144 -24.576  1.00 22.04           O
ATOM   3312  OE2 GLU P   5      20.781  13.651 -24.433  1.00 24.57           O
ATOM   3313  N   TYR P   6      19.464   9.979 -20.145  1.00 31.01           N
ATOM   3314  CA  TYR P   6      18.409  10.410 -19.228  1.00 34.43           C
ATOM   3315  C   TYR P   6      18.354  11.920 -19.206  1.00 33.94           C
ATOM   3316  O   TYR P   6      19.419  12.506 -18.972  1.00 35.68           O
ATOM   3317  CB  TYR P   6      18.737   9.834 -17.839  1.00 34.87           C
ATOM   3318  CG  TYR P   6      17.785  10.332 -16.768  1.00 36.73           C
ATOM   3319  CD1 TYR P   6      16.458   9.895 -16.735  1.00 36.01           C
ATOM   3320  CD2 TYR P   6      18.234  11.220 -15.814  1.00 35.99           C
ATOM   3321  CE1 TYR P   6      15.596  10.343 -15.743  1.00 35.81           C
ATOM   3322  CE2 TYR P   6      17.367  11.686 -14.813  1.00 36.01           C
ATOM   3323  CZ  TYR P   6      16.070  11.225 -14.792  1.00 37.50           C
ATOM   3324  OH  TYR P   6      15.232  11.666 -13.801  1.00 37.68           O
ATOM   3325  N   GLU P   7      17.200  12.518 -19.441  1.00 33.03           N
ATOM   3326  CA  GLU P   7      17.093  13.964 -19.391  1.00 38.47           C
ATOM   3327  C   GLU P   7      16.164  14.469 -18.288  1.00 44.24           C
ATOM   3328  O   GLU P   7      16.286  15.663 -17.967  1.00 46.55           O
ATOM   3329  CB  GLU P   7      16.582  14.559 -20.709  1.00 34.66           C
ATOM   3330  CG  GLU P   7      17.504  14.149 -21.826  1.00 31.31           C
ATOM   3331  CD  GLU P   7      17.202  14.901 -23.139  1.00 32.19           C
ATOM   3332  OE1 GLU P   7      15.956  14.941 -23.401  1.00 29.65           O
ATOM   3333  OXT GLU P   7      15.672  13.468 -17.541  1.00 45.81           O
""",
  "linking_test_over_valence.pdb" : """
ATOM   3451  C27 SEI L   1      65.362  58.026  53.718  1.00 20.00           C
ATOM   3463  O39 SEI L   1      65.091  58.948  52.663  1.00 20.00           O
ATOM   3495  H39 SEI L   1      65.115  58.519  51.936  1.00 20.00           H
ATOM   2595  N   GLY A 193      65.965  60.594  51.166  1.00  2.00           N
ATOM   2596  CA  GLY A 193      64.824  61.465  50.949  1.00  2.00           C
ATOM   2597  C   GLY A 193      63.671  60.712  50.330  1.00  3.58           C
ATOM   2598  O   GLY A 193      62.604  61.282  50.115  1.00  6.18           O
ATOM   2599  H   GLY A 193      65.893  60.083  51.841  1.00  2.00           H
ATOM   2600  HA2 GLY A 193      64.531  61.829  51.799  1.00  2.00           H
ATOM   2601  HA3 GLY A 193      65.062  62.202  50.366  1.00  2.00           H
""",
  "sei.cif" : """
data_comp_SEI
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.charge
SEI C27    C 'CH1 '  0.000 0
SEI O39    O 'OH1 '  0.000 0
SEI H39    H 'HCH1'  0.000 0
""",
  "linking_test_c2_c6.pdb" : """
HETATM15124  C1  BMA B1954      12.904 -11.942 -55.687  1.00 47.00           C
HETATM15125  C2  BMA B1954      12.394 -13.167 -54.934  1.00 46.45           C
HETATM15126  C3  BMA B1954      13.532 -14.149 -54.722  1.00 47.94           C
HETATM15127  C4  BMA B1954      14.248 -14.423 -56.039  1.00 48.98           C
HETATM15128  C5  BMA B1954      14.607 -13.131 -56.783  1.00 50.00           C
HETATM15129  C6  BMA B1954      15.208 -13.371 -58.158  1.00 49.36           C
HETATM15130  O2  BMA B1954      11.341 -13.818 -55.637  1.00 43.37           O
HETATM15131  O3  BMA B1954      12.997 -15.377 -54.260  1.00 50.40           O
HETATM15132  O4  BMA B1954      15.420 -15.170 -55.759  1.00 49.22           O
HETATM15133  O5  BMA B1954      13.436 -12.314 -56.945  1.00 47.68           O
HETATM15134  O6  BMA B1954      14.350 -14.201 -58.934  1.00 50.59           O
HETATM15157  C1  MAN B1957      13.587 -14.646 -60.954  1.00 26.90           C
HETATM15158  C2  MAN B1957      14.671 -14.846 -59.882  1.00 25.95           C
HETATM15159  C3  MAN B1957      15.018 -16.290 -60.094  1.00 25.04           C
HETATM15160  C4  MAN B1957      13.866 -17.159 -59.640  1.00 24.45           C
HETATM15161  C5  MAN B1957      12.714 -16.962 -60.596  1.00 24.51           C
HETATM15162  C6  MAN B1957      11.427 -17.605 -60.044  1.00 26.73           C
HETATM15163  O3  MAN B1957      16.232 -16.635 -59.360  1.00 26.05           O
HETATM15164  O4  MAN B1957      14.341 -18.518 -59.750  1.00 26.46           O
HETATM15165  O5  MAN B1957      12.440 -15.560 -60.696  1.00 26.50           O
HETATM15166  O6  MAN B1957      10.269 -17.445 -60.950  1.00 27.70           O
""",
  "linking_test_ccp4_other.pdb" : """
HETATM 7642  O3G DCT I   1      28.285 -27.771 -23.129  1.00 52.29           O
HETATM 7643  PG  DCT I   1      27.597 -26.624 -23.822  1.00 59.65           P
HETATM 7644  O1G DCT I   1      28.519 -25.410 -23.888  1.00 63.83           O
HETATM 7645  O2G DCT I   1      26.360 -26.264 -23.025  1.00 63.59           O
HETATM 7646  O3B DCT I   1      27.136 -27.016 -25.315  1.00 64.88           O
HETATM 7647  PB  DCT I   1      25.604 -26.763 -25.757  1.00 70.90           P
HETATM 7648  O1B DCT I   1      25.443 -27.147 -27.205  1.00 57.04           O
HETATM 7649  O2B DCT I   1      24.673 -27.588 -24.895  1.00 61.96           O
HETATM 7650  O3A DCT I   1      25.433 -25.167 -25.518  1.00 78.96           O
HETATM 7651  PA  DCT I   1      24.084 -24.455 -24.963  1.00 91.36           P
HETATM 7652  O1A DCT I   1      24.435 -23.045 -24.498  1.00 68.41           O
HETATM 7653  O2A DCT I   1      23.403 -25.253 -23.850  1.00 60.82           O
HETATM 7654  O5' DCT I   1      23.173 -24.412 -26.309  1.00 69.24           O
HETATM 7655  C5' DCT I   1      22.624 -25.634 -26.758  1.00 56.94           C
HETATM 7656  C4' DCT I   1      21.239 -25.488 -27.366  1.00 54.90           C
HETATM 7657  C3' DCT I   1      21.192 -26.185 -28.716  1.00 49.77           C
HETATM 7658  C2' DCT I   1      21.038 -25.070 -29.730  1.00 40.06           C
HETATM 7659  O4' DCT I   1      20.924 -24.128 -27.578  1.00 65.15           O
HETATM 7660  C1' DCT I   1      20.528 -23.892 -28.917  1.00 52.31           C
HETATM 7661  N1  DCT I   1      21.157 -22.640 -29.382  1.00 53.66           N
HETATM 7662  C6  DCT I   1      22.021 -21.959 -28.559  1.00 49.31           C
HETATM 7663  C5  DCT I   1      22.612 -20.776 -28.995  1.00 41.72           C
HETATM 7664  C4  DCT I   1      22.312 -20.298 -30.266  1.00 43.10           C
HETATM 7665  N4  DCT I   1      22.872 -19.170 -30.699  1.00 43.04           N
HETATM 7666  N3  DCT I   1      21.431 -20.987 -31.082  1.00 43.68           N
HETATM 7667  C2  DCT I   1      20.850 -22.158 -30.641  1.00 44.64           C
HETATM 7668  O2  DCT I   1      20.070 -22.778 -31.367  1.00 44.43           O
ATOM    167  P    Cd C  13      19.617 -19.479 -22.016  1.00 46.14           P
ATOM    168  O1P  Cd C  13      19.618 -19.829 -20.585  1.00 37.17           O
ATOM    169  O2P  Cd C  13      20.813 -18.844 -22.664  1.00 37.84           O
ATOM    170  O5*  Cd C  13      19.190 -20.751 -22.929  1.00 43.68           O
ATOM    171  C1*  Cd C  13      18.014 -20.805 -26.805  1.00 40.99           C
ATOM    172  C2*  Cd C  13      18.861 -22.059 -26.630  1.00 36.03           C
ATOM    173  C3*  Cd C  13      18.539 -22.506 -25.197  1.00 45.57           C
ATOM    174  C4*  Cd C  13      17.601 -21.426 -24.638  1.00 41.26           C
ATOM    175  C5*  Cd C  13      17.815 -21.006 -23.191  1.00 34.63           C
ATOM    176  O4*  Cd C  13      17.790 -20.288 -25.514  1.00 41.01           O
ATOM    177  N1   Cd C  13      18.634 -19.740 -27.625  1.00 40.70           N
ATOM    178  C2   Cd C  13      17.989 -19.394 -28.817  1.00 38.52           C
ATOM    179  N3   Cd C  13      18.516 -18.420 -29.604  1.00 29.24           N
ATOM    180  C4   Cd C  13      19.648 -17.820 -29.229  1.00 38.13           C
ATOM    181  C5   Cd C  13      20.322 -18.161 -28.007  1.00 42.18           C
ATOM    182  C6   Cd C  13      19.791 -19.126 -27.237  1.00 40.17           C
ATOM    183  O2   Cd C  13      16.945 -20.000 -29.087  1.00 37.42           O
ATOM    184  N4   Cd C  13      20.140 -16.868 -30.033  1.00 34.83           N
ATOM    246  P    Gd D   5      22.253 -15.271 -40.882  1.00 71.71           P
ATOM    247  O1P  Gd D   5      21.628 -14.785 -42.128  1.00 43.88           O
ATOM    248  O2P  Gd D   5      22.140 -14.525 -39.624  1.00 61.19           O
ATOM    249  O5*  Gd D   5      21.647 -16.686 -40.522  1.00 48.09           O
ATOM    250  O3*  Gd D   5      18.486 -19.251 -41.040  1.00 26.87           O
ATOM    251  C1*  Gd D   5      20.240 -18.959 -38.638  1.00 30.71           C
ATOM    252  C2*  Gd D   5      19.325 -17.910 -39.252  1.00 33.86           C
ATOM    253  C3*  Gd D   5      19.402 -18.250 -40.735  1.00 30.14           C
ATOM    254  C4*  Gd D   5      20.806 -18.812 -40.918  1.00 32.12           C
ATOM    255  C5*  Gd D   5      21.805 -17.790 -41.392  1.00 43.37           C
ATOM    256  O4*  Gd D   5      21.221 -19.274 -39.612  1.00 31.50           O
ATOM    257  N9   Gd D   5      20.898 -18.510 -37.405  1.00 39.95           N
ATOM    258  C8   Gd D   5      21.748 -17.431 -37.226  1.00 40.77           C
ATOM    259  N7   Gd D   5      22.170 -17.293 -35.994  1.00 41.36           N
ATOM    260  C5   Gd D   5      21.566 -18.346 -35.305  1.00 35.61           C
ATOM    261  C4   Gd D   5      20.784 -19.102 -36.161  1.00 37.92           C
ATOM    262  N1   Gd D   5      20.867 -19.840 -33.654  1.00 42.49           N
ATOM    263  C2   Gd D   5      20.119 -20.525 -34.581  1.00 37.44           C
ATOM    264  N3   Gd D   5      20.036 -20.196 -35.871  1.00 30.36           N
ATOM    265  C6   Gd D   5      21.637 -18.715 -33.939  1.00 33.17           C
ATOM    266  O6   Gd D   5      22.276 -18.200 -33.026  1.00 33.11           O
ATOM    267  N2   Gd D   5      19.447 -21.579 -34.104  1.00 37.54           N
""",
  "linking_test_two_ASN-NAG.pdb" : """
ATOM  10643  N   ASN C 281     -19.478  10.820  35.398  1.00 81.69      UNK  N
ATOM  10644  CA  ASN C 281     -19.872  11.829  34.423  1.00 81.69      UNK  C
ATOM  10645  C   ASN C 281     -18.678  12.767  34.352  1.00 81.69      UNK  C
ATOM  10646  O   ASN C 281     -17.651  12.498  34.967  1.00 81.69      UNK  O
ATOM  10647  CB  ASN C 281     -21.101  12.584  34.900  1.00 20.00      UNK  C
ATOM  10648  CG  ASN C 281     -20.973  12.990  36.338  1.00 20.00      UNK  C
ATOM  10649  OD1 ASN C 281     -19.993  13.641  36.709  1.00 20.00      UNK  O
ATOM  10650  ND2 ASN C 281     -21.942  12.609  37.167  1.00 20.00      UNK  N
ATOM  10662  N   ASN C 284     -17.638  14.353  38.268  1.00  1.00           N
ATOM  10663  CA  ASN C 284     -17.098  13.672  39.449  1.00  1.00           C
ATOM  10664  C   ASN C 284     -15.831  12.857  39.226  1.00  1.00           C
ATOM  10665  O   ASN C 284     -15.399  12.126  40.112  1.00  1.00           O
ATOM  10666  CB  ASN C 284     -18.176  12.838  40.159  1.00 20.00           C
ATOM  10667  CG  ASN C 284     -19.570  13.376  39.919  1.00 20.00           C
ATOM  10668  OD1 ASN C 284     -19.960  14.412  40.455  1.00 20.00           O
ATOM  10669  ND2 ASN C 284     -20.321  12.672  39.100  1.00 20.00           N
HETATM  135  C1  NAG C 751     -21.757  12.993  38.644  1.00 88.86           C
HETATM  136  C2  NAG C 751     -22.945  12.574  39.656  1.00 88.86           C
HETATM  137  C3  NAG C 751     -23.602  13.737  40.389  1.00 88.86           C
HETATM  138  C4  NAG C 751     -24.038  14.792  39.375  1.00 88.86           C
HETATM  139  C5  NAG C 751     -23.012  15.048  38.266  1.00 88.86           C
HETATM  140  C6  NAG C 751     -22.855  16.548  38.038  1.00 88.86           C
HETATM  141  C7  NAG C 751     -23.929  10.591  38.674  1.00 88.86           C
HETATM  142  C8  NAG C 751     -22.762  10.103  37.868  1.00 88.86           C
HETATM  143  N2  NAG C 751     -23.970  11.898  38.896  1.00 88.86           N
HETATM  144  O3  NAG C 751     -22.675  14.320  41.308  1.00 88.86           O
HETATM  145  O4  NAG C 751     -25.267  14.372  38.772  1.00 88.86           O
HETATM  146  O5  NAG C 751     -21.747  14.436  38.578  1.00 88.86           O
HETATM  147  O6  NAG C 751     -22.602  16.781  36.649  1.00 88.86           O
HETATM  148  O7  NAG C 751     -24.791   9.842  39.102  1.00 88.86           O
""",
  "linking_test_SER_C.pdb" : """
ATOM      1  P   C     643     -75.543  -1.303 308.744  1.00 91.01      D16S P
ATOM      2  C5' C     643     -77.315  -3.206 308.385  1.00 99.83      D16S C
ATOM      3  O5' C     643     -75.992  -2.834 308.736  1.00 92.18      D16S O
ATOM      4  C4' C     643     -77.466  -4.699 308.231  1.00 91.83      D16S C
ATOM      5  O4' C     643     -76.541  -5.198 307.222  1.00 86.11      D16S O
ATOM      6  C3' C     643     -77.142  -5.522 309.466  1.00 83.89      D16S C
ATOM      7  O3' C     643     -78.192  -5.550 310.411  1.00 90.22      D16S O
ATOM      8  C2' C     643     -76.796  -6.877 308.870  1.00 75.13      D16S C
ATOM      9  O2' C     643     -77.966  -7.558 308.435  1.00 66.17      D16S O
ATOM     10  C1' C     643     -76.015  -6.448 307.629  1.00 71.55      D16S C
ATOM     11  N1  C     643     -74.568  -6.277 307.938  1.00 72.53      D16S N
ATOM     12  C2  C     643     -73.785  -7.404 308.299  1.00 72.36      D16S C
ATOM     13  O2  C     643     -74.298  -8.549 308.341  1.00 82.06      D16S O
ATOM     14  N3  C     643     -72.467  -7.201 308.589  1.00 67.68      D16S N
ATOM     15  C4  C     643     -71.937  -5.966 308.532  1.00 71.55      D16S C
ATOM     16  N4  C     643     -70.646  -5.814 308.832  1.00 68.18      D16S N
ATOM     17  C5  C     643     -72.701  -4.825 308.168  1.00 79.88      D16S C
ATOM     18  C6  C     643     -73.988  -5.025 307.883  1.00 75.82      D16S C
ATOM     19  OP1 C     643     -76.597  -0.455 309.339  1.00 76.88      D16S O
ATOM     20  OP2 C     643     -74.175  -1.305 309.279  1.00 82.47      D16S O1-
ATOM      1  N   SER   115     -76.499  -0.597 309.288  1.00 76.74      DS08 N
ATOM      2  CA  SER   115     -75.220   0.093 309.453  1.00 83.12      DS08 C
ATOM      3  C   SER   115     -75.112   1.417 308.701  1.00 78.72      DS08 C
ATOM      4  O   SER   115     -74.545   2.372 309.209  1.00 72.36      DS08 O
ATOM      5  CB  SER   115     -74.096  -0.829 309.007  1.00 81.98      DS08 C
ATOM      6  OG  SER   115     -74.613  -2.114 308.687  1.00 93.10      DS08 O
""",
  "linking_test_ASN-NAG-not-THR.pdb" : """
ATOM  11100  N   ASN C 281     -20.823  10.632  36.650  1.00128.11           N
ATOM  11101  CA  ASN C 281     -21.542  11.723  36.020  1.00126.84           C
ATOM  11102  C   ASN C 281     -20.321  12.535  35.687  1.00133.47           C
ATOM  11103  O   ASN C 281     -19.289  12.313  36.291  1.00141.10           O
ATOM  11104  CB  ASN C 281     -22.554  12.488  36.879  1.00132.74           C
ATOM  11105  CG  ASN C 281     -22.348  12.326  38.361  1.00136.41           C
ATOM  11106  OD1 ASN C 281     -21.585  11.488  38.806  1.00137.40           O
ATOM  11107  ND2 ASN C 281     -23.058  13.149  39.138  1.00142.79           N
ATOM  11108  N   GLY C 282     -20.400  13.479  34.769  1.00132.24           N
ATOM  11109  CA  GLY C 282     -19.205  14.196  34.360  1.00136.73           C
ATOM  11110  C   GLY C 282     -18.431  14.831  35.495  1.00138.96           C
ATOM  11111  O   GLY C 282     -17.216  14.977  35.444  1.00139.29           O
ATOM  11112  N   THR C 283     -19.149  15.231  36.515  1.00133.75           N
ATOM  11113  CA  THR C 283     -18.571  15.870  37.669  1.00134.08           C
ATOM  11114  C   THR C 283     -17.581  15.120  38.564  1.00137.98           C
ATOM  11115  O   THR C 283     -16.604  15.691  39.022  1.00143.56           O
ATOM  11116  CB  THR C 283     -19.722  16.166  38.613  1.00144.68           C
ATOM  11117  OG1 THR C 283     -20.222  14.917  39.104  1.00143.45           O
ATOM  11118  CG2 THR C 283     -20.836  16.883  37.874  1.00141.28           C
HETATM14227  C1  NAG C1281     -23.045  13.197  40.578  1.00155.05           C
HETATM14228  C2  NAG C1281     -24.038  13.231  41.733  1.00155.85           C
HETATM14229  C3  NAG C1281     -24.245  14.639  42.281  1.00159.74           C
HETATM14230  C4  NAG C1281     -22.956  15.416  42.476  1.00164.53           C
HETATM14231  C5  NAG C1281     -21.978  15.164  41.334  1.00162.17           C
HETATM14232  C6  NAG C1281     -20.583  15.654  41.678  1.00163.10           C
HETATM14233  C7  NAG C1281     -26.167  13.389  40.586  1.00151.28           C
HETATM14234  C8  NAG C1281     -26.553  12.787  39.271  1.00144.48           C
HETATM14235  N2  NAG C1281     -25.308  12.684  41.309  1.00153.73           N
HETATM14236  O3  NAG C1281     -24.905  14.555  43.545  1.00169.63           O
HETATM14237  O4  NAG C1281     -23.259  16.819  42.509  1.00164.51           O
HETATM14238  O5  NAG C1281     -21.850  13.775  41.067  1.00160.58           O
HETATM14239  O6  NAG C1281     -19.653  14.829  40.962  1.00160.49           O
HETATM14240  O7  NAG C1281     -26.606  14.456  40.975  1.00157.16           O
""",
  "linking_test_HEM_TYR.pdb" : """
ATOM      6  CG  TYR A 140      98.997   4.611  10.805  1.00 38.35           C
ATOM      7  CD1 TYR A 140      99.820   3.648  11.373  1.00 36.78           C
ATOM      8  CD2 TYR A 140      98.847   5.821  11.476  1.00 37.30           C
ATOM      9  CE1 TYR A 140     100.465   3.853  12.607  1.00 37.29           C
ATOM     10  CE2 TYR A 140      99.502   6.047  12.739  1.00 42.05           C
ATOM     11  CZ  TYR A 140     100.316   5.053  13.295  1.00 39.27           C
ATOM     12  OH  TYR A 140     101.024   5.185  14.508  1.00 47.51           O
ATOM     13  N   CYS A 470     106.717   6.301  20.283  1.00 34.42           N
ATOM     14  CA  CYS A 470     107.171   7.557  19.724  1.00 35.14           C
ATOM     15  C   CYS A 470     108.672   7.540  19.477  1.00 37.62           C
ATOM     16  O   CYS A 470     109.451   7.026  20.285  1.00 39.29           O
ATOM     17  CB  CYS A 470     106.833   8.698  20.693  1.00 30.41           C
ATOM     18  SG  CYS A 470     107.139  10.328  20.022  1.00 34.85           S
ATOM     19  CAA HEM A 601     102.373   6.235  20.462  1.00 29.04           C
ATOM     20  CAB HEM A 601     105.667  15.514  22.495  1.00 32.15           C
ATOM     21  CAC HEM A 601     108.356  13.649  15.330  1.00 33.91           C
ATOM     22  CAD HEM A 601     105.102   6.382  15.766  1.00 29.31           C
ATOM     23  NA  HEM A 601     103.895   9.666  19.977  1.00 30.01           N
ATOM     24  CBA HEM A 601     101.025   6.147  19.746  1.00 31.00           C
ATOM     25  CBB HEM A 601     105.408  16.127  23.646  1.00 39.39           C
ATOM     26  CBC HEM A 601     109.441  14.426  15.313  1.00 37.97           C
ATOM     27  CBD HEM A 601     103.624   6.099  15.536  1.00 31.13           C
ATOM     28  NB  HEM A 601     104.842  12.376  20.622  1.00 28.66           N
ATOM     29  CGA HEM A 601     100.437   4.775  19.956  1.00 31.11           C
ATOM     30  CGD HEM A 601     103.294   4.677  15.216  1.00 36.63           C
ATOM     31  ND  HEM A 601     105.110   9.805  17.433  1.00 28.27           N
ATOM     32  CHA HEM A 601     104.002   7.823  18.350  1.00 28.03           C
ATOM     33  CHB HEM A 601     103.388  10.956  22.005  1.00 29.20           C
ATOM     34  CHC HEM A 601     106.256  14.228  19.862  1.00 31.30           C
ATOM     35  CHD HEM A 601     106.404  11.313  15.985  1.00 28.55           C
ATOM     36  CMA HEM A 601     101.920   8.258  22.900  1.00 30.34           C
ATOM     37  CMB HEM A 601     103.702  13.324  24.099  1.00 31.20           C
ATOM     38  CMC HEM A 601     108.238  15.628  17.974  1.00 30.92           C
ATOM     39  CMD HEM A 601     106.699   8.724  14.156  1.00 27.61           C
ATOM     40  C1A HEM A 601     103.617   8.383  19.542  1.00 31.73           C
ATOM     41  C1B HEM A 601     104.130  12.102  21.772  1.00 30.51           C
ATOM     42  C1C HEM A 601     106.579  13.710  18.629  1.00 31.98           C
ATOM     43  C1D HEM A 601     105.856  10.084  16.301  1.00 28.65           C
ATOM     44  O1A HEM A 601      99.301   4.511  19.494  1.00 36.78           O
ATOM     45  O1D HEM A 601     104.283   3.883  15.263  1.00 33.51           O
ATOM     46  C2A HEM A 601     102.863   7.699  20.570  1.00 33.39           C
ATOM     47  C2B HEM A 601     104.328  13.232  22.690  1.00 33.80           C
ATOM     48  C2C HEM A 601     107.499  14.297  17.699  1.00 30.13           C
ATOM     49  C2D HEM A 601     105.949   8.853  15.516  1.00 36.74           C
ATOM     50  O2A HEM A 601     101.111   3.928  20.599  1.00 35.99           O
ATOM     51  O2D HEM A 601     102.032   4.472  14.888  1.00 47.04           O
ATOM     52  C3A HEM A 601     102.693   8.556  21.587  1.00 33.17           C
ATOM     53  C3B HEM A 601     105.108  14.136  22.069  1.00 30.61           C
ATOM     54  C3C HEM A 601     107.530  13.499  16.617  1.00 32.10           C
ATOM     55  C3D HEM A 601     105.302   7.882  16.179  1.00 35.98           C
ATOM     56  C4A HEM A 601     103.337   9.814  21.234  1.00 28.53           C
ATOM     57  C4B HEM A 601     105.443  13.603  20.772  1.00 26.80           C
ATOM     58  C4C HEM A 601     106.652  12.369  16.852  1.00 31.29           C
ATOM     59  C4D HEM A 601     104.750   8.469  17.397  1.00 29.94           C
ATOM     60  NC  HEM A 601     106.086  12.539  18.099  1.00 32.94           N
ATOM     61 FE   HEM A 601     105.139  11.046  19.090  1.00 32.47          Fe
""",
  "linking_test_HEM_TYR.params" : """
pdb_interpretation {
  automatic_linking {
    automatic_linking {
      link_metals = True
      link_residues = True
      inter_residue_bond_cutoff = 2.5
    }
  }
}
geometry_restraints.edits {
  bond {
    atom_selection_1 = chain 'A' and resid ' 601 ' and name 'FE  ' and \
                       altloc ' '
    atom_selection_2 = chain 'A' and resid ' 470 ' and name ' SG ' and \
                       altloc ' '
    distance_ideal = 2.33
    sigma = 0.02
  }
}
""",
  "linking_test_MAN-before-ASN.pdb" : """
HETATM   32  C1  MAN A  13     -23.598  63.781-180.576  1.00 20.00           C
HETATM   33  O5  MAN A  13     -24.198  64.994-181.028  1.00 20.00           O
HETATM   34  C5  MAN A  13     -25.269  65.301-180.140  1.00 20.00           C
HETATM   35  C6  MAN A  13     -26.064  66.487-180.693  1.00 20.00           C
HETATM   36  O6  MAN A  13     -26.588  66.150-181.978  1.00 20.00           O
HETATM   37  C4  MAN A  13     -24.712  65.661-178.762  1.00 20.00           C
HETATM   38  O4  MAN A  13     -25.786  65.996-177.882  1.00 20.00           O
HETATM   39  C3  MAN A  13     -23.949  64.456-178.202  1.00 20.00           C
HETATM   40  O3  MAN A  13     -23.296  64.820-176.985  1.00 20.00           O
HETATM   41  C2  MAN A  13     -22.904  64.015-179.233  1.00 20.00           C
HETATM   42  O2  MAN A  13     -21.915  65.035-179.378  1.00 20.00           O
ATOM     43  O   ASN A 210     -64.131  35.488 -16.482  1.00 77.95           O
ATOM     44  N   ASN A 210     -63.707  32.434 -14.675  1.00 35.76           N
ATOM     45  CA  ASN A 210     -64.179  33.784 -14.816  1.00 64.69           C
ATOM     46  C   ASN A 210     -63.605  34.499 -16.022  1.00 72.64           C
ATOM     47  CB  ASN A 210     -63.830  34.555 -13.567  1.00 77.34           C
""",
  "linking_test_partial_alt_loc_glyco.pdb" : """
HETATM 3615  C1 ANDG B 401      11.587 187.498  -2.208  0.50 17.33           C
HETATM 3616  C2 ANDG B 401      10.855 188.231  -1.091  0.50 17.55           C
HETATM 3617  C3 ANDG B 401      10.531 189.647  -1.547  0.50 16.69           C
HETATM 3618  C4 ANDG B 401      11.775 190.352  -2.092  0.50 16.31           C
HETATM 3619  C5 ANDG B 401      12.356 189.491  -3.197  0.50 16.42           C
HETATM 3620  C6 ANDG B 401      13.546 190.144  -3.891  0.50 17.65           C
HETATM 3621  C7 ANDG B 401       9.454 186.869   0.419  0.50 17.55           C
HETATM 3622  C8 ANDG B 401       8.079 186.314   0.643  0.50 16.99           C
HETATM 3623  O  ANDG B 401      12.733 188.238  -2.631  0.50 16.58           O
HETATM 3624  O3 ANDG B 401      10.019 190.380  -0.453  0.50 17.04           O
HETATM 3625  O4 ANDG B 401      11.386 191.615  -2.628  0.50 15.06           O
HETATM 3626  O6 ANDG B 401      14.616 190.314  -2.961  0.50 17.54           O
HETATM 3627  O7 ANDG B 401      10.336 186.702   1.248  0.50 18.38           O
HETATM 3628  N2 ANDG B 401       9.620 187.548  -0.722  0.50 17.55           N
HETATM 3629  O1LANDG B 401      10.701 187.349  -3.322  0.50 17.86           O
HETATM 3630  C1 BNAG B 402      11.696 187.522  -2.281  0.50 19.30           C
HETATM 3631  C2 BNAG B 402      10.842 188.218  -1.233  0.50 19.73           C
HETATM 3632  C3 BNAG B 402      10.507 189.651  -1.615  0.50 18.14           C
HETATM 3633  C4 BNAG B 402      11.754 190.377  -2.114  0.50 17.57           C
HETATM 3634  C5 BNAG B 402      12.374 189.564  -3.233  0.50 17.69           C
HETATM 3635  C6 BNAG B 402      13.544 190.278  -3.896  0.50 19.10           C
HETATM 3636  C7 BNAG B 402       9.520 186.549  -0.070  0.50 22.13           C
HETATM 3637  C8 BNAG B 402       8.215 185.818   0.005  0.50 22.97           C
HETATM 3638  N2 BNAG B 402       9.630 187.444  -1.049  0.50 20.88           N
HETATM 3639  O1 BNAG B 402      12.174 186.314  -1.691  0.50 19.78           O
HETATM 3640  O3 BNAG B 402      10.007 190.337  -0.478  0.50 17.72           O
HETATM 3641  O4 BNAG B 402      11.350 191.649  -2.605  0.50 15.73           O
HETATM 3642  O5 BNAG B 402      12.803 188.319  -2.692  0.50 18.12           O
HETATM 3643  O6 BNAG B 402      14.675 190.273  -3.024  0.50 18.91           O
HETATM 3644  O7 BNAG B 402      10.422 186.331   0.724  0.50 24.09           O
HETATM 3645  C1 AFUC B 403       8.585 190.434  -0.438  0.50 18.44           C
HETATM 3646  C1 BFUC B 403       8.574 190.425  -0.442  0.50 18.71           C
HETATM 3647  C2  FUC B 403       8.154 190.791   0.979  1.00 19.23           C
HETATM 3648  C3  FUC B 403       8.688 192.169   1.357  1.00 19.53           C
HETATM 3649  C4  FUC B 403       8.143 193.181   0.356  1.00 20.34           C
HETATM 3650  C5  FUC B 403       8.484 192.740  -1.077  1.00 20.67           C
HETATM 3651  C6  FUC B 403       7.889 193.667  -2.123  1.00 22.65           C
HETATM 3652  O2  FUC B 403       8.618 189.805   1.908  1.00 20.29           O
HETATM 3653  O3  FUC B 403       8.325 192.496   2.712  1.00 19.86           O
HETATM 3654  O4  FUC B 403       6.720 193.304   0.464  1.00 20.42           O
HETATM 3655  O5  FUC B 403       8.050 191.390  -1.356  1.00 18.33           O
HETATM 3656  C1 AGAL B 404      12.208 192.729  -2.260  0.50 14.68           C
HETATM 3657  C1 BGAL B 404      12.218 192.731  -2.262  0.50 15.08           C
HETATM 3658  C2 AGAL B 404      11.977 193.809  -3.313  0.50 14.65           C
HETATM 3659  C2 BGAL B 404      11.974 193.813  -3.310  0.50 14.87           C
HETATM 3660  C3  GAL B 404      12.714 195.089  -2.980  1.00 14.36           C
HETATM 3661  C4 AGAL B 404      12.326 195.524  -1.566  0.50 14.68           C
HETATM 3662  C4 BGAL B 404      12.331 195.517  -1.566  0.50 14.89           C
HETATM 3663  C5 AGAL B 404      12.545 194.397  -0.567  0.50 14.42           C
HETATM 3664  C5 BGAL B 404      12.656 194.395  -0.596  0.50 14.74           C
HETATM 3665  C6 AGAL B 404      12.057 194.809   0.817  0.50 14.38           C
HETATM 3666  C6 BGAL B 404      12.314 194.806   0.831  0.50 14.69           C
HETATM 3667  O2 AGAL B 404      12.381 193.340  -4.606  0.50 14.14           O
HETATM 3668  O2 BGAL B 404      12.367 193.348  -4.608  0.50 14.30           O
HETATM 3669  O3  GAL B 404      12.351 196.082  -3.971  1.00 15.47           O
HETATM 3670  O4 AGAL B 404      10.951 195.945  -1.547  0.50 14.29           O
HETATM 3671  O4 BGAL B 404      10.932 195.846  -1.506  0.50 14.61           O
HETATM 3672  O5 AGAL B 404      11.840 193.232  -0.984  0.50 15.24           O
HETATM 3673  O5 BGAL B 404      11.921 193.232  -0.964  0.50 15.68           O
HETATM 3674  O6 AGAL B 404      12.320 193.790   1.794  0.50 13.65           O
HETATM 3675  O6 BGAL B 404      12.433 193.702   1.739  0.50 14.16           O
""",
  'linking_test_CM-SO4.pdb' : '''
HETATM 5816  N26  CM C   4      11.872  46.521  11.694  1.00 15.75      A    N
HETATM 6475  S   SO4 D   4      12.593  45.477   7.849  1.00 30.00           S
HETATM 6476  O1  SO4 D   4      13.522  46.224   6.981  1.00 30.00           O
HETATM 6477  O2  SO4 D   4      11.558  44.829   7.018  1.00 30.00           O
HETATM 6478  O3  SO4 D   4      13.343  44.439   8.588  1.00 30.00           O
HETATM 6479  O4  SO4 D   4      11.949  46.390   8.816  1.00 30.00           O
''',
  "CM.cif" : """
data_comp_CM
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
_chem_comp_atom.charge
CM N26    N N  0.000 0
""",
  'linking_test_exclusion_SO4.pdb' : '''
CRYST1  126.007  546.508  585.856  90.00  90.00  90.00 P 2 21 21
ATOM   5781  N   ARG A 782     -68.913-438.467 717.772  1.00 96.17           N
ATOM   5782  CA  ARG A 782     -68.051-437.794 718.737  1.00 96.10           C
ATOM   5783  C   ARG A 782     -68.057-438.454 720.114  1.00 96.19           C
ATOM   5784  O   ARG A 782     -67.489-437.889 721.054  1.00 96.06           O
ATOM   5785  CB  ARG A 782     -68.457-436.323 718.866  1.00 95.82           C
ATOM   5786  CG  ARG A 782     -67.625-435.372 718.021  1.00 95.68           C
ATOM   5787  CD  ARG A 782     -66.340-434.995 718.736  1.00 95.65           C
ATOM   5788  NE  ARG A 782     -65.486-434.139 717.918  1.00 95.53           N
ATOM   5789  CZ  ARG A 782     -64.442-434.578 717.222  1.00 95.65           C
ATOM   5790  NH1 ARG A 782     -64.123-435.865 717.246  1.00 95.88           N1+
ATOM   5791  NH2 ARG A 782     -63.715-433.731 716.504  1.00 95.52           N
ATOM   4791  N   LYS A 660     -54.507-433.397 715.702  1.00 89.00           N
ATOM   4792  CA  LYS A 660     -55.052-434.689 715.302  1.00 89.16           C
ATOM   4793  C   LYS A 660     -53.950-435.729 715.157  1.00 89.32           C
ATOM   4794  O   LYS A 660     -54.180-436.916 715.414  1.00 89.49           O
ATOM   4795  CB  LYS A 660     -55.827-434.552 713.989  1.00 89.09           C
ATOM   4796  CG  LYS A 660     -57.191-433.881 714.110  1.00 88.95           C
ATOM   4797  CD  LYS A 660     -58.164-434.700 714.947  1.00 89.07           C
ATOM   4798  CE  LYS A 660     -59.455-433.926 715.203  1.00 88.92           C
ATOM   4799  NZ  LYS A 660     -60.489-434.725 715.925  1.00 89.03           N1+
HETATM 6886  O1  SO4 A 923     -60.519-433.073 715.370  0.00 86.34           O
HETATM 6887  O2  SO4 A 923     -62.688-432.381 714.671  1.00 86.19           O
HETATM 6888  O3  SO4 A 923     -62.430-433.668 716.666  0.10 86.40           O
HETATM 6889  O4  SO4 A 923     -61.695-431.404 716.606  1.00 86.05           O
HETATM 6890  S   SO4 A 923     -61.834-432.632 715.827  1.00 86.25           S
  ''',
  'so4_exclude.phil' : '''
pdb_interpretation {
  exclude_from_automatic_linking {
    selection_1 = "resname SO4"
    selection_2 = all
  }
}
  ''',
  'linking_test_cyclic_d_amino_acid.pdb' : '''
CRYST1   62.592   88.672   79.544  90.00  90.00  90.00 P 1           0
SCALE1      0.015976  0.000000  0.000000        0.00000
SCALE2      0.000000  0.011278  0.000000        0.00000
SCALE3      0.000000  0.000000  0.012572        0.00000
ATOM      1  N   ALA A   1     139.541 141.181  96.260  1.00 20.00           N
ATOM      2  CA  ALA A   1     138.637 140.412  97.146  1.00 20.00           C
ATOM      3  C   ALA A   1     139.314 139.161  97.689  1.00 20.00           C
ATOM      4  O   ALA A   1     140.095 138.531  96.963  1.00 20.00           O
ATOM      5  CB  ALA A   1     137.356 140.101  96.417  1.00 20.00           C
ATOM      6  O   ALA A   2     137.325 136.971  99.601  1.00 30.00           O
ATOM      7  N   ALA A   2     139.043 138.844  98.952  1.00 30.00           N
ATOM      8  CA  ALA A   2     139.595 137.635  99.598  1.00 30.00           C
ATOM      9  C   ALA A   2     138.408 136.872 100.186  1.00 30.00           C
ATOM     10  CB  ALA A   2     140.622 138.035 100.627  1.00 30.00           C
ATOM     11  O   GLY A   3     137.830 133.641 100.447  1.00 30.00           O
ATOM     12  N   GLY A   3     138.559 136.152 101.290  1.00 30.00           N
ATOM     13  CA  GLY A   3     137.375 135.482 101.846  1.00 30.00           C
ATOM     14  C   GLY A   3     136.942 134.331 100.969  1.00 30.00           C
ATOM     15  O   ALA A   4     136.226 132.849  97.746  1.00 30.00           O
ATOM     16  N   ALA A   4     135.635 134.109 100.811  1.00 30.00           N
ATOM     17  CA  ALA A   4     135.185 133.051  99.882  1.00 30.00           C
ATOM     18  C   ALA A   4     135.449 133.514  98.449  1.00 30.00           C
ATOM     19  CB  ALA A   4     133.736 132.690 100.110  1.00 30.00           C
ATOM     20  O   ALA A   5     136.629 136.668  96.001  1.00 30.00           O
ATOM     21  N   ALA A   5     134.856 134.630  98.036  1.00 30.00           N
ATOM     22  CA  ALA A   5     134.990 135.102  96.644  1.00 30.00           C
ATOM     23  C   ALA A   5     136.429 135.521  96.406  1.00 30.00           C
ATOM     24  CB  ALA A   5     134.055 136.265  96.409  1.00 30.00           C
ATOM     25  O   ALA A   6     139.834 134.403  94.458  1.00 30.00           O
ATOM     26  N   ALA A   6     137.389 134.638  96.643  1.00 30.00           N
ATOM     27  CA  ALA A   6     138.809 135.028  96.518  1.00 30.00           C
ATOM     28  C   ALA A   6     139.123 135.234  95.043  1.00 30.00           C
ATOM     29  CB  ALA A   6     139.676 133.969  97.130  1.00 30.00           C
ATOM     30  O   ALA A   7     137.722 138.645  93.058  1.00 30.00           O
ATOM     31  N   ALA A   7     138.555 136.270  94.446  1.00 30.00           N
ATOM     32  CA  ALA A   7     138.715 136.491  92.997  1.00 30.00           C
ATOM     33  C   ALA A   7     138.736 137.996  92.776  1.00 30.00           C
ATOM     34  CB  ALA A   7     137.571 135.846  92.258  1.00 30.00           C
ATOM     35  N   DPR A   8     139.914 138.630  92.485  1.00 20.00           N
ATOM     36  CA  DPR A   8     139.816 140.123  92.327  1.00 20.00           C
ATOM     37  CB  DPR A   8     140.573 140.423  91.018  1.00 20.00           C
ATOM     38  CG  DPR A   8     141.601 139.328  90.913  1.00 20.00           C
ATOM     39  CD  DPR A   8     141.393 138.483  92.150  1.00 20.00           C
ATOM     40  C   DPR A   8     140.000 141.439  93.121  1.00 20.00           C
ATOM     41  O   DPR A   8     139.549 142.442  92.591  1.00 20.00           O
ATOM     42  N   DAL A   9     140.609 141.365  93.957  1.00 20.00           N
ATOM     43  CA  DAL A   9     140.665 142.742  94.522  1.00 20.00           C
ATOM     44  CB  DAL A   9     142.190 142.550  94.599  1.00 20.00           C
ATOM     45  C   DAL A   9     140.049 142.247  95.882  1.00 20.00           C
ATOM     46  O   DAL A   9     139.571 143.185  96.612  1.00 20.00           O
END
  ''',
        }

links = {
  "linking_test_ASN-NAG.pdb" : [21, 22],
  "linking_test_ASN-NAG-altloc1.pdb" : [35, 37],
  "linking_test_ASN-NAG-altloc2.pdb" : [42, 44],
  "linking_test_ASN-NAG-altloc3.pdb" : [38, 40],
  "linking_test_ASN-NAG-altloc4.pdb" : [43, 43],
  "linking_test_NAG-FU4.pdb" : [24, 25],
  "linking_test_NAG-NAG.pdb" : [35, 37],
  'linking_test_BGC-BGC.pdb' : [0, 1],
  'linking_test_CYS-VSP.pdb' : [66, 68],
  "linking_test_LEU-CSY-VAL.pdb" : [13,15],
  'linking_test_DT-4MF-DA.pdb' : [69, 69],
  'linking_test_LYS-ABA-GLU.pdb' : [23,23],
  "linking_test_ALA-ALA-ALA.pdb" : [32,32],
  "linking_test_MAN-SER.pdb" : [17, 17],
  "linking_test_3g2j-LLP.pdb" : [44, 44],
  "linking_test_ASN_A-NAG_B.pdb" : [21, 22],
  "linking_test_nstd_rna_dna_h_bond.pdb" : [0,0],
  "linking_test_nstd_rna_dna.pdb" : [0,1],
  "linking_test_Mg_HOH.pdb" : [0,6],                 #6], # metal coordination
  "linking_test_Mg_HOH_CRYST1.pdb" : [0,6],
  "linking_test_Mg_EDT.pdb" : [19,25],               #25],
  "linking_test_1jbe_ALA-SNN-ACY-ALA.pdb" : [10,13],
  "linking_test_3gmq_NAG-FUC.pdb" : [24,25],
  "linking_test_CD_GHE_A_B.pdb" : [0,0],             #4],
  "linking_test_XYP_XYP.pdb" : [18,19],
  "linking_test_ALY_MCM.pdb" : [11,12], # links AA with quasi-AA
  # cross link not working - DL: I think it works now
  "linking_test_cyclic_side_chain.pdb" : [67,68], # side chain cross link
  "linking_test_cyclic_main_chain.pdb" : [34,35], # main chain cyclic
  "linking_test_over_valence.pdb" : [6,6],
  "linking_test_c2_c6.pdb" : [21,22],
  "linking_test_ccp4_other.pdb" : [71,71],
  "linking_test_SER_C.pdb" : [26,27],
  #
  "linking_test_two_ASN-NAG.pdb" : [28,29,29,29],
  #
  "linking_test_ASN-NAG-not-THR.pdb" : [32,33],
  #
  "linking_test_HEM_TYR.pdb" : [63,63],
  #
  "linking_test_CYS_CYS_alt_loc.pdb" : [17,17],
  #
  "linking_test_MAN-before-ASN.pdb" : [15,15],
  "linking_test_partial_alt_loc_glyco.pdb" : [63,67],
  'linking_test_CM-SO4.pdb' : [4,4],
  'linking_test_exclusion_SO4.pdb' :[24,22],
  'linking_test_cyclic_d_amino_acid.pdb' : [46,47],
  }

def run_and_test(cmd, pdb, i, skip_links=False):
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert (result.return_code == 0)
  for line in result.stdout_lines :
    if ("Write PDB file" in line):
      break
  else :
    raise RuntimeError("Missing expected log output")
  print("OK")
  # test .geo
  f=open(pdb.replace(".pdb", "_minimized.geo"), "r")
  lines = f.read()
  f.close()
  if pdb=='linking_test_cyclic_main_chain.pdb':
    if i==1:
      assert lines.find('link_TRANS | restraints: 1')>-1
  bonds = 0
  for line in lines.splitlines():
    for bond_like in ["Bond | covalent geometry | restraints:",
                      'Bond-like restraints:',
                      'Bond | Metal coordination | restraints',
                      'Bond | User supplied | restraints',
                      'Bond | Custom Glycosidic | restraints',
                      'Bond | link_BETA1-6 | restraints',
                      'Bond | link_NAG-ASN | restraints',
                      'Bond | link_BETA1-3 | restraints',
                      'Bond | link_BETA1-4 | restraints',
                      'Bond | link_ALPHA2-6 | restraints',
                      'Bond | Misc. | restraints',
                      'Bond | Disulphide bridge | restraints',
                      'Bond | link_TRANS | restraints',
                      ]:
      if line.find(bond_like)>-1:
        print('line',line)
        print('bond_like',bond_like)
        print('Adding %s for "%s"' % (int(line.split()[-1]), line))
        bonds += int(line.split()[-1])
    if line.find('Bond angle')>-1: break
  assert bonds == links[pdb][i], "found %d bonds but expected %s! File: %s" % (
    bonds,
    links[pdb][i],
    pdb)
  new_geo = pdb.replace(".pdb", "_minimized_%d.geo" % i)
  if (os.path.isfile(new_geo)):
    os.remove(new_geo)
  os.rename(pdb.replace(".pdb", "_minimized.geo"), new_geo)
  print("OK")
  if skip_links: return
  number_of_links=0
  fname = pdb.replace(".pdb", "_minimized.pdb")
  f=open(fname, "r")
  lines = f.readlines()
  f.close()
  for line in lines:
    if line.find("LINK")>-1:
      number_of_links+=1
  if i==0:
    expected = 0
  else:
    expected = links[pdb][i]-links[pdb][0]
    if pdb in ['linking_test_LEU-CSY-VAL.pdb']:
      expected -= 2 # peptide-like link
  if pdb in ["linking_test_HEM_TYR.pdb"]: # unset defined edits
    expected += 1
  #
  assert number_of_links == expected, "found %d LINK but expected %s! File: %s" % (
    number_of_links,
    expected,
    fname)
  if 0:
    cmd = "phenix.start_coot --no-guano --pdb %s" % pdb.replace(".pdb",
                                                                "_minimized.pdb"
                                                                )
    os.system(cmd)

def run(only_i=None):
  try: only_i=int(only_i)
  except ValueError: only_i=None
  except TypeError: only_i=None
  longer_tests = [
    "linking_test_two_ASN-NAG.pdb",
    ]
  # test exclusion
  if only_i is not None and only_i==1:
    cifs = ""
    for pdb in pdbs:
      f=open(pdb, "w")
      f.write(pdbs[pdb])
      f.close()
      if pdb.endswith(".phil"): cifs += " %s" % pdb
    for pdb in pdbs:
      if pdb=='linking_test_exclusion_SO4.pdb':
        cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
        cmd += ' link_small_molecules=True'
        print(cmd)
        run_and_test(cmd, pdb, 0, True)
        cmd += "  %s" % (cifs)
        print(cmd)
        run_and_test(cmd, pdb, 1, True)
  #
  cifs = ""
  for pdb in pdbs:
    f=open(pdb, "w")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
  j=0
  for k, pdb in enumerate(sorted(pdbs)):
    #break
    if pdb.endswith(".cif"): continue
    if pdb.endswith(".params"): continue
    if pdb.endswith(".phil"): continue
    if pdb in longer_tests: continue
    if pdb.find("CD_GHE")>-1: continue
    if pdb.find("partial")>-1: continue
    #if pdb.find('SO4')==-1: continue
    if pdb in ['linking_test_exclusion_SO4.pdb']: continue
    print('pdb '*10,k-10,pdb)
    j+=1
    if only_i is not None and only_i!=j: continue
    for i in range(2):
      log_filename = "%s_%d.log" % (pdb, i)
      cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
      cmd += " link_all=%d link_carbohydrate=%d %s" % (i, i, cifs)
      cmd += " link_ligand=%d" % i
      if pdb in ["linking_test_ccp4_other.pdb"]:
        if i: cmd += " secondary_structure.enabled=1"
      if pdb in ['linking_test_cyclic_d_amino_acid.pdb',
                 'linking_test_cyclic_main_chain.pdb',
                 'linking_test_cyclic_side_chain.pdb',
                 ]:
        if not i: cmd += ' link_residues=False'
      if pdb.replace(".pdb", ".params") in pdbs:
        cmd += " %s" % pdb.replace(".pdb", ".params")
      print('='*80)
      print("test number: %s (%s)" % (j,i))
      print('='*80)
      print(cmd)
      run_and_test(cmd, pdb, i)

  for pdb in sorted(pdbs):
    if pdb not in longer_tests: continue
    j+=1
    if only_i is not None and only_i!=j: continue
    k=0
    for i in range(2):
      for ii in range(2):
        log_filename = "%s_%d_%d.log" % (pdb, i, ii)
        cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
        cmd += " link_all=%d link_carbohydrate=%d" % (i, ii)
        print(cmd)
        run_and_test(cmd, pdb,k)
        k+=1


if __name__=="__main__":
  import sys
  run(*tuple(sys.argv[1:]))
