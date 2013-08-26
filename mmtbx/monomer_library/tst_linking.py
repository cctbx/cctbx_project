from __future__ import division
import os
from libtbx import easy_run

pdbs = {"linking_test_NAG-NAG.pdb" : """
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
HETATM 4156  O19 VSP B 600     -13.564  46.514 -32.785  1.00 21.41           O  """,
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
#
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
ATOM   5294  N   MET A 680      23.170  19.173  20.583  1.00 27.07           N
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
  "linking_test_3g2j-LLP.pdb" : [43, 43],
  }

def run():
  cifs = ""
  for pdb in pdbs:
    f=file(pdb, "wb")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
  for pdb in sorted(pdbs):
    if pdb.endswith(".cif"): continue
    #if pdb.find("LLP")>-1: continue
    for i in range(2):
      log_filename = "%s_%d.log" % (pdb, i)
      cmd = "phenix.geometry_minimization %s write_geo_file=True" % pdb
      cmd += " intra_chain=%d %s > %s" % (i, cifs, log_filename)
      print cmd
      easy_run.call(cmd)
      f=file(log_filename, "rb")
      lines = f.read()
      assert lines.find("Write PDB file")>-1
      print "OK"
      f=file(pdb.replace(".pdb", "_minimized.geo"), "rb")
      lines = f.readlines()
      f.close()
      for line in lines:
        if line.find("Bond restraints:")>-1:
          bonds = int(line.split()[2])
          break
      print bonds
      assert bonds == links[pdb][i]
      new_geo = pdb.replace(".pdb", "_minimized_%d.geo" % i)
      if (os.path.isfile(new_geo)) :
        os.remove(new_geo)
      os.rename(pdb.replace(".pdb", "_minimized.geo"), new_geo)
      print "OK"

if __name__=="__main__":
  run()
