import os, sys

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
        }

def run():
  cifs = ""
  for pdb in pdbs:
    f=file(pdb, "wb")
    f.write(pdbs[pdb])
    f.close()
    if pdb.endswith(".cif"): cifs += " %s" % pdb
  for pdb in pdbs:
    if pdb.endswith(".cif"): continue
    for i in range(2):
      log_filename = "%s_%d.log" % (pdb, i)
      cmd = "phenix.pdbtools %s" % pdb
      cmd += " intra=%d %s > %s" % (i, cifs, log_filename)
      print cmd
      easy_run.call(cmd)
      f=file(log_filename, "rb")
      lines = f.read()
      assert lines.find("Output model")>-1
      print "OK"

if __name__=="__main__":
  run()
