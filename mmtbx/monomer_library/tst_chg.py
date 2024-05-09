
from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
from six.moves import cStringIO as StringIO

l1r_pdb = """
HETATM 3890  C1 AL1R A1247      21.777  10.761  29.793  0.50  4.10           C
HETATM 3891  C2 AL1R A1247      21.624  12.183  29.278  0.50 23.99           C
HETATM 3892  C3 AL1R A1247      20.238  12.698  29.622  0.50 30.37           C
HETATM 3893  O4 AL1R A1247      22.518  13.081  29.986  0.50 26.04           O
HETATM 3894  C5 AL1R A1247      23.787  12.666  30.304  0.50 24.33           C
HETATM 3895  C6 AL1R A1247      24.205  13.396  31.411  0.50 20.37           C
HETATM 3896  C7 AL1R A1247      25.448  13.174  31.967  0.50 22.78           C
HETATM 3897  C8 AL1R A1247      26.303  12.246  31.406  0.50 24.56           C
HETATM 3898  C9 AL1R A1247      25.914  11.515  30.291  0.50 23.51           C
HETATM 3899  N10AL1R A1247      26.803  10.567  29.731  0.50 19.01           N
HETATM 3900  C11AL1R A1247      27.021  10.367  28.421  0.50 17.27           C
HETATM 3901  O12AL1R A1247      26.512  11.039  27.552  0.50 19.42           O
HETATM 3902  C13AL1R A1247      27.975   9.276  28.052  0.50  8.77           C
HETATM 3903  C14AL1R A1247      28.790   9.523  26.954  0.50  9.85           C
HETATM 3904  C15AL1R A1247      29.707   8.574  26.543  0.50  8.72           C
HETATM 3905 CL1 AL1R A1247      30.719   8.896  25.168  0.50 15.30          CL
HETATM 3906  C17AL1R A1247      29.805   7.356  27.241  0.50 12.61           C
HETATM 3907  O18AL1R A1247      30.698   6.376  26.912  0.50 12.55           O
HETATM 3908  C19AL1R A1247      31.814   6.192  27.805  0.50 17.96           C
HETATM 3909  C20AL1R A1247      32.392   4.807  27.590  0.50 16.10           C
HETATM 3910  N21AL1R A1247      33.124   4.698  26.299  0.50 10.39           N
HETATM 3911  C22AL1R A1247      28.996   7.117  28.347  0.50 14.87           C
HETATM 3912 CL2 AL1R A1247      29.141   5.627  29.220  0.50 17.63          CL
HETATM 3913  C24AL1R A1247      28.081   8.079  28.757  0.50  6.06           C
HETATM 3914  C25AL1R A1247      24.657  11.738  29.738  0.50 26.75           C
HETATM 3915 H1C1AL1R A1247      23.622  13.862  34.406  0.50  4.10           H
HETATM 3916 H1C2AL1R A1247      23.741  15.372  34.937  0.50  4.10           H
HETATM 3917 H1C3AL1R A1247      24.624  14.213  35.610  0.50  4.10           H
HETATM 3918  H2 AL1R A1247      24.984  15.262  32.978  0.50 23.99           H
HETATM 3919 H3C1AL1R A1247      26.894  15.122  35.102  0.50 30.37           H
HETATM 3920 H3C2AL1R A1247      27.151  15.812  33.676  0.50 30.37           H
HETATM 3921 H3C3AL1R A1247      26.119  16.472  34.712  0.50 30.37           H
HETATM 3922  H6 AL1R A1247      24.037  14.353  31.774  0.50 20.37           H
HETATM 3923  H7 AL1R A1247      23.526  13.546  29.621  0.50 22.78           H
HETATM 3924  H25AL1R A1247      27.038  11.581  32.111  0.50 26.75           H
HETATM 3925  H8 AL1R A1247      24.760  11.752  28.709  0.50 24.56           H
HETATM 3926  H10AL1R A1247      26.981   9.745  30.294  0.50 19.01           H
HETATM 3927  H14AL1R A1247      29.071  10.449  26.709  0.50  9.85           H
HETATM 3928  H24AL1R A1247      27.305   7.828  29.328  0.50  6.06           H
HETATM 3929 H191AL1R A1247      32.651   6.813  27.086  0.50 17.96           H
HETATM 3930 H192AL1R A1247      31.921   6.649  28.507  0.50 17.96           H
HETATM 3931 H201AL1R A1247      31.734   4.274  27.962  0.50 16.10           H
HETATM 3932 H202AL1R A1247      33.194   4.782  28.360  0.50 16.10           H
HETATM 3933 H211AL1R A1247      33.935   4.283  26.431  0.50 10.39           H
HETATM 3934 H212AL1R A1247      32.687   3.476  26.242  0.50 10.39           H
HETATM 3935  C1 BL1R A1247      24.235  14.539  34.768  0.50 24.67           C
HETATM 3936  C2 BL1R A1247      25.347  14.803  33.765  0.50 20.23           C
HETATM 3937  C3 BL1R A1247      26.478  15.625  34.367  0.50 10.96           C
HETATM 3938  O4 BL1R A1247      25.904  13.526  33.376  0.50 18.91           O
HETATM 3939  C5 BL1R A1247      25.578  13.036  32.139  0.50 15.56           C
HETATM 3940  C6 BL1R A1247      24.554  13.612  31.397  0.50 15.71           C
HETATM 3941  C7 BL1R A1247      24.251  13.135  30.140  0.50 13.06           C
HETATM 3942  C8 BL1R A1247      24.978  12.089  29.604  0.50 11.07           C
HETATM 3943  C9 BL1R A1247      26.005  11.502  30.330  0.50 14.33           C
HETATM 3944  N10BL1R A1247      26.744  10.445  29.752  0.50 18.41           N
HETATM 3945  C11BL1R A1247      27.125  10.393  28.461  0.50 22.57           C
HETATM 3946  O12BL1R A1247      26.764  11.204  27.639  0.50 21.58           O
HETATM 3947  C13BL1R A1247      28.054   9.286  28.075  0.50 23.49           C
HETATM 3948  C14BL1R A1247      29.022   9.556  27.113  0.50 24.22           C
HETATM 3949  C15BL1R A1247      29.919   8.573  26.720  0.50 22.20           C
HETATM 3950 CL1 BL1R A1247      31.123   8.934  25.520  0.50 22.46          CL
HETATM 3951  C17BL1R A1247      29.846   7.299  27.295  0.50 21.60           C
HETATM 3952  O18BL1R A1247      30.717   6.303  26.941  0.50 17.30           O
HETATM 3953  C19BL1R A1247      32.007   6.279  27.601  0.50 18.57           C
HETATM 3954  C20BL1R A1247      32.486   4.840  27.681  0.50 17.56           C
HETATM 3955  N21BL1R A1247      33.019   4.318  26.392  0.50 14.88           N
HETATM 3956  C22BL1R A1247      28.877   7.029  28.259  0.50 22.09           C
HETATM 3957 CL2 BL1R A1247      28.784   5.457  28.974  0.50 17.57          CL
HETATM 3958  C24BL1R A1247      27.986   8.021  28.649  0.50 22.43           C
HETATM 3959  C25BL1R A1247      26.313  11.991  31.592  0.50 14.30           C
HETATM 3960 H1C1BL1R A1247      23.622  13.862  34.406  0.50 24.67           H
HETATM 3961 H1C2BL1R A1247      23.741  15.372  34.937  0.50 24.67           H
HETATM 3962 H1C3BL1R A1247      24.624  14.213  35.610  0.50 24.67           H
HETATM 3963  H2 BL1R A1247      24.984  15.262  32.978  0.50 20.23           H
HETATM 3964 H3C1BL1R A1247      26.894  15.122  35.102  0.50 10.96           H
HETATM 3965 H3C2BL1R A1247      27.151  15.812  33.676  0.50 10.96           H
HETATM 3966 H3C3BL1R A1247      26.119  16.472  34.712  0.50 10.96           H
HETATM 3967  H6 BL1R A1247      24.037  14.353  31.774  0.50 15.71           H
HETATM 3968  H7 BL1R A1247      23.526  13.546  29.621  0.50 13.06           H
HETATM 3969  H25BL1R A1247      27.038  11.581  32.111  0.50 14.30           H
HETATM 3970  H8 BL1R A1247      24.760  11.752  28.709  0.50 11.07           H
HETATM 3971  H10BL1R A1247      26.981   9.745  30.294  0.50 18.41           H
HETATM 3972  H14BL1R A1247      29.071  10.449  26.709  0.50 24.22           H
HETATM 3973  H24BL1R A1247      27.305   7.828  29.328  0.50 22.43           H
HETATM 3974 H191BL1R A1247      32.651   6.813  27.086  0.50 18.57           H
HETATM 3975 H192BL1R A1247      31.921   6.649  28.507  0.50 18.57           H
HETATM 3976 H201BL1R A1247      31.734   4.274  27.962  0.50 17.56           H
HETATM 3977 H202BL1R A1247      33.194   4.782  28.360  0.50 17.56           H
HETATM 3978 H211BL1R A1247      33.935   4.283  26.431  0.50 14.88           H
HETATM 3979 H212BL1R A1247      32.687   3.476  26.242  0.50 14.88           H
"""
l1r_cif = """
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
L1R        L1R 'Unknown                  ' ligand 45 25 .
#
data_comp_L1R
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
L1R         C1     C   CH3   .          2.8974    3.5183   -6.1638
L1R         C2     C   CH1   .          2.8051    2.8353   -4.8033
L1R         C3     C   CH3   .          3.0983    3.8453   -3.7109
L1R         O4     O   O2    .          1.4865    2.3109   -4.6215
L1R         C5     C   CR6   .          1.4498    1.0843   -3.9638
L1R         C6     C   CR16  .          1.9122   -0.0757   -4.6162
L1R         C7     C   CR16  .          1.8001   -1.3045   -4.0047
L1R         C8     C   CR16  .          1.2312   -1.4031   -2.7447
L1R         C9     C   CR6   .          0.7688   -0.2476   -2.0911
L1R         N10    N   NH1   .          0.0935   -0.3502   -0.8144
L1R         C11    C   C     .          0.2928   -1.4890    0.0136
L1R         O12    O   O     .          1.0040   -2.3511   -0.3340
L1R         C13    C   CR6   .         -0.4260   -1.5883    1.3457
L1R         C14    C   CR16  .         -1.4607   -2.5265    1.5148
L1R         C15    C   CR6   .         -2.0425   -2.6929    2.7524
L1R        CL1     CL  CL    .         -3.3279   -3.8717    2.9672
L1R         C17    C   CR6   .         -1.6068   -1.9326    3.8358
L1R         O18    O   O2    .         -2.1493   -2.1530    5.1031
L1R         C19    C   CH2   .         -3.2986   -1.3409    5.3928
L1R         C20    C   CH2   .         -3.8899   -1.7593    6.7556
L1R         N21    N   NH2   .         -4.0324   -0.5996    7.5943
L1R         C22    C   CR6   .         -0.5745   -0.9968    3.6662
L1R        CL2     CL  CL    .          0.0312   -0.0872    5.0466
L1R         C24    C   CR16  .          0.0084   -0.8313    2.4222
L1R         C25    C   CR16  .          0.8832    0.9906   -2.7102
L1R        H1C1    H   HCH3  .          3.7452    3.3549   -6.5396
L1R        H1C2    H   HCH3  .          2.7708    4.4782   -6.0540
L1R        H1C3    H   HCH3  .          2.1773    3.1498   -6.7790
L1R         H2     H   HCH1  .          3.4506    2.1144   -4.7592
L1R        H3C1    H   HCH3  .          2.2868    4.1686   -3.3576
L1R        H3C2    H   HCH3  .          3.6009    4.5753   -4.0769
L1R        H3C3    H   HCH3  .          3.6107    3.4212   -3.0058
L1R         H6     H   HCR6  .          2.2983   -0.0118   -5.4727
L1R         H7     H   HCR6  .          2.1289   -2.1229   -4.4694
L1R         H25    H   HCR6  .          0.5581    1.8029   -2.2527
L1R         H8     H   HCR6  .          1.1516   -2.2771   -2.3119
L1R         H10    H   HNH1  .         -0.5730    0.3075   -0.5714
L1R         H14    H   HCR6  .         -1.6952   -3.1409    0.7736
L1R         H24    H   HCR6  .          0.7227   -0.1851    2.3037
L1R        H191    H   HCH2  .         -3.0377   -0.4105    5.4269
L1R        H192    H   HCH2  .         -3.9800   -1.4670    4.6814
L1R        H201    H   HCH2  .         -4.7783   -2.1791    6.6150
L1R        H202    H   HCH2  .         -3.2999   -2.3908    7.1809
L1R        H211    H   HNH2  .         -3.2822   -0.0557    7.5009
L1R        H212    H   HNH2  .         -4.8097   -0.1153    7.3360
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
L1R   C1      C2    single        1.525 0.020
L1R   C2      C3    single        1.516 0.020
L1R   C2      O4    single        1.431 0.020
L1R   O4      C5    single        1.392 0.020
L1R   C5      C6    aromatic      1.409 0.020
L1R   C5      C25   aromatic      1.379 0.020
L1R   C6      C7    aromatic      1.377 0.020
L1R   C7      C8    aromatic      1.386 0.020
L1R   C7      H7    single        0.997 0.020
L1R   C8      C9    aromatic      1.406 0.020
L1R   C9      N10   single        1.448 0.020
L1R   C9      C25   aromatic      1.389 0.020
L1R   N10     C11   single        1.422 0.020
L1R   C11     O12   double        1.170 0.020
L1R   C11     C13   single        1.517 0.020
L1R   C13     C14   aromatic      1.407 0.020
L1R   C13     C24   aromatic      1.386 0.020
L1R   C14     C15   aromatic      1.378 0.020
L1R   C15    CL1    single        1.757 0.020
L1R   C15     C17   aromatic      1.393 0.020
L1R   C17     O18   single        1.396 0.020
L1R   C17     C22   aromatic      1.404 0.020
L1R   O18     C19   single        1.437 0.020
L1R   C19     C20   single        1.543 0.020
L1R   C20     N21   single        1.438 0.020
L1R   C22    CL2    single        1.761 0.020
L1R   C22     C24   aromatic      1.384 0.020
L1R   C1     H1C1   single        0.942 0.020
L1R   C1     H1C2   single        0.974 0.020
L1R   C1     H1C3   single        1.016 0.020
L1R   C2      H2    single        0.969 0.020
L1R   C3     H3C1   single        0.942 0.020
L1R   C3     H3C2   single        0.959 0.020
L1R   C3     H3C3   single        0.969 0.020
L1R   C6      H6    single        0.942 0.020
L1R   C25     H25   single        0.987 0.020
L1R   C8      H8    single        0.979 0.020
L1R   N10     H10   single        0.967 0.020
L1R   C14     H14   single        0.991 0.020
L1R   C24     H24   single        0.970 0.020
L1R   C19    H191   single        0.967 0.020
L1R   C19    H192   single        0.993 0.020
L1R   C20    H201   single        0.993 0.020
L1R   C20    H202   single        0.963 0.020
L1R   N21    H211   single        0.931 0.020
L1R   N21    H212   single        0.952 0.020
#
loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
L1R  H1C3     C1     H1C2         109.48 3.000
L1R  H1C3     C1     H1C1         109.48 3.000
L1R  H1C2     C1     H1C1         109.45 3.000
L1R  H1C3     C1      C2          109.55 3.000
L1R  H1C2     C1      C2          109.43 3.000
L1R  H1C1     C1      C2          109.44 3.000
L1R   H2      C2      O4          109.61 3.000
L1R   H2      C2      C3          109.52 3.000
L1R   O4      C2      C3          109.32 3.000
L1R   H2      C2      C1          109.48 3.000
L1R   O4      C2      C1          109.47 3.000
L1R   C3      C2      C1          109.43 3.000
L1R  H3C3     C3     H3C2         109.50 3.000
L1R  H3C3     C3     H3C1         109.43 3.000
L1R  H3C2     C3     H3C1         109.47 3.000
L1R  H3C3     C3      C2          109.56 3.000
L1R  H3C2     C3      C2          109.48 3.000
L1R  H3C1     C3      C2          109.40 3.000
L1R   C5      O4      C2          114.03 3.000
L1R   C25     C5      C6          119.99 3.000
L1R   C25     C5      O4          120.01 3.000
L1R   C6      C5      O4          119.87 3.000
L1R   H6      C6      C7          119.86 3.000
L1R   H6      C6      C5          119.99 3.000
L1R   C7      C6      C5          120.15 3.000
L1R   H7      C7      C8          120.05 3.000
L1R   H7      C7      C6          119.91 3.000
L1R   C8      C7      C6          120.04 3.000
L1R   H8      C8      C9          120.12 3.000
L1R   H8      C8      C7          119.93 3.000
L1R   C9      C8      C7          119.95 3.000
L1R   C25     C9      N10         119.64 3.000
L1R   C25     C9      C8          119.89 3.000
L1R   N10     C9      C8          120.34 3.000
L1R   H10     N10     C11         119.65 3.000
L1R   H10     N10     C9          119.65 3.000
L1R   C11     N10     C9          120.32 3.000
L1R   C13     C11     O12         120.04 3.000
L1R   C13     C11     N10         119.82 3.000
L1R   O12     C11     N10         120.14 3.000
L1R   C24     C13     C14         120.09 3.000
L1R   C24     C13     C11         119.86 3.000
L1R   C14     C13     C11         119.85 3.000
L1R   H14     C14     C15         119.81 3.000
L1R   H14     C14     C13         119.84 3.000
L1R   C15     C14     C13         119.94 3.000
L1R   C17     C15    CL1          119.98 3.000
L1R   C17     C15     C14         120.03 3.000
L1R  CL1      C15     C14         119.99 3.000
L1R   C22     C17     O18         120.05 3.000
L1R   C22     C17     C15         119.99 3.000
L1R   O18     C17     C15         119.88 3.000
L1R   C19     O18     C17         113.87 3.000
L1R  H192     C19    H191         109.42 3.000
L1R  H192     C19     C20         109.59 3.000
L1R  H191     C19     C20         109.46 3.000
L1R  H192     C19     O18         109.43 3.000
L1R  H191     C19     O18         109.58 3.000
L1R   C20     C19     O18         109.35 3.000
L1R  H202     C20    H201         109.48 3.000
L1R  H202     C20     N21         109.38 3.000
L1R  H201     C20     N21         109.56 3.000
L1R  H202     C20     C19         109.45 3.000
L1R  H201     C20     C19         109.42 3.000
L1R   N21     C20     C19         109.53 3.000
L1R  H212     N21    H211         109.48 3.000
L1R  H212     N21     C20         109.45 3.000
L1R  H211     N21     C20         109.43 3.000
L1R   C24     C22    CL2          119.88 3.000
L1R   C24     C22     C17         119.88 3.000
L1R  CL2      C22     C17         120.18 3.000
L1R   H24     C24     C22         119.96 3.000
L1R   H24     C24     C13         119.97 3.000
L1R   C22     C24     C13         120.07 3.000
L1R   H25     C25     C9          119.98 3.000
L1R   H25     C25     C5          120.05 3.000
L1R   C9      C25     C5          119.97 3.000
"""

def run():
  if (not libtbx.env.has_module("phenix")):
    print("phenix not configured, skipping")
    return
  elif (not libtbx.env.find_in_repositories("chem_data")):
    print("chem_data not configured, skipping")
    return
  f=open("l1r.pdb", "w")
  f.write(l1r_pdb)
  f.close()
  f=open("l1r.cif", "w")
  f.write(l1r_cif)
  f.close()

  cmd = 'phenix.geometry_minimization'
  cmd += " l1r.pdb l1r.cif use_neutron_distances=False"
  cmd += ' selection="element H or element D" write_geo_file=False'
  cmd += ' correct_hydrogens=True'
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO()
  ero.show_stdout(out=err)

  cmd = "phenix.fmodel high_res=4.5 format=mtz label=FOBS type=real r_free=0.1 l1r.pdb generate_fake_p1_symmetry=1"
  assert not easy_run.call(cmd)
  cmd = 'phenix.refine'
  cmd += " l1r.pdb l1r.pdb.mtz l1r.cif"
  cmd += " main.number_of_macro_cycles=1"
  cmd += ' correct_hydrogens=True'
  print(cmd)
  ero = easy_run.fully_buffered(command=cmd)
  err = StringIO()
  ero.show_stdout(out=err)
  print("OK")


if __name__=="__main__":
  run()
