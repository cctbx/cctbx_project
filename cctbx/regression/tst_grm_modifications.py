from __future__ import division
from __future__ import print_function
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cctbx.array_family import flex
import mmtbx.model
from StringIO import StringIO
#import libtbx.load_env
from cctbx import geometry_restraints
from libtbx.test_utils import show_diff
from libtbx.utils import null_out
import iotbx

raw_records1 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.016447  0.009496  0.000000        0.00000
SCALE2      0.000000  0.018992  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010309        0.00000
ATOM   1050  N   LYS A 135      31.992  14.930  -7.233  1.00  9.47           N
ATOM   1051  CA  LYS A 135      31.388  16.216  -7.637  1.00 12.89           C
ATOM   1052  C   LYS A 135      30.807  16.840  -6.406  1.00  6.47           C
ATOM   1053  O   LYS A 135      29.583  16.869  -6.191  1.00 15.74           O
ATOM   1054  CB  LYS A 135      30.263  16.059  -8.655  1.00 13.51           C
ATOM   1055  CG  LYS A 135      30.742  15.277  -9.843  1.00 16.23           C
ATOM   1056  CD  LYS A 135      29.612  15.131 -10.835  1.00 28.55           C
ATOM   1057  CE  LYS A 135      30.173  14.812 -12.216  1.00 34.52           C
ATOM   1058  NZ  LYS A 135      29.396  13.756 -12.899  1.00 46.18           N
TER    1294      LYS A 162
END

""".splitlines()

raw_records2 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
HETATM 1406  O   HOH A 282      32.366  19.942  24.727  1.00 38.09           O
END
""".splitlines()

# connect to O, distance 2.9A with symop -y+1,x-y,z-1/3
raw_records3 = """\
CRYST1   60.800   60.800   97.000  90.00  90.00 120.00 P 32 2 1      6
HETATM 1406  O   HOH A 282      32.366  19.942  24.727  1.00 38.09           O
HETATM 1407  O   HOH A 283      33.366  18.942  23.727  1.00 38.09           O
END
""".splitlines()

raw_records4 = """\
CRYST1   15.775   12.565   13.187  90.00  90.00  90.00 P 1
ATOM      1  N   MET A   1       9.821   6.568   5.000  1.00 66.07           N
ATOM      2  CA  MET A   1       9.946   7.171   6.357  1.00 66.55           C
ATOM      3  C   MET A   1      10.571   6.157   7.305  1.00 64.57           C
ATOM      4  O   MET A   1      10.775   5.000   6.933  1.00 66.25           O
ATOM      5  CB  MET A   1       8.570   7.565   6.897  1.00 69.08           C
ATOM      6  CG  MET A   1       7.723   6.373   7.299  1.00 71.37           C
ATOM      7  SD  MET A   1       6.247   6.862   8.187  1.00 76.22           S
ATOM      8  CE  MET A   1       5.000   6.694   6.892  1.00 74.93           C
END
""".splitlines()

raw_records5 = """\
CRYST1  258.687  258.687   47.103  90.00  90.00 120.00 P 63          6
ATOM    213  N   ILE A  78      87.236 -55.209   0.578  1.00179.51           N
ATOM    321  O   LYS A  93      81.801 -49.470  26.164  1.00197.87           O
"""

raw_records6 = """\
CRYST1  101.940  101.370  203.540  90.00  90.00  90.00 C 1 2 1
ATOM      1  N   GLU F 440      51.717  35.260  96.810  1.00173.09           N
ATOM      2  CA  GLU F 440      50.591  34.319  97.027  1.00158.07           C
ATOM      3  C   GLU F 440      49.903  34.410  98.438  1.00181.12           C
ATOM      4  O   GLU F 440      48.809  34.940  98.522  1.00178.50           O
ATOM      5  N   ILE F 441      50.585  34.020  99.528  1.00197.83           N
ATOM      6  CA  ILE F 441      51.029  34.946 100.496  1.00205.24           C
ATOM      7  C   ILE F 441      50.630  36.474 100.471  1.00215.52           C
ATOM      8  O   ILE F 441      50.805  37.162 101.476  1.00228.70           O
ATOM      9  N   GLY F 442      50.138  37.042  99.430  1.00196.69           N
ATOM     10  CA  GLY F 442      49.342  38.214  98.975  1.00156.00           C
ATOM     11  C   GLY F 442      48.218  38.275  97.877  1.00171.08           C
ATOM     12  O   GLY F 442      47.220  38.941  98.131  1.00175.64           O
ATOM     13  N   ILE F 443      48.268  37.635  96.688  1.00184.50           N
ATOM     14  CA  ILE F 443      47.069  37.544  95.818  1.00185.10           C
ATOM     15  C   ILE F 443      46.104  36.600  96.497  1.00180.57           C
ATOM     16  O   ILE F 443      44.887  36.627  96.262  1.00172.35           O
ATOM     17  N   LEU F 444      46.637  35.791  97.397  1.00186.46           N
ATOM     18  CA  LEU F 444      45.893  35.148  98.444  1.00184.47           C
ATOM     19  C   LEU F 444      45.942  35.880  99.806  1.00195.40           C
ATOM     20  O   LEU F 444      44.848  35.788 100.444  1.00193.90           O
"""
raw_records7 = """\
CRYST1  101.940  101.370  203.540  90.00  90.00  90.00 C 1 2 1
ATOM      1  N   ILE F 437      55.794  33.676  96.336  1.00160.31           N
ATOM      2  CA  ILE F 437      54.818  33.116  97.305  1.00161.41           C
ATOM      3  C   ILE F 437      54.482  34.092  98.461  1.00184.09           C
ATOM      4  O   ILE F 437      53.584  33.767  99.360  1.00201.20           O
ATOM      5  N   ASN F 438      55.255  35.246  98.488  1.00184.69           N
ATOM      6  CA  ASN F 438      55.094  36.686  98.818  1.00185.32           C
ATOM      7  C   ASN F 438      53.768  37.307  98.414  1.00190.11           C
ATOM      8  O   ASN F 438      52.882  37.698  99.302  1.00206.13           O
ATOM      9  N   SER F 439      53.683  37.260  97.042  1.00171.17           N
ATOM     10  CA  SER F 439      52.551  37.517  96.151  1.00180.82           C
ATOM     11  C   SER F 439      51.492  36.418  96.171  1.00176.83           C
ATOM     12  O   SER F 439      50.456  36.620  95.555  1.00163.02           O
ATOM     13  N   GLU F 440      51.717  35.260  96.810  1.00173.09           N
ATOM     14  CA  GLU F 440      50.591  34.319  97.027  1.00158.07           C
ATOM     15  C   GLU F 440      49.903  34.410  98.438  1.00181.12           C
ATOM     16  O   GLU F 440      48.809  34.940  98.522  1.00178.50           O
ATOM     17  N   ILE F 441      50.585  34.020  99.528  1.00197.83           N
ATOM     18  CA  ILE F 441      51.029  34.946 100.496  1.00205.24           C
ATOM     19  C   ILE F 441      50.630  36.474 100.471  1.00215.52           C
ATOM     20  O   ILE F 441      50.805  37.162 101.476  1.00228.70           O
"""
# excluded:
#ATOM      1  N   MET A   1       9.821   6.568   5.000  1.00 66.07           N

raw_records8 = """\
CRYST1   52.400   52.400   30.700  90.00  90.00 120.00 P 6           6
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.019084  0.011018  0.000000        0.00000
SCALE2      0.000000  0.022036  0.000000        0.00000
SCALE3      0.000000  0.000000  0.032573        0.00000
ATOM    159  N   GLU A  21     -19.011 -14.048  -5.068  1.00 30.55           N
ATOM    160  CA  GLU A  21     -19.766 -14.465  -3.891  1.00 30.92           C
ATOM    161  C   GLU A  21     -19.105 -15.666  -3.211  1.00 35.35           C
ATOM    162  O   GLU A  21     -19.038 -15.734  -1.983  1.00 35.38           O
ATOM    163  CB  GLU A  21     -21.195 -14.821  -4.305  1.00 32.92           C
ATOM    164  CG  GLU A  21     -22.217 -14.782  -3.174  1.00 46.59           C
ATOM    165  CD  GLU A  21     -23.578 -15.309  -3.606  1.00 71.69           C
ATOM    166  OE1 GLU A  21     -24.218 -16.017  -2.800  1.00 69.37           O
ATOM    167  OE2 GLU A  21     -24.009 -15.022  -4.747  1.00 69.89           O
ATOM    168  N   GLU A  22     -18.619 -16.598  -4.030  1.00 31.53           N
ATOM    169  CA  GLU A  22     -17.953 -17.814  -3.566  1.00 30.28           C
ATOM    170  C   GLU A  22     -16.654 -17.504  -2.822  1.00 31.94           C
ATOM    171  O   GLU A  22     -16.408 -18.048  -1.744  1.00 30.30           O
ATOM    172  CB  GLU A  22     -17.674 -18.742  -4.753  1.00 31.45           C
ATOM    173  CG  GLU A  22     -17.177 -20.134  -4.374  1.00 41.37           C
ATOM    174  CD  GLU A  22     -16.558 -20.880  -5.546  1.00 59.74           C
ATOM    175  OE1 GLU A  22     -17.028 -20.710  -6.694  1.00 42.11           O
ATOM    176  OE2 GLU A  22     -15.599 -21.646  -5.317  1.00 54.25           O
ATOM    177  N   LEU A  23     -15.839 -16.626  -3.409  1.00 28.44           N
ATOM    178  CA  LEU A  23     -14.552 -16.219  -2.841  1.00 27.53           C
ATOM    179  C   LEU A  23     -14.712 -15.658  -1.429  1.00 29.94           C
ATOM    180  O   LEU A  23     -13.955 -16.017  -0.523  1.00 28.82           O
ATOM    181  CB  LEU A  23     -13.875 -15.185  -3.755  1.00 28.11           C
ATOM    182  CG  LEU A  23     -12.371 -14.864  -3.672  1.00 33.41           C
ATOM    183  CD1 LEU A  23     -11.984 -14.133  -2.387  1.00 35.92           C
ATOM    184  CD2 LEU A  23     -11.515 -16.114  -3.866  1.00 34.00           C
ATOM    185  N   LEU A  24     -15.700 -14.781  -1.255  1.00 26.30           N
ATOM    186  CA  LEU A  24     -16.000 -14.189   0.047  1.00 25.29           C
ATOM    187  C   LEU A  24     -16.376 -15.244   1.082  1.00 29.65           C
ATOM    188  O   LEU A  24     -16.052 -15.103   2.261  1.00 29.77           O
ATOM    189  CB  LEU A  24     -17.118 -13.147  -0.071  1.00 25.33           C
ATOM    190  CG  LEU A  24     -16.839 -11.842  -0.825  1.00 29.60           C
ATOM    191  CD1 LEU A  24     -18.092 -10.984  -0.842  1.00 30.42           C
ATOM    192  CD2 LEU A  24     -15.675 -11.066  -0.212  1.00 30.06           C
ATOM    193  N   ARG A  25     -17.051 -16.298   0.629  1.00 25.80           N
ATOM    194  CA  ARG A  25     -17.470 -17.388   1.506  1.00 25.11           C
ATOM    195  C   ARG A  25     -16.307 -18.309   1.874  1.00 27.79           C
ATOM    196  O   ARG A  25     -16.268 -18.843   2.983  1.00 27.92           O
ATOM    197  CB  ARG A  25     -18.628 -18.170   0.882  1.00 23.23           C
ATOM    198  CG  ARG A  25     -19.892 -17.339   0.708  1.00 33.72           C
ATOM    199  CD  ARG A  25     -21.050 -18.167   0.193  1.00 39.36           C
ATOM    200  NE  ARG A  25     -22.246 -17.357  -0.034  1.00 49.67           N
ATOM    201  CZ  ARG A  25     -23.186 -17.118   0.878  1.00 67.61           C
ATOM    202  NH1 ARG A  25     -23.087 -17.626   2.100  1.00 57.30           N
ATOM    203  NH2 ARG A  25     -24.233 -16.367   0.565  1.00 54.80           N
ATOM    204  N   VAL A  26     -15.368 -18.482   0.942  1.00 24.42           N
ATOM    205  CA  VAL A  26     -14.107 -19.183   1.205  1.00 24.67           C
ATOM    206  C   VAL A  26     -13.336 -18.449   2.304  1.00 28.49           C
ATOM    207  O   VAL A  26     -12.730 -19.071   3.179  1.00 28.78           O
ATOM    208  CB  VAL A  26     -13.223 -19.272  -0.070  1.00 28.72           C
ATOM    209  CG1 VAL A  26     -11.867 -19.887   0.249  1.00 28.52           C
ATOM    210  CG2 VAL A  26     -13.925 -20.074  -1.157  1.00 28.70           C
ATOM    211  N   ALA A  27     -13.375 -17.120   2.242  1.00 24.71           N
ATOM    212  CA  ALA A  27     -12.739 -16.260   3.236  1.00 24.02           C
ATOM    213  C   ALA A  27     -13.451 -16.320   4.585  1.00 27.75           C
ATOM    214  O   ALA A  27     -12.811 -16.200   5.633  1.00 28.05           O
ATOM    215  CB  ALA A  27     -12.668 -14.827   2.731  1.00 24.38           C
ATOM    216  N   VAL A  28     -14.773 -16.497   4.550  1.00 22.56           N
ATOM    217  CA  VAL A  28     -15.562 -16.707   5.766  1.00 22.49           C
ATOM    218  C   VAL A  28     -15.168 -18.028   6.426  1.00 28.43           C
ATOM    219  O   VAL A  28     -14.905 -18.062   7.634  1.00 28.37           O
ATOM    220  CB  VAL A  28     -17.094 -16.637   5.504  1.00 25.59           C
ATOM    221  CG1 VAL A  28     -17.887 -17.115   6.721  1.00 25.73           C
ATOM    222  CG2 VAL A  28     -17.504 -15.223   5.145  1.00 24.61           C
ATOM    223  N  AASN A  29     -15.101 -19.081   5.664  0.73 26.11           N
ATOM    224  CA AASN A  29     -14.718 -20.307   6.284  0.73 27.31           C
ATOM    225  C  AASN A  29     -13.257 -20.393   6.647  0.73 32.18           C
ATOM    226  O  AASN A  29     -12.947 -21.066   7.577  0.73 32.90           O
ATOM    227  CB AASN A  29     -15.331 -21.509   5.631  0.73 28.51           C
ATOM    228  CG AASN A  29     -16.814 -21.654   6.039  0.73 51.96           C
ATOM    229  OD1AASN A  29     -17.167 -22.363   6.966  0.73 47.90           O
ATOM    230  ND2AASN A  29     -17.663 -20.915   5.365  0.73 42.10           N
ATOM    231  N  BASN A  29     -15.140 -19.071   5.632  0.27 25.87           N
ATOM    232  CA BASN A  29     -14.710 -20.410   5.973  0.27 26.86           C
ATOM    233  C  BASN A  29     -13.308 -20.386   6.582  0.27 31.86           C
ATOM    234  O  BASN A  29     -13.079 -21.020   7.567  0.27 32.37           O
ATOM    235  CB BASN A  29     -14.680 -21.206   4.683  0.27 28.33           C
ATOM    236  CG BASN A  29     -14.977 -22.655   4.875  0.27 51.98           C
ATOM    237  OD1BASN A  29     -15.533 -23.029   5.867  0.27 45.83           O
ATOM    238  ND2BASN A  29     -14.576 -23.486   3.918  0.27 44.66           N
ATOM    239  N   ALA A  30     -12.413 -19.597   5.972  1.00 28.45           N
ATOM    240  CA  ALA A  30     -10.994 -19.492   6.328  1.00 29.07           C
ATOM    241  C   ALA A  30     -10.801 -18.726   7.633  1.00 34.39           C
ATOM    242  O   ALA A  30      -9.817 -18.931   8.344  1.00 35.36           O
ATOM    243  CB  ALA A  30     -10.211 -18.827   5.204  1.00 29.27           C
END
"""

def make_initial_grm(mon_lib_srv, ener_lib, records):
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = records,
    force_symmetry = True)

  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = True,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  return geometry, xrs

def show_sorted_geometry(geometry, xrs, file_name):
  out = open(file_name,'w')
  geometry.show_sorted(
      sites_cart=xrs.sites_cart(),
      site_labels=xrs.scatterers().extract_labels(),
      f=out)
  out.close()

def exercise_add_new_bond_restraint_in_place(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records4)

  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(0,3),
    distance_ideal=2.0,
    weight=3000)
  assert not geometry.is_bonded_atoms(0,3)
  assert not geometry.is_bonded_atoms(3,0)
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart())
  assert geometry.is_bonded_atoms(0,3)
  assert geometry.is_bonded_atoms(3,0)
  assert geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  # That's the way to get them:
  simple, asu = geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 8
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0

def exercise_add_super_long_bond(mon_lib_srv, ener_lib):
  # distance between the two is 26A, they are not added because of
  # max_distance_between_connecting_atoms=5 is default.
  # Inspired by 4c8q, atoms are from offending SHEET with wrong parallel/
  # antiparallel definition.
  #
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records5)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(0,1),
    distance_ideal=2.0,
    weight=3000)
  assert not geometry.is_bonded_atoms(0,1)
  assert not geometry.is_bonded_atoms(1,0)
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=5)
  assert not geometry.is_bonded_atoms(0,1)
  assert not geometry.is_bonded_atoms(1,0)
  # !!! This will fail, but should not. Left for future investigation.
  # geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
  #     max_distance_between_connecting_atoms=30)

def exercise_bond_near_symmetry(mon_lib_srv, ener_lib):
  """ Making bond near symmetry mate:
  " N   LEU F 444 " -- " O   GLU F 440 " """
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records6)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(16,3),
    distance_ideal=2.9,
    weight=3000)
  assert not geometry.is_bonded_atoms(16,3)
  assert not geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'before_xercise_bond_near_symmetry.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 19
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=10)
  assert geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'after_xercise_bond_near_symmetry.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 20  # assert 20
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0  # assert 0

def exercise_bond_near_symmetry2(mon_lib_srv, ener_lib):
  """ Same as previous, other atoms were still failing
  Making bond near symmetry mate:
  " N   LEU F 441 " -- " O   GLU F 437 " """
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records7)
  proxy = geometry_restraints.bond_simple_proxy(
    i_seqs=(16,3),
    distance_ideal=2.9,
    weight=3000)
  assert not geometry.is_bonded_atoms(16,3)
  assert not geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'before_xercise_bond_near_symmetry2.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 19
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  geometry.add_new_bond_restraints_in_place([proxy], xrs.sites_cart(),
      max_distance_between_connecting_atoms=10)
  assert geometry.is_bonded_atoms(3,16)
  # show_sorted_geometry(geometry, xrs, 'after_xercise_bond_near_symmetry2.geo')
  assert geometry.pair_proxies().bond_proxies.simple.size() == 20  # assert 20
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0  # assert 0

def exercise_single_atom(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)

  # output for debugging!!!
  # show_sorted_geometry(geometry, xrs, 'before_exersice_single_atoms.geo')

  xrs_add = iotbx.pdb.input(source_info=None, lines=raw_records2) \
      .xray_structure_simple()
  proxy1 = geometry_restraints.bond_simple_proxy(
      i_seqs=(3,9),
      distance_ideal=2.0,
      weight=3000,
      origin_id=2) # just a test

  new_xrs = xrs.concatenate(xrs_add)
  all_sites_cart = new_xrs.sites_cart()
  number_of_new_atoms = len(xrs_add.sites_cart())
  new_geometry = geometry.new_included_bonded_atoms(
      proxies=[proxy1],
      sites_cart=all_sites_cart,
      site_symmetry_table=xrs_add.site_symmetry_table(),
      nonbonded_types=flex.std_string(["OH2"]*number_of_new_atoms),
      nonbonded_charges=flex.int(number_of_new_atoms, 0),
      skip_max_proxy_distance_calculation=True)
  # output for debugging!!!
  # show_sorted_geometry(new_geometry, new_xrs,
  #     'after_exersice_single_atoms.geo')
  assert new_geometry.is_bonded_atoms(3,9)

  assert new_geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert new_geometry.pair_proxies().bond_proxies.asu.size() == 1
  # That's the way to get them:
  simple, asu = new_geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8
  simple, asu = new_geometry.get_all_bond_proxies()
  assert simple.size() + asu.size() == 9, "%d, %d" % (simple.size(), asu.size())
  assert new_geometry.pair_proxies().nonbonded_proxies.simple.size() == 10
  assert new_geometry.pair_proxies().nonbonded_proxies.asu.size() ==2
  assert new_geometry.get_hbond_proxies_iseqs() == [(3, 9)]
  simple, asu = new_geometry.get_covalent_bond_proxies()
  assert simple.size() + asu.size() == 8, "%d" % (simple.size() + asu.size())

def exercise_multiple_atoms(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)

  # output for debugging!!!
  # show_sorted_geometry(geometry, xrs, 'before_exersice_multiple_atoms.geo')

  xrs_add = iotbx.pdb.input(source_info=None, lines=raw_records3) \
      .xray_structure_simple()
  proxy1 = geometry_restraints.bond_simple_proxy(
      i_seqs=(3,9),
      distance_ideal=2.0,
      weight=3000)
  proxy2 = geometry_restraints.bond_simple_proxy(
      i_seqs=(4,10),
      distance_ideal=2.0,
      weight=3000)

  new_xrs = xrs.concatenate(xrs_add)
  all_sites_cart = new_xrs.sites_cart()
  number_of_new_atoms = len(xrs_add.sites_cart())
  assert geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert geometry.pair_proxies().bond_proxies.asu.size() == 0
  assert geometry.pair_proxies().nonbonded_proxies.simple.size() == 10
  assert geometry.pair_proxies().nonbonded_proxies.asu.size() == 0
  new_geometry = geometry.new_included_bonded_atoms(
      proxies=[proxy1, proxy2],
      sites_cart=all_sites_cart,
      site_symmetry_table=xrs_add.site_symmetry_table(),
      nonbonded_types=flex.std_string(["OH2"]*number_of_new_atoms),
      nonbonded_charges=flex.int(number_of_new_atoms, 0),
      skip_max_proxy_distance_calculation=True)
  # output for debugging!!!
  # show_sorted_geometry(new_geometry, new_xrs,
  #     'after_exersice_multiple_atoms.geo')
  assert new_geometry.pair_proxies().bond_proxies.simple.size() == 8
  assert new_geometry.pair_proxies().bond_proxies.asu.size() == 2
  assert new_geometry.pair_proxies().nonbonded_proxies.simple.size() == 11
  assert new_geometry.pair_proxies().nonbonded_proxies.asu.size() == 4

def exercise_is_bonded_atoms(mon_lib_srv, ener_lib):
  geometry, xrs = make_initial_grm(mon_lib_srv, ener_lib, raw_records1)
  # show_sorted_geometry(geometry, xrs, 'blabla.geo')
  linked_atoms = [(0,1), (1,2), (1,4), (2,3), (4,5),(5,6), (6,7), (7,8)]
  for i in range (9):
    for j in range(9):
      assert geometry.is_bonded_atoms(i,j) == geometry.is_bonded_atoms(j,i)
      assert geometry.is_bonded_atoms(i,j) == (tuple(sorted([i,j])) in linked_atoms)

def exercise_bond_near_symmetry3(mon_lib_srv, ener_lib):
  """ Since neighbors_fast_pair_generator for non-symmetry interactions
  provides only (i,j) pair and not (j,i), and there's no sorting involved
  (no way to guess what to check), there was a bug where non-symmetry interaction
  was missed in add_new_bond_restraints_in_place.
  Actually testing that both bonds are added without symmetry operators.
  """

  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=raw_records8)
  model = mmtbx.model.manager(model_input = pdb_inp, log=null_out())
  grm = model.get_restraints_manager().geometry
  h = model.get_hierarchy()
  proxy = geometry_restraints.bond_simple_proxy(
      i_seqs=(64,37),
      distance_ideal=2.9,
      weight=400,
      origin_id=origin_ids.get_origin_id('hydrogen bonds'))
  proxy2 = geometry_restraints.bond_simple_proxy(
      i_seqs=(72,37),
      distance_ideal=2.9,
      weight=400,
      origin_id=origin_ids.get_origin_id('hydrogen bonds'))
  grm.add_new_hbond_restraints_in_place(
      proxies=[proxy, proxy2],
      sites_cart=h.atoms().extract_xyz(),
      max_distance_between_connecting_atoms=10)
  sites_cart = h.atoms().extract_xyz()
  site_labels = model.get_xray_structure().scatterers().extract_labels()
  pair_proxies = grm.pair_proxies(flags=None, sites_cart=sites_cart)

  out = StringIO()
  pair_proxies.bond_proxies.show_sorted(
      by_value="residual",
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=out,
      prefix="",
      origin_id=origin_ids.get_origin_id('hydrogen bonds'))
  outtxt = out.getvalue()
  assert not show_diff(outtxt,
"""\
Bond restraints: 2
Sorted by residual:
bond pdb=" O   ARG A  25 "
     pdb=" N  AASN A  29 "
  ideal  model  delta    sigma   weight residual
  2.900  2.934 -0.034 5.00e-02 4.00e+02 4.53e-01
bond pdb=" O   ARG A  25 "
     pdb=" N  BASN A  29 "
  ideal  model  delta    sigma   weight residual
  2.900  2.888  0.012 5.00e-02 4.00e+02 5.59e-02
"""
)

def exercise():
  mon_lib_srv = None
  ener_lib = None
  try:
    mon_lib_srv = monomer_library.server.server()
    ener_lib = monomer_library.server.ener_lib()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None and ener_lib is not None:
    exercise_add_new_bond_restraint_in_place(mon_lib_srv, ener_lib)
    exercise_single_atom(mon_lib_srv, ener_lib)
    exercise_multiple_atoms(mon_lib_srv, ener_lib)
    exercise_is_bonded_atoms(mon_lib_srv, ener_lib)
    exercise_add_super_long_bond(mon_lib_srv, ener_lib)
    exercise_bond_near_symmetry(mon_lib_srv, ener_lib)
    exercise_bond_near_symmetry2(mon_lib_srv, ener_lib)
    exercise_bond_near_symmetry3(mon_lib_srv, ener_lib)

if (__name__ == "__main__"):
  exercise()
