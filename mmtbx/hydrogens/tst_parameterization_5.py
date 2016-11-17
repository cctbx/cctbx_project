from __future__ import division
import time

import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from cctbx import geometry_restraints
from mmtbx.hydrogens import riding_h
from mmtbx.hydrogens import parameterization

#-----------------------------------------------------------------------------
# This test checks the parameterization of H atoms for multiple conformations
# These examples are from pdb and initially failed
#-----------------------------------------------------------------------------

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = None,
    raw_records    = pdb_str,
    force_symmetry = True)
  pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xray_structure = processed_pdb_file.xray_structure()

  geometry_restraints = processed_pdb_file.geometry_restraints_manager(
    show_energies = False)

  sites_cart = xray_structure.sites_cart()
  atoms = pdb_hierarchy.atoms()

  riding_h_manager = riding_h.create_riding_h_manager(
    hierarchy           = pdb_hierarchy,
    geometry_restraints = geometry_restraints,
    crystal_symmetry    = xray_structure.crystal_symmetry())

  h_parameterization = riding_h_manager.h_parameterization
  diagnostics = parameterization.diagnostics_parameterization(
    connectivity_obj   = riding_h_manager.connectivity_obj,
    h_parameterization = h_parameterization,
    sites_cart         = sites_cart,
    threshold          = 0.05)

  h_distances        = diagnostics.h_distances
  unk_list           = diagnostics.unk_list

# There are 53 H atoms in the pdb_string, check if all of them are recognized
# one H atom is refused (VAL 7 HA), as it is bound to two CA atoms at once)
  assert (len(h_parameterization.keys()) == 52), 'Not all H atoms are parameterized'

  type_list = []
  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    hp = h_parameterization[ih]
    type_list.append(hp.htype)
    assert (h_distances[ih] < 0.2), 'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (hp.htype, atoms[ih].name, ih, labels.resseq.strip())

  assert(len(unk_list) == 0), 'Some H atoms are not recognized'

  for type1, type2 in zip(type_list, type_list_known):
    assert (type1 == type2)

# Ideal amino acids
pdb_str = """\
CRYST1   27.805   30.931   25.453  90.00  90.00  90.00 P 1
SCALE1      0.035965  0.000000  0.000000        0.00000
SCALE2      0.000000  0.032330  0.000000        0.00000
SCALE3      0.000000  0.000000  0.039288        0.00000
ATOM      1  CB  SER A   1      10.251  13.658  12.560  1.00  6.41           C
ANISOU    1  CB  SER A   1      638    591   1207   -162    148    106       C
ATOM      2  OG  SER A   1       9.057  14.277  12.222  1.00  7.60           O
ANISOU    2  OG  SER A   1      680    955   1251   -244   -164    339       O
ATOM      3  HB2 SER A   1      10.809  13.573  11.771  1.00  6.41           H
ATOM      4  HB3 SER A   1      10.073  12.759  12.879  1.00  6.41           H
ATOM      5  HG  SER A   1       8.756  13.939  11.514  1.00  7.60           H
ATOM      6  N  ASER A   1      11.430  15.721  13.108  0.40  5.10           N
ANISOU    6  N  ASER A   1      500    632    808   -107     58    104       N
ATOM      7  CA ASER A   1      11.003  14.453  13.643  0.40  5.95           C
ANISOU    7  CA ASER A   1      582    671   1007   -168    135    158       C
ATOM      8  C  ASER A   1      12.173  13.691  14.274  0.40  5.79           C
ANISOU    8  C  ASER A   1      661    588    952   -116    206    135       C
ATOM      9  O  ASER A   1      12.067  13.130  15.360  0.40  6.83           O
ANISOU    9  O  ASER A   1      874    705   1015    -39    314    258       O
ATOM     10  HA ASER A   1      10.323  14.627  14.313  0.40  5.95           H
ATOM     11  N  BSER A   1      11.530  15.721  13.108  0.60  5.10           N
ANISOU   11  N  BSER A   1      500    632    808   -107     58    104       N
ATOM     12  CA BSER A   1      11.103  14.453  13.643  0.60  5.95           C
ANISOU   12  CA BSER A   1      582    671   1007   -168    135    158       C
ATOM     13  C  BSER A   1      12.273  13.691  14.274  0.60  5.79           C
ANISOU   13  C  BSER A   1      661    588    952   -116    206    135       C
ATOM     14  O  BSER A   1      12.167  13.130  15.360  0.60  6.83           O
ANISOU   14  O  BSER A   1      874    705   1015    -39    314    258       O
ATOM     15  HA BSER A   1      10.423  14.627  14.313  0.60  5.95           H
ATOM     16  CA  SER A   2       7.228  18.726  13.388  1.00  5.95           C
ANISOU   16  CA  SER A   2      582    671   1007   -168    135    158       C
ATOM     17  OG  SER A   2       5.298  18.567  11.943  1.00  7.60           O
ANISOU   17  OG  SER A   2      680    955   1251   -244   -164    339       O
ATOM     18  HA  SER A   2       6.548  18.900  14.058  1.00  5.95           H
ATOM     19  HG  SER A   2       5.000  18.226  11.236  1.00  7.60           H
ATOM     20  N  ASER A   2       7.555  19.994  12.853  0.40  5.10           N
ANISOU   20  N  ASER A   2      500    632    808   -107     58    104       N
ATOM     21  C  ASER A   2       8.298  17.964  14.019  0.40  5.79           C
ANISOU   21  C  ASER A   2      661    588    952   -116    206    135       C
ATOM     22  O  ASER A   2       8.192  17.403  15.105  0.40  6.83           O
ANISOU   22  O  ASER A   2      874    705   1015    -39    314    258       O
ATOM     23  CB ASER A   2       6.376  17.931  12.305  0.40  6.41           C
ANISOU   23  CB ASER A   2      638    591   1207   -162    148    106       C
ATOM     24  HB2ASER A   2       6.934  17.846  11.516  0.40  6.41           H
ATOM     25  HB3ASER A   2       6.198  17.032  12.624  0.40  6.41           H
ATOM     26  N  BSER A   2       7.655  19.994  12.853  0.60  5.10           N
ANISOU   26  N  BSER A   2      500    632    808   -107     58    104       N
ATOM     27  C  BSER A   2       8.398  17.964  14.019  0.60  5.79           C
ANISOU   27  C  BSER A   2      661    588    952   -116    206    135       C
ATOM     28  O  BSER A   2       8.292  17.403  15.105  0.60  6.83           O
ANISOU   28  O  BSER A   2      874    705   1015    -39    314    258       O
ATOM     29  CB BSER A   2       6.476  17.931  12.305  0.60  6.41           C
ANISOU   29  CB BSER A   2      638    591   1207   -162    148    106       C
ATOM     30  HB2BSER A   2       7.043  17.826  11.525  0.60  6.41           H
ATOM     31  HB3BSER A   2       6.278  17.039  12.633  0.60  6.41           H
ATOM     32  OG  SER A   3       8.321  21.302  15.330  1.00  7.60           O
ANISOU   32  OG  SER A   3      680    955   1251   -244   -164    339       O
ATOM     33  HG  SER A   3       8.020  20.964  14.622  1.00  7.60           H
ATOM     34  N  ASER A   3      10.694  22.746  16.216  0.40  5.10           N
ANISOU   34  N  ASER A   3      500    632    808   -107     58    104       N
ATOM     35  CA ASER A   3      10.267  21.478  16.751  0.40  5.95           C
ANISOU   35  CA ASER A   3      582    671   1007   -168    135    158       C
ATOM     36  C  ASER A   3      11.437  20.716  17.382  0.40  5.79           C
ANISOU   36  C  ASER A   3      661    588    952   -116    206    135       C
ATOM     37  O  ASER A   3      11.331  20.155  18.468  0.40  6.83           O
ANISOU   37  O  ASER A   3      874    705   1015    -39    314    258       O
ATOM     38  CB ASER A   3       9.515  20.683  15.668  0.40  6.41           C
ANISOU   38  CB ASER A   3      638    591   1207   -162    148    106       C
ATOM     39  HA ASER A   3       9.587  21.652  17.421  0.40  5.95           H
ATOM     40  HB2ASER A   3      10.073  20.598  14.879  0.40  6.41           H
ATOM     41  HB3ASER A   3       9.337  19.784  15.987  0.40  6.41           H
ATOM     42  N  BSER A   3      10.794  22.746  16.216  0.60  5.10           N
ANISOU   42  N  BSER A   3      500    632    808   -107     58    104       N
ATOM     43  CA BSER A   3      10.367  21.478  16.751  0.60  5.95           C
ANISOU   43  CA BSER A   3      582    671   1007   -168    135    158       C
ATOM     44  C  BSER A   3      11.537  20.716  17.382  0.60  5.79           C
ANISOU   44  C  BSER A   3      661    588    952   -116    206    135       C
ATOM     45  O  BSER A   3      11.431  20.155  18.468  0.60  6.83           O
ANISOU   45  O  BSER A   3      874    705   1015    -39    314    258       O
ATOM     46  CB BSER A   3       9.615  20.683  15.668  0.60  6.41           C
ANISOU   46  CB BSER A   3      638    591   1207   -162    148    106       C
ATOM     47  HA BSER A   3       9.687  21.652  17.421  0.60  5.95           H
ATOM     48  HB2BSER A   3      10.182  20.578  14.888  0.60  6.41           H
ATOM     49  HB3BSER A   3       9.417  19.791  15.996  0.60  6.41           H
ATOM     50  CA  TYR A   4      22.263   7.688   7.671  1.00 13.00           C
ANISOU   50  CA  TYR A   4     2514   1244   1180    -28    116    139       C
ATOM     51  C   TYR A   4      22.047   6.341   6.986  1.00 11.81           C
ANISOU   51  C   TYR A   4     2332   1026   1130   -105     88    190       C
ATOM     52  O   TYR A   4      21.511   5.414   7.583  1.00 12.59           O
ANISOU   52  O   TYR A   4     2489   1154   1140    -72    111    219       O
ATOM     53  N   VAL A   5      22.490   6.252   5.739  1.00 10.80           N
ANISOU   53  N   VAL A   5     2053    884   1167    -29     80    167       N
ATOM     54  H   VAL A   5      22.805   6.918   5.296  1.00 12.96           H
ATOM     55  CA AVAL A   5      22.489   5.000   5.000  0.40 10.37           C
ANISOU   55  CA AVAL A   5     1783    902   1255     67   -148    -14       C
ATOM     56  CA BVAL A   5      22.589   5.000   5.000  0.60 10.37           C
ANISOU   56  CA BVAL A   5     1783    902   1255     67   -148    -14       C
ATOM     57  N   SER A   6      17.731  16.841  12.524  1.00  4.20           N
ANISOU   57  N   SER A   6      578    469    547    -59     10   -119       N
ATOM     58  CA  SER A   6      17.665  17.225  11.128  1.00  4.34           C
ANISOU   58  CA  SER A   6      576    472    600     -2      0      0       C
ATOM     59  C   SER A   6      18.893  16.719  10.356  1.00  3.90           C
ANISOU   59  C   SER A   6      554    387    543    -50    -15    -34       C
ATOM     60  O   SER A   6      18.771  16.316   9.200  1.00  4.24           O
ANISOU   60  O   SER A   6      634    443    534    -77    -67    -59       O
ATOM     61  H   SER A   6      17.744  17.490  13.088  1.00  4.20           H
ATOM     62  HA  SER A   6      16.878  16.752  10.817  1.00  4.34           H
ATOM     63  CB ASER A   6      17.624  18.766  11.122  0.60  5.87           C
ANISOU   63  CB ASER A   6      843    568    821    229    -39    -32       C
ATOM     64  OG ASER A   6      17.467  19.197   9.817  0.60  6.95           O
ANISOU   64  OG ASER A   6     1107    702    833    257    200    132       O
ATOM     65  HB2ASER A   6      16.892  19.086  11.673  0.60  4.80           H
ATOM     66  HB3ASER A   6      18.442  19.126  11.500  0.60  4.80           H
ATOM     67  HG ASER A   6      16.940  19.851   9.799  0.60  4.95           H
ATOM     68  CB BSER A   6      17.207  18.629  10.754  0.20  4.89           C
ANISOU   68  CB BSER A   6      429    509    922     98    134    102       C
ATOM     69  OG BSER A   6      18.185  19.487  11.305  0.20  5.81           O
ANISOU   69  OG BSER A   6      757    477    971     72     15   -142       O
ATOM     70  HB2BSER A   6      17.148  18.735   9.792  0.20  4.80           H
ATOM     71  HB3BSER A   6      16.327  18.820  11.115  0.20  4.80           H
ATOM     72  HG BSER A   6      17.993  20.284  11.124  0.20  4.95           H
ATOM     73  CB CSER A   6      17.650  18.745  10.952  0.20  4.80           C
ANISOU   73  CB CSER A   6      432    322   1070   -132    173      4       C
ATOM     74  OG CSER A   6      16.465  19.192  11.576  0.20  4.95           O
ANISOU   74  OG CSER A   6      511    563    806    114    -22    -75       O
ATOM     75  HB2CSER A   6      18.432  19.150  11.359  0.20  4.80           H
ATOM     76  HB3CSER A   6      17.661  18.987  10.013  0.20  4.80           H
ATOM     77  HG CSER A   6      16.593  19.261  12.403  0.20  4.95           H
ATOM     78  N   VAL A   7      17.914  10.234  15.459  1.00 10.80           N
ANISOU   78  N   VAL A   7     2053    884   1167    -29     80    167       N
ATOM     79  C   VAL A   7      19.220   8.954  13.918  1.00  9.56           C
ANISOU   79  C   VAL A   7     1874    685   1073      3   -185     85       C
ATOM     80  O   VAL A   7      19.644   9.970  13.400  1.00 10.25           O
ANISOU   80  O   VAL A   7     2013    634   1247     25      7     95       O
ATOM     82  HA  VAL A   7      17.893   8.230  15.321  0.44 12.11           H
ATOM     83  CA AVAL A   7      17.913   8.982  14.720  0.44 10.37           C
ANISOU   83  CA AVAL A   7     1783    902   1255     67   -148    -14       C
ATOM     84  CB AVAL A   7      16.666   8.901  13.831  0.44 11.26           C
ANISOU   84  CB AVAL A   7     1571    963   1743     96   -176     62       C
ATOM     85  CG1AVAL A   7      16.663  10.050  12.901  0.44 11.06           C
ANISOU   85  CG1AVAL A   7     1704   1072   1425    125   -105     37       C
ATOM     86  CG2AVAL A   7      16.581   7.572  13.075  0.44 10.02           C
ANISOU   86  CG2AVAL A   7     1510    782   1516    -92     19    232       C
ATOM     87  HB AVAL A   7      15.878   8.972  14.392  0.44 13.51           H
ATOM     88 HG11AVAL A   7      16.650  10.872  13.416  0.44 13.27           H
ATOM     89 HG12AVAL A   7      15.875   9.999  12.339  0.44 13.27           H
ATOM     90 HG13AVAL A   7      17.463  10.015  12.354  0.44 13.27           H
ATOM     91 HG21AVAL A   7      15.762   7.118  13.327  0.44 12.03           H
ATOM     92 HG22AVAL A   7      17.348   7.026  13.309  0.44 12.03           H
ATOM     93 HG23AVAL A   7      16.581   7.751  12.121  0.44 12.03           H
ATOM     94  CA BVAL A   7      17.908   8.985  14.702  0.56 10.10           C
ANISOU   94  CA BVAL A   7     1914    779   1143    -16      5    185       C
ATOM     95  CB BVAL A   7      16.719   8.890  13.693  0.56 11.08           C
ANISOU   95  CB BVAL A   7     1662   1363   1184     84    -25     58       C
ATOM     96  CG1BVAL A   7      16.672   7.503  13.046  0.56 12.75           C
ANISOU   96  CG1BVAL A   7     1479   1600   1766    255   -607   -238       C
ATOM     97  CG2BVAL A   7      15.371   9.207  14.360  0.56 11.66           C
ANISOU   97  CG2BVAL A   7     1508   1829   1094    214   -319     28       C
ATOM     98  HB BVAL A   7      16.858   9.540  12.987  0.56 13.29           H
ATOM     99 HG11BVAL A   7      17.505   7.348  12.573  0.56 15.30           H
ATOM    100 HG12BVAL A   7      15.927   7.469  12.426  0.56 15.30           H
ATOM    101 HG13BVAL A   7      16.556   6.835  13.740  0.56 15.30           H
ATOM    102 HG21BVAL A   7      14.972   9.972  13.917  0.56 13.99           H
ATOM    103 HG22BVAL A   7      15.523   9.409  15.296  0.56 13.99           H
ATOM    104 HG23BVAL A   7      14.790   8.435  14.277  0.56 13.99           H
ATOM    105  N   SER A   8      14.499  15.875  18.201  1.00  5.10           N
ANISOU  105  N   SER A   8      500    632    808   -107     58    104       N
ATOM    106  CA  SER A   8      14.072  14.607  18.736  1.00  5.95           C
ANISOU  106  CA  SER A   8      582    671   1007   -168    135    158       C
ATOM    107  C   SER A   8      15.242  13.845  19.367  1.00  5.79           C
ANISOU  107  C   SER A   8      661    588    952   -116    206    135       C
ATOM    108  O   SER A   8      15.136  13.284  20.453  1.00  6.83           O
ANISOU  108  O   SER A   8      874    705   1015    -39    314    258       O
ATOM    110  HA  SER A   8      13.392  14.781  19.406  1.00  5.95           H
ATOM    111  CB ASER A   8      13.543  13.720  17.561  0.52 10.56           C
ANISOU  111  CB ASER A   8     1137    993   1883   -635   -478    170       C
ATOM    112  OG ASER A   8      13.113  12.490  18.030  0.52 13.20           O
ANISOU  112  OG ASER A   8     1612   1164   2238   -898    145    -29       O
ATOM    113  HB2ASER A   8      12.813  14.173  17.111  0.52  6.41           H
ATOM    114  HB3ASER A   8      14.245  13.593  16.904  0.52  6.41           H
ATOM    115  HG ASER A   8      12.320  12.554  18.300  0.52  7.60           H
ATOM    116  CB BSER A   8      13.320  13.812  17.653  0.47  6.41           C
ANISOU  116  CB BSER A   8      638    591   1207   -162    148    106       C
ATOM    117  OG BSER A   8      12.126  14.431  17.315  0.47  7.60           O
ANISOU  117  OG BSER A   8      680    955   1251   -244   -164    339       O
ATOM    118  HB2BSER A   8      13.878  13.727  16.864  0.47  6.41           H
ATOM    119  HB3BSER A   8      13.142  12.913  17.972  0.47  6.41           H
ATOM    120  HG BSER A   8      11.825  14.093  16.607  0.47  7.60           H
ATOM      5  C  AASN A   9      19.177   3.788  12.902  0.46 10.01           C
ANISOU    5  C  AASN A   9     1311   1084   1408    -15   -324    -40       C
ATOM      6  C  BASN A   9      19.096   3.776  12.973  0.54  9.67           C
ANISOU    6  C  BASN A   9     1289   1143   1241   -122   -106    160       C
ATOM     15  N   SER A  10      19.894   2.942  13.626  1.00  9.28           N
ANISOU   15  N   SER A  10     1274    902   1350    -21    -87    -74       N
ATOM     16  CA ASER A  10      19.300   1.835  14.341  0.63  9.94           C
ANISOU   16  CA ASER A  10     1406    774   1598     87   -335   -196       C
ATOM     17  CA BSER A  10      19.317   1.815  14.347  0.37  9.47           C
ANISOU   17  CA BSER A  10     1420    764   1413     89    320   -137       C
ATOM     22  H   SER A  10      20.749   2.997  13.694  1.00 11.14           H
TER
END
"""

type_list_known = ['2tetra', '2tetra', 'alg1b', '3neigbs', '3neigbs',
  '3neigbs', 'alg1b', '2tetra', '2tetra', '2tetra', '2tetra', 'alg1b',
  '3neigbs', '2tetra', '2tetra', '3neigbs', '2tetra', '2tetra',
  'flat_2neigbs', 'alg1b', '3neigbs', '2tetra', '2tetra', 'alg1b',
  '2tetra', '2tetra', 'alg1b', '2tetra', '2tetra', 'alg1b', '3neigbs',
  'prop', 'prop', 'prop', 'prop', 'prop', 'prop', '3neigbs', 'prop',
  'prop', 'prop', 'prop', 'prop', 'prop', '3neigbs', '2tetra', '2tetra',
  'alg1b', '2tetra', '2tetra', 'alg1b', 'flat_2neigbs']

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "OK. Time: %8.3f"%(time.time()-t0)
