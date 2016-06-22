from __future__ import division
import time

import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
from cctbx import geometry_restraints
import hydrogen_connectivity
import hydrogen_parametrization

#----------------------------------------------------
# This test checks the parameterization of hydrogen atoms
# for amino acids
# Steps:
# 1) determine parameterization
# 2) Compare calculated position of H from parameterization
# to input position
# test fails if distance is > 0.001 A (=precision of coordinates)
#----------------------------------------------------

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
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints,
    normalization = False)

  sites_cart = xray_structure.sites_cart()
  names = list(pdb_hierarchy.atoms().extract_name())
  atoms_list = list(pdb_hierarchy.atoms_with_labels())

  bond_proxies_simple, asu = restraints_manager.geometry.get_all_bond_proxies(
      sites_cart = sites_cart)
  angle_proxies = restraints_manager.geometry.get_all_angle_proxies()

  connectivity = hydrogen_connectivity.determine_H_neighbors(
    bond_proxies   = bond_proxies_simple,
    angle_proxies  = angle_proxies,
    xray_structure = xray_structure)

#--------------------------------------
# This is useful to keep for debugging: human readable output of connectivity
#--------------------------------------
#  for ih in connectivity.keys():
#    if(len(connectivity[ih])==3):
#      string = (" ".join([names[p.iseq] for p in connectivity[ih][2]]))
#    else:
#      string = 'n/a'
#    print  names[ih],': ', names[(connectivity[ih][0]).iseq], \
#      ',', (" ".join([names[p.iseq] for p in connectivity[ih][1]])), ',', string
#--------------------------------------

  h_parameterization = hydrogen_parametrization.get_h_parameterization(
    connectivity   = connectivity,
    sites_cart     = sites_cart,
    names          = names,
    atoms_list     = atoms_list)

# There are 152 H atoms in the pdb_string, check if all of them are recognized
  assert (len(h_parameterization.keys()) == 152), 'Not all H atoms are parameterized'

# For each H atom, check if distance compared to input model is not changed
  n_unk = 0
  for ih in h_parameterization.keys():
    residue = atoms_list[ih].resseq
    hp = h_parameterization[ih]
    h_obj = hydrogen_parametrization.generate_H_positions(
      sites_cart        = sites_cart,
      ih                = ih,
      para_info         = hp)
    assert (h_obj.distance < 0.001), 'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (hp.htype, names[ih], ih, residue)
    if(hp.htype == 'unk'):
      n_unk = n_unk + 1

  assert(n_unk == 0), 'Some H atoms are not recognized'

#--------------------------------------
# This is useful to keep for debugging
#--------------------------------------
    #if(h_obj.distance is not None):
    #  print hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', \
    #    residue, 'distance:', h_obj.distance
    #else:
    #  print hp.htype, 'atom:', names[ih]+' ('+str(ih)+ ') residue:', residue
#--------------------------------------

# Ideal amino acids
pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
SCALE1      0.033333  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033333        0.00000
ATOM      1  N   ARG A   2       3.250  19.540  25.752  1.00  0.00           N
ATOM      2  CA  ARG A   2       3.908  20.348  24.733  1.00  0.00           C
ATOM      3  C   ARG A   2       5.395  20.503  25.033  1.00  0.00           C
ATOM      4  O   ARG A   2       5.777  21.119  26.029  1.00  0.00           O
ATOM      5  CB  ARG A   2       3.246  21.724  24.629  1.00  0.00           C
ATOM      6  CG  ARG A   2       3.853  22.627  23.568  1.00  0.00           C
ATOM      7  CD  ARG A   2       3.152  23.975  23.522  1.00  0.00           C
ATOM      8  NE  ARG A   2       3.720  24.853  22.503  1.00  0.00           N
ATOM      9  CZ  ARG A   2       3.291  26.086  22.257  1.00  0.00           C
ATOM     10  NH1 ARG A   2       2.285  26.593  22.956  1.00  0.00           N1+
ATOM     11  NH2 ARG A   2       3.868  26.814  21.310  1.00  0.00           N
ATOM      0  HA  ARG A   2       3.814  19.890  23.883  1.00  0.00           H
ATOM      0  HB2 ARG A   2       2.303  21.605  24.437  1.00  0.00           H
ATOM      0  HB3 ARG A   2       3.306  22.167  25.490  1.00  0.00           H
ATOM      0  HG2 ARG A   2       4.796  22.757  23.752  1.00  0.00           H
ATOM      0  HG3 ARG A   2       3.790  22.198  22.700  1.00  0.00           H
ATOM      0  HD2 ARG A   2       2.208  23.842  23.343  1.00  0.00           H
ATOM      0  HD3 ARG A   2       3.219  24.403  24.390  1.00  0.00           H
ATOM      0  HE  ARG A   2       4.374  24.553  22.032  1.00  0.00           H
ATOM      0 HH11 ARG A   2       1.909  26.123  23.570  1.00  0.00           H
ATOM      0 HH12 ARG A   2       2.009  27.391  22.795  1.00  0.00           H
ATOM      0 HH21 ARG A   2       4.521  26.488  20.855  1.00  0.00           H
ATOM      0 HH22 ARG A   2       3.590  27.612  21.152  1.00  0.00           H
ATOM     12  N   HIS A   3       6.232  19.940  24.165  1.00  0.00           N
ATOM     13  CA  HIS A   3       7.680  19.999  24.307  1.00  0.00           C
ATOM     14  C   HIS A   3       8.283  20.588  23.041  1.00  0.00           C
ATOM     15  O   HIS A   3       7.916  20.190  21.930  1.00  0.00           O
ATOM     16  CB  HIS A   3       8.266  18.611  24.580  1.00  0.00           C
ATOM     17  CG  HIS A   3       7.752  17.974  25.834  1.00  0.00           C
ATOM     18  ND1 HIS A   3       8.265  18.264  27.080  1.00  0.00           N
ATOM     19  CD2 HIS A   3       6.771  17.062  26.033  1.00  0.00           C
ATOM     20  CE1 HIS A   3       7.622  17.558  27.993  1.00  0.00           C
ATOM     21  NE2 HIS A   3       6.711  16.820  27.384  1.00  0.00           N
ATOM      0  H   HIS A   3       5.969  19.508  23.469  1.00  0.00           H
ATOM      0  HA  HIS A   3       7.897  20.564  25.065  1.00  0.00           H
ATOM      0  HB2 HIS A   3       8.067  18.032  23.828  1.00  0.00           H
ATOM      0  HB3 HIS A   3       9.232  18.682  24.636  1.00  0.00           H
ATOM      0  HD2 HIS A   3       6.238  16.673  25.378  1.00  0.00           H
ATOM      0  HE1 HIS A   3       7.783  17.577  28.909  1.00  0.00           H
ATOM     22  N   LYS A   4       9.206  21.532  23.211  1.00  0.00           N
ATOM     23  CA  LYS A   4       9.863  22.178  22.082  1.00  0.00           C
ATOM     24  C   LYS A   4      11.288  22.530  22.477  1.00  0.00           C
ATOM     25  O   LYS A   4      11.504  23.213  23.482  1.00  0.00           O
ATOM     26  CB  LYS A   4       9.103  23.436  21.645  1.00  0.00           C
ATOM     27  CG  LYS A   4       9.729  24.159  20.463  1.00  0.00           C
ATOM     28  CD  LYS A   4       8.928  25.393  20.084  1.00  0.00           C
ATOM     29  CE  LYS A   4       9.556  26.119  18.906  1.00  0.00           C
ATOM     30  NZ  LYS A   4       8.775  27.325  18.516  1.00  0.00           N1+
ATOM      0  H   LYS A   4       9.467  21.814  23.981  1.00  0.00           H
ATOM      0  HA  LYS A   4       9.872  21.567  21.329  1.00  0.00           H
ATOM      0  HB2 LYS A   4       8.193  23.190  21.416  1.00  0.00           H
ATOM      0  HB3 LYS A   4       9.051  24.048  22.396  1.00  0.00           H
ATOM      0  HG2 LYS A   4      10.638  24.416  20.683  1.00  0.00           H
ATOM      0  HG3 LYS A   4       9.780  23.558  19.703  1.00  0.00           H
ATOM      0  HD2 LYS A   4       8.020  25.136  19.861  1.00  0.00           H
ATOM      0  HD3 LYS A   4       8.874  25.992  20.845  1.00  0.00           H
ATOM      0  HE2 LYS A   4      10.462  26.381  19.134  1.00  0.00           H
ATOM      0  HE3 LYS A   4       9.617  25.515  18.149  1.00  0.00           H
ATOM      0  HZ1 LYS A   4       9.170  27.724  17.826  1.00  0.00           H
ATOM      0  HZ2 LYS A   4       7.950  27.083  18.286  1.00  0.00           H
ATOM      0  HZ3 LYS A   4       8.738  27.890  19.203  1.00  0.00           H
ATOM     31  N   ASP A   5      12.251  22.063  21.685  1.00  0.00           N
ATOM     32  CA  ASP A   5      13.662  22.321  21.940  1.00  0.00           C
ATOM     33  C   ASP A   5      14.368  22.572  20.617  1.00  0.00           C
ATOM     34  O   ASP A   5      14.175  21.821  19.656  1.00  0.00           O
ATOM     35  CB  ASP A   5      14.318  21.148  22.678  1.00  0.00           C
ATOM     36  CG  ASP A   5      15.782  21.396  22.982  1.00  0.00           C
ATOM     37  OD1 ASP A   5      16.154  22.567  23.210  1.00  0.00           O
ATOM     38  OD2 ASP A   5      16.562  20.421  22.993  1.00  0.00           O1-
ATOM      0  H   ASP A   5      12.102  21.587  20.984  1.00  0.00           H
ATOM      0  HA  ASP A   5      13.739  23.104  22.508  1.00  0.00           H
ATOM      0  HB2 ASP A   5      13.842  20.985  23.507  1.00  0.00           H
ATOM      0  HB3 ASP A   5      14.234  20.345  22.141  1.00  0.00           H
ATOM     39  N   GLU A   6      15.182  23.622  20.573  1.00  0.00           N
ATOM     40  CA  GLU A   6      15.919  23.969  19.364  1.00  0.00           C
ATOM     41  C   GLU A   6      17.320  24.468  19.701  1.00  0.00           C
ATOM     42  O   GLU A   6      17.684  24.583  20.871  1.00  0.00           O
ATOM     43  CB  GLU A   6      15.162  25.028  18.560  1.00  0.00           C
ATOM     44  CG  GLU A   6      15.859  25.448  17.276  1.00  0.00           C
ATOM     45  CD  GLU A   6      15.084  26.502  16.509  1.00  0.00           C
ATOM     46  OE1 GLU A   6      13.994  26.896  16.975  1.00  0.00           O
ATOM     47  OE2 GLU A   6      15.564  26.937  15.442  1.00  0.00           O1-
ATOM      0  H   GLU A   6      15.322  24.149  21.238  1.00  0.00           H
ATOM      0  HA  GLU A   6      16.003  23.167  18.825  1.00  0.00           H
ATOM      0  HB2 GLU A   6      14.281  24.686  18.341  1.00  0.00           H
ATOM      0  HB3 GLU A   6      15.032  25.811  19.117  1.00  0.00           H
ATOM      0  HG2 GLU A   6      16.741  25.791  17.488  1.00  0.00           H
ATOM      0  HG3 GLU A   6      15.985  24.670  16.711  1.00  0.00           H
TER
ATOM     48  N   SER B   2      10.369   5.142  24.245  1.00  0.00           N
ATOM     49  CA  SER B   2      11.201   6.333  24.366  1.00  0.00           C
ATOM     50  C   SER B   2      12.675   5.962  24.488  1.00  0.00           C
ATOM     51  O   SER B   2      13.050   5.137  25.320  1.00  0.00           O
ATOM     52  CB  SER B   2      10.767   7.171  25.571  1.00  0.00           C
ATOM     53  OG  SER B   2      10.895   6.435  26.775  1.00  0.00           O
ATOM      0  HA  SER B   2      11.085   6.860  23.560  1.00  0.00           H
ATOM      0  HB2 SER B   2      11.307   7.975  25.621  1.00  0.00           H
ATOM      0  HB3 SER B   2       9.846   7.454  25.457  1.00  0.00           H
ATOM      0  HG  SER B   2      10.808   5.615  26.613  1.00  0.00           H
ATOM     54  N   THR B   3      13.508   6.579  23.652  1.00  0.00           N
ATOM     55  CA  THR B   3      14.944   6.326  23.652  1.00  0.00           C
ATOM     56  C   THR B   3      15.671   7.653  23.507  1.00  0.00           C
ATOM     57  O   THR B   3      15.545   8.321  22.477  1.00  0.00           O
ATOM     58  CB  THR B   3      15.343   5.371  22.522  1.00  0.00           C
ATOM     59  OG1 THR B   3      14.664   4.120  22.689  1.00  0.00           O
ATOM     60  CG2 THR B   3      16.846   5.132  22.533  1.00  0.00           C
ATOM      0  H   THR B   3      13.254   7.157  23.068  1.00  0.00           H
ATOM      0  HA  THR B   3      15.191   5.901  24.488  1.00  0.00           H
ATOM      0  HB  THR B   3      15.094   5.771  21.674  1.00  0.00           H
ATOM      0  HG1 THR B   3      14.882   3.598  22.068  1.00  0.00           H
ATOM      0 HG21 THR B   3      17.083   4.526  21.813  1.00  0.00           H
ATOM      0 HG22 THR B   3      17.309   5.976  22.411  1.00  0.00           H
ATOM      0 HG23 THR B   3      17.106   4.741  23.382  1.00  0.00           H
ATOM     61  N   ASN B   4      16.427   8.031  24.535  1.00  0.00           N
ATOM     62  CA  ASN B   4      17.183   9.276  24.541  1.00  0.00           C
ATOM     63  C   ASN B   4      18.633   8.977  24.886  1.00  0.00           C
ATOM     64  O   ASN B   4      18.908   8.244  25.841  1.00  0.00           O
ATOM     65  CB  ASN B   4      16.594  10.278  25.541  1.00  0.00           C
ATOM     66  CG  ASN B   4      15.152  10.631  25.231  1.00  0.00           C
ATOM     67  OD1 ASN B   4      14.224  10.007  25.745  1.00  0.00           O
ATOM     68  ND2 ASN B   4      14.957  11.638  24.387  1.00  0.00           N
ATOM      0  H   ASN B   4      16.515   7.567  25.254  1.00  0.00           H
ATOM      0  HA  ASN B   4      17.132   9.677  23.659  1.00  0.00           H
ATOM      0  HB2 ASN B   4      16.648   9.907  26.436  1.00  0.00           H
ATOM      0  HB3 ASN B   4      17.130  11.087  25.537  1.00  0.00           H
ATOM      0 HD21 ASN B   4      14.157  11.877  24.180  1.00  0.00           H
ATOM      0 HD22 ASN B   4      15.631  12.051  24.047  1.00  0.00           H
ATOM     69  N   GLN B   5      19.552   9.543  24.112  1.00  0.00           N
ATOM     70  CA  GLN B   5      20.978   9.337  24.339  1.00  0.00           C
ATOM     71  C   GLN B   5      21.781  10.569  23.935  1.00  0.00           C
ATOM     72  O   GLN B   5      21.698  11.030  22.797  1.00  0.00           O
ATOM     73  CB  GLN B   5      21.472   8.111  23.569  1.00  0.00           C
ATOM     74  CG  GLN B   5      22.953   7.818  23.751  1.00  0.00           C
ATOM     75  CD  GLN B   5      23.407   6.598  22.974  1.00  0.00           C
ATOM     76  OE1 GLN B   5      24.582   6.232  23.001  1.00  0.00           O
ATOM     77  NE2 GLN B   5      22.474   5.961  22.275  1.00  0.00           N
ATOM      0  H   GLN B   5      19.369  10.054  23.445  1.00  0.00           H
ATOM      0  HA  GLN B   5      21.110   9.185  25.288  1.00  0.00           H
ATOM      0  HB2 GLN B   5      20.962   7.336  23.852  1.00  0.00           H
ATOM      0  HB3 GLN B   5      21.291   8.240  22.625  1.00  0.00           H
ATOM      0  HG2 GLN B   5      23.469   8.589  23.466  1.00  0.00           H
ATOM      0  HG3 GLN B   5      23.139   7.684  24.694  1.00  0.00           H
ATOM      0 HE21 GLN B   5      21.663   6.247  22.280  1.00  0.00           H
ATOM      0 HE22 GLN B   5      22.681   5.263  21.817  1.00  0.00           H
TER
ATOM     78  N   CYS C   2      14.570  22.019   7.108  1.00  0.00           N
ATOM     79  CA  CYS C   2      15.123  23.368   7.079  1.00  0.00           C
ATOM     80  C   CYS C   2      16.641  23.340   7.225  1.00  0.00           C
ATOM     81  O   CYS C   2      17.175  22.720   8.144  1.00  0.00           O
ATOM     82  CB  CYS C   2      14.504  24.226   8.184  1.00  0.00           C
ATOM     83  SG  CYS C   2      15.103  25.931   8.231  1.00  0.00           S
ATOM      0  HA  CYS C   2      14.906  23.761   6.219  1.00  0.00           H
ATOM      0  HB2 CYS C   2      13.541  24.236   8.069  1.00  0.00           H
ATOM      0  HB3 CYS C   2      14.683  23.808   9.041  1.00  0.00           H
ATOM      0  HG  CYS C   2      14.561  26.521   9.124  1.00  0.00           H
ATOM     84  N   GLY C   3      17.331  24.017   6.311  1.00  0.00           N
ATOM     85  CA  GLY C   3      18.772  24.074   6.332  1.00  0.00           C
ATOM     86  C   GLY C   3      19.290  25.483   6.537  1.00  0.00           C
ATOM     87  O   GLY C   3      18.548  26.464   6.432  1.00  0.00           O
ATOM      0  H   GLY C   3      16.969  24.454   5.664  1.00  0.00           H
ATOM      0  HA2 GLY C   3      19.106  23.503   7.041  1.00  0.00           H
ATOM      0  HA3 GLY C   3      19.120  23.722   5.497  1.00  0.00           H
ATOM     88  N   PRO C   4      20.591  25.610   6.837  1.00  0.00           N
ATOM     89  CA  PRO C   4      21.227  26.912   7.061  1.00  0.00           C
ATOM     90  C   PRO C   4      21.373  27.723   5.777  1.00  0.00           C
ATOM     91  O   PRO C   4      21.372  27.134   4.696  1.00  0.00           O
ATOM     92  CB  PRO C   4      22.598  26.536   7.627  1.00  0.00           C
ATOM     93  CG  PRO C   4      22.876  25.192   7.052  1.00  0.00           C
ATOM     94  CD  PRO C   4      21.544  24.497   6.989  1.00  0.00           C
ATOM      0  HA  PRO C   4      20.702  27.478   7.648  1.00  0.00           H
ATOM      0  HB2 PRO C   4      23.276  27.179   7.368  1.00  0.00           H
ATOM      0  HB3 PRO C   4      22.586  26.511   8.597  1.00  0.00           H
ATOM      0  HG2 PRO C   4      23.274  25.265   6.170  1.00  0.00           H
ATOM      0  HG3 PRO C   4      23.501  24.698   7.605  1.00  0.00           H
ATOM      0  HD2 PRO C   4      21.497  23.880   6.242  1.00  0.00           H
ATOM      0  HD3 PRO C   4      21.370  23.983   7.793  1.00  0.00           H
TER
ATOM     95  N   ALA D   2       2.568   5.353   8.137  1.00  0.00           N
ATOM     96  CA  ALA D   2       3.387   6.519   8.442  1.00  0.00           C
ATOM     97  C   ALA D   2       4.855   6.132   8.587  1.00  0.00           C
ATOM     98  O   ALA D   2       5.215   5.338   9.456  1.00  0.00           O
ATOM     99  CB  ALA D   2       2.890   7.200   9.709  1.00  0.00           C
ATOM      0  HA  ALA D   2       3.311   7.143   7.703  1.00  0.00           H
ATOM      0  HB1 ALA D   2       3.444   7.973   9.899  1.00  0.00           H
ATOM      0  HB2 ALA D   2       1.971   7.484   9.585  1.00  0.00           H
ATOM      0  HB3 ALA D   2       2.938   6.578  10.451  1.00  0.00           H
ATOM    100  N   LEU D   3       5.699   6.699   7.729  1.00  0.00           N
ATOM    101  CA  LEU D   3       7.136   6.436   7.731  1.00  0.00           C
ATOM    102  C   LEU D   3       7.858   7.764   7.938  1.00  0.00           C
ATOM    103  O   LEU D   3       8.085   8.514   6.984  1.00  0.00           O
ATOM    104  CB  LEU D   3       7.568   5.760   6.433  1.00  0.00           C
ATOM    105  CG  LEU D   3       6.916   4.412   6.120  1.00  0.00           C
ATOM    106  CD1 LEU D   3       7.444   3.850   4.809  1.00  0.00           C
ATOM    107  CD2 LEU D   3       7.140   3.429   7.258  1.00  0.00           C
ATOM      0  H   LEU D   3       5.450   7.254   7.121  1.00  0.00           H
ATOM      0  HA  LEU D   3       7.364   5.826   8.450  1.00  0.00           H
ATOM      0  HB2 LEU D   3       7.383   6.366   5.698  1.00  0.00           H
ATOM      0  HB3 LEU D   3       8.529   5.634   6.461  1.00  0.00           H
ATOM      0  HG  LEU D   3       5.961   4.552   6.025  1.00  0.00           H
ATOM      0 HD11 LEU D   3       7.020   2.997   4.628  1.00  0.00           H
ATOM      0 HD12 LEU D   3       7.246   4.468   4.088  1.00  0.00           H
ATOM      0 HD13 LEU D   3       8.404   3.726   4.873  1.00  0.00           H
ATOM      0 HD21 LEU D   3       6.720   2.582   7.041  1.00  0.00           H
ATOM      0 HD22 LEU D   3       8.092   3.294   7.387  1.00  0.00           H
ATOM      0 HD23 LEU D   3       6.751   3.783   8.073  1.00  0.00           H
ATOM    108  N   ILE D   4       8.217   8.049   9.185  1.00  0.00           N
ATOM    109  CA  ILE D   4       8.916   9.275   9.552  1.00  0.00           C
ATOM    110  C   ILE D   4      10.371   8.898   9.802  1.00  0.00           C
ATOM    111  O   ILE D   4      10.733   8.457  10.899  1.00  0.00           O
ATOM    112  CB  ILE D   4       8.289   9.950  10.777  1.00  0.00           C
ATOM    113  CG1 ILE D   4       6.807  10.237  10.527  1.00  0.00           C
ATOM    114  CG2 ILE D   4       9.030  11.236  11.115  1.00  0.00           C
ATOM    115  CD1 ILE D   4       6.102  10.879  11.702  1.00  0.00           C
ATOM      0  H   ILE D   4       8.059   7.528   9.851  1.00  0.00           H
ATOM      0  HA  ILE D   4       8.849   9.926   8.836  1.00  0.00           H
ATOM      0  HB  ILE D   4       8.363   9.346  11.532  1.00  0.00           H
ATOM      0 HG12 ILE D   4       6.725  10.818   9.755  1.00  0.00           H
ATOM      0 HG13 ILE D   4       6.358   9.406  10.306  1.00  0.00           H
ATOM      0 HG21 ILE D   4       8.622  11.651  11.891  1.00  0.00           H
ATOM      0 HG22 ILE D   4       9.959  11.033  11.308  1.00  0.00           H
ATOM      0 HG23 ILE D   4       8.982  11.845  10.361  1.00  0.00           H
ATOM      0 HD11 ILE D   4       5.171  11.033  11.477  1.00  0.00           H
ATOM      0 HD12 ILE D   4       6.155  10.291  12.472  1.00  0.00           H
ATOM      0 HD13 ILE D   4       6.527  11.725  11.912  1.00  0.00           H
ATOM    116  N   MET D   5      11.211   9.068   8.783  1.00  0.00           N
ATOM    117  CA  MET D   5      12.632   8.753   8.871  1.00  0.00           C
ATOM    118  C   MET D   5      13.427   9.929   8.326  1.00  0.00           C
ATOM    119  O   MET D   5      13.353  10.230   7.130  1.00  0.00           O
ATOM    120  CB  MET D   5      12.966   7.474   8.097  1.00  0.00           C
ATOM    121  CG  MET D   5      12.265   6.229   8.618  1.00  0.00           C
ATOM    122  SD  MET D   5      12.713   4.743   7.700  1.00  0.00           S
ATOM    123  CE  MET D   5      11.715   3.513   8.537  1.00  0.00           C
ATOM      0  H   MET D   5      10.969   9.372   8.016  1.00  0.00           H
ATOM      0  HA  MET D   5      12.867   8.598   9.799  1.00  0.00           H
ATOM      0  HB2 MET D   5      12.728   7.601   7.165  1.00  0.00           H
ATOM      0  HB3 MET D   5      13.925   7.330   8.128  1.00  0.00           H
ATOM      0  HG2 MET D   5      12.487   6.105   9.554  1.00  0.00           H
ATOM      0  HG3 MET D   5      11.305   6.358   8.568  1.00  0.00           H
ATOM      0  HE1 MET D   5      11.864   2.644   8.132  1.00  0.00           H
ATOM      0  HE2 MET D   5      11.962   3.479   9.474  1.00  0.00           H
ATOM      0  HE3 MET D   5      10.777   3.749   8.459  1.00  0.00           H
ATOM    124  N   PHE D   6      14.184  10.590   9.200  1.00  0.00           N
ATOM    125  CA  PHE D   6      15.007  11.731   8.816  1.00  0.00           C
ATOM    126  C   PHE D   6      16.370  11.591   9.473  1.00  0.00           C
ATOM    127  O   PHE D   6      16.463  11.528  10.703  1.00  0.00           O
ATOM    128  CB  PHE D   6      14.348  13.054   9.219  1.00  0.00           C
ATOM    129  CG  PHE D   6      13.042  13.318   8.524  1.00  0.00           C
ATOM    130  CD1 PHE D   6      11.844  12.926   9.099  1.00  0.00           C
ATOM    131  CD2 PHE D   6      13.013  13.959   7.297  1.00  0.00           C
ATOM    132  CE1 PHE D   6      10.642  13.168   8.461  1.00  0.00           C
ATOM    133  CE2 PHE D   6      11.814  14.204   6.655  1.00  0.00           C
ATOM    134  CZ  PHE D   6      10.627  13.808   7.238  1.00  0.00           C
ATOM      0  H   PHE D   6      14.234  10.387  10.034  1.00  0.00           H
ATOM      0  HA  PHE D   6      15.105  11.742   7.851  1.00  0.00           H
ATOM      0  HB2 PHE D   6      14.200  13.053  10.178  1.00  0.00           H
ATOM      0  HB3 PHE D   6      14.960  13.782   9.027  1.00  0.00           H
ATOM      0  HD1 PHE D   6      11.849  12.495   9.923  1.00  0.00           H
ATOM      0  HD2 PHE D   6      13.810  14.228   6.900  1.00  0.00           H
ATOM      0  HE1 PHE D   6       9.844  12.900   8.855  1.00  0.00           H
ATOM      0  HE2 PHE D   6      11.807  14.635   5.831  1.00  0.00           H
ATOM      0  HZ  PHE D   6       9.819  13.972   6.808  1.00  0.00           H
ATOM    135  N   TRP D   7      17.421  11.543   8.658  1.00  0.00           N
ATOM    136  CA  TRP D   7      18.790  11.407   9.142  1.00  0.00           C
ATOM    137  C   TRP D   7      19.627  12.537   8.562  1.00  0.00           C
ATOM    138  O   TRP D   7      19.653  12.729   7.342  1.00  0.00           O
ATOM    139  CB  TRP D   7      19.381  10.047   8.758  1.00  0.00           C
ATOM    140  CG  TRP D   7      18.622   8.883   9.319  1.00  0.00           C
ATOM    141  CD1 TRP D   7      17.587   8.219   8.729  1.00  0.00           C
ATOM    142  CD2 TRP D   7      18.839   8.246  10.584  1.00  0.00           C
ATOM    143  NE1 TRP D   7      17.146   7.207   9.547  1.00  0.00           N
ATOM    144  CE2 TRP D   7      17.898   7.202  10.692  1.00  0.00           C
ATOM    145  CE3 TRP D   7      19.738   8.456  11.634  1.00  0.00           C
ATOM    146  CZ2 TRP D   7      17.831   6.371  11.809  1.00  0.00           C
ATOM    147  CZ3 TRP D   7      19.669   7.630  12.742  1.00  0.00           C
ATOM    148  CH2 TRP D   7      18.722   6.601  12.820  1.00  0.00           C
ATOM      0  H   TRP D   7      17.358  11.589   7.802  1.00  0.00           H
ATOM      0  HA  TRP D   7      18.793  11.459  10.111  1.00  0.00           H
ATOM      0  HB2 TRP D   7      19.403   9.974   7.791  1.00  0.00           H
ATOM      0  HB3 TRP D   7      20.300  10.003   9.066  1.00  0.00           H
ATOM      0  HD1 TRP D   7      17.231   8.422   7.894  1.00  0.00           H
ATOM      0  HE1 TRP D   7      16.503   6.665   9.369  1.00  0.00           H
ATOM      0  HE3 TRP D   7      20.370   9.137  11.589  1.00  0.00           H
ATOM      0  HZ2 TRP D   7      17.204   5.686  11.864  1.00  0.00           H
ATOM      0  HZ3 TRP D   7      20.262   7.760  13.447  1.00  0.00           H
ATOM      0  HH2 TRP D   7      18.699   6.062  13.577  1.00  0.00           H
ATOM    149  N   TYR D   8      20.307  13.279   9.434  1.00  0.00           N
ATOM    150  CA  TYR D   8      21.159  14.394   9.030  1.00  0.00           C
ATOM    151  C   TYR D   8      22.486  14.274   9.766  1.00  0.00           C
ATOM    152  O   TYR D   8      22.530  14.383  10.996  1.00  0.00           O
ATOM    153  CB  TYR D   8      20.489  15.737   9.326  1.00  0.00           C
ATOM    154  CG  TYR D   8      19.184  15.946   8.592  1.00  0.00           C
ATOM    155  CD1 TYR D   8      19.164  16.475   7.308  1.00  0.00           C
ATOM    156  CD2 TYR D   8      17.972  15.615   9.183  1.00  0.00           C
ATOM    157  CE1 TYR D   8      17.974  16.668   6.632  1.00  0.00           C
ATOM    158  CE2 TYR D   8      16.776  15.805   8.515  1.00  0.00           C
ATOM    159  CZ  TYR D   8      16.783  16.331   7.241  1.00  0.00           C
ATOM    160  OH  TYR D   8      15.596  16.522   6.572  1.00  0.00           O
ATOM      0  H   TYR D   8      20.286  13.148  10.284  1.00  0.00           H
ATOM      0  HA  TYR D   8      21.309  14.358   8.072  1.00  0.00           H
ATOM      0  HB2 TYR D   8      20.328  15.804  10.280  1.00  0.00           H
ATOM      0  HB3 TYR D   8      21.100  16.452   9.090  1.00  0.00           H
ATOM      0  HD1 TYR D   8      19.966  16.703   6.896  1.00  0.00           H
ATOM      0  HD2 TYR D   8      17.964  15.260  10.042  1.00  0.00           H
ATOM      0  HE1 TYR D   8      17.976  17.023   5.772  1.00  0.00           H
ATOM      0  HE2 TYR D   8      15.971  15.579   8.923  1.00  0.00           H
ATOM      0  HH  TYR D   8      14.957  16.274   7.058  1.00  0.00           H
ATOM    161  N   VAL D   9      23.561  14.050   9.017  1.00  0.00           N
ATOM    162  CA  VAL D   9      24.888  13.916   9.604  1.00  0.00           C
ATOM    163  C   VAL D   9      25.718  15.163   9.324  1.00  0.00           C
ATOM    164  O   VAL D   9      25.820  15.609   8.181  1.00  0.00           O
ATOM    165  CB  VAL D   9      25.604  12.659   9.082  1.00  0.00           C
ATOM    166  CG1 VAL D   9      26.992  12.547   9.693  1.00  0.00           C
ATOM    167  CG2 VAL D   9      24.780  11.416   9.383  1.00  0.00           C
ATOM      0  H   VAL D   9      23.542  13.972   8.161  1.00  0.00           H
ATOM      0  HA  VAL D   9      24.785  13.820  10.564  1.00  0.00           H
ATOM      0  HB  VAL D   9      25.701  12.734   8.120  1.00  0.00           H
ATOM      0 HG11 VAL D   9      27.431  11.751   9.355  1.00  0.00           H
ATOM      0 HG12 VAL D   9      27.514  13.329   9.456  1.00  0.00           H
ATOM      0 HG13 VAL D   9      26.917  12.490  10.658  1.00  0.00           H
ATOM      0 HG21 VAL D   9      25.243  10.632   9.049  1.00  0.00           H
ATOM      0 HG22 VAL D   9      24.656  11.334  10.342  1.00  0.00           H
ATOM      0 HG23 VAL D   9      23.915  11.489   8.951  1.00  0.00           H
TER
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "OK. Time: %8.3f"%(time.time()-t0)

