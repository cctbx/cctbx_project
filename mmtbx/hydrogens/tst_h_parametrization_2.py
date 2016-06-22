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
# for nucleic acids
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

#-------------------------------------------------------------
# This is useful to keep for debugging: human readable output of connectivity
#-------------------------------------------------------------
#  for ih in connectivity.keys():
#    if(len(connectivity[ih])==3):
#      string = (" ".join([names[p.iseq] for p in connectivity[ih][2]]))
#    else:
#      string = 'n/a'
#    print  names[ih],': ', names[(connectivity[ih][0]).iseq], \
#      ',', (" ".join([names[p.iseq] for p in connectivity[ih][1]])), ',', string
#-------------------------------------------------------------

  h_parameterization = hydrogen_parametrization.get_h_parameterization(
    connectivity   = connectivity,
    sites_cart     = sites_cart,
    names          = names,
    atoms_list     = atoms_list)

# There are 90 H atoms in the pdb_string, check if all of them are recognized
  assert (len(h_parameterization.keys()) == 90), 'Not all H atoms are parameterized'

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


# DNA and RNA nucleic acids
pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
SCALE1      0.033333  0.000000  0.000000        0.00000
SCALE2      0.000000  0.033333  0.000000        0.00000
SCALE3      0.000000  0.000000  0.033333        0.00000
ATOM      1  P    DA A   1       4.636   5.591  10.123  1.00 76.88           P
ATOM      2  OP1  DA A   1       4.480   4.227  10.674  1.00 78.59           O
ATOM      3  OP2  DA A   1       3.818   6.703  10.656  1.00 75.20           O
ATOM      4  O5'  DA A   1       6.175   6.019  10.221  1.00 78.76           O
ATOM      5  C5'  DA A   1       7.164   5.047  10.540  1.00 77.93           C
ATOM      6  C4'  DA A   1       8.559   5.611  10.340  1.00 78.66           C
ATOM      7  O4'  DA A   1       8.765   5.904   8.934  1.00 79.49           O
ATOM      8  C3'  DA A   1       8.832   6.920  11.064  1.00 78.23           C
ATOM      9  O3'  DA A   1       9.313   6.657  12.377  1.00 77.17           O
ATOM     10  C2'  DA A   1       9.909   7.550  10.190  1.00 79.34           C
ATOM     11  C1'  DA A   1       9.458   7.132   8.792  1.00 78.22           C
ATOM     12  N9   DA A   1       8.566   8.096   8.153  1.00 20.00           N
ATOM     13  C8   DA A   1       7.201   8.043   8.094  1.00 20.00           C
ATOM     14  N7   DA A   1       6.660   9.051   7.452  1.00 20.00           N
ATOM     15  C5   DA A   1       7.745   9.818   7.063  1.00 20.00           C
ATOM     16  C6   DA A   1       7.844  11.026   6.343  1.00 20.00           C
ATOM     17  N6   DA A   1       6.784  11.692   5.873  1.00 20.00           N
ATOM     18  N1   DA A   1       9.079  11.524   6.125  1.00 20.00           N
ATOM     19  C2   DA A   1      10.137  10.855   6.598  1.00 20.00           C
ATOM     20  N3   DA A   1      10.169   9.714   7.286  1.00 20.00           N
ATOM     21  C4   DA A   1       8.928   9.244   7.487  1.00 20.00           C
ATOM      0  H5'  DA A   1       7.045   4.263   9.982  1.00 77.93           H
ATOM      0 H5''  DA A   1       7.056   4.759  11.460  1.00 77.93           H
ATOM      0  H4'  DA A   1       9.149   4.929  10.697  1.00 78.66           H
ATOM      0  H3'  DA A   1       8.054   7.488  11.178  1.00 78.23           H
ATOM      0  H2'  DA A   1       9.941   8.514  10.289  1.00 79.34           H
ATOM      0 H2''  DA A   1      10.794   7.214  10.403  1.00 79.34           H
ATOM      0  H1'  DA A   1      10.244   7.066   8.227  1.00 78.22           H
ATOM      0  H8   DA A   1       6.702   7.356   8.474  1.00 20.00           H
ATOM      0  H61  DA A   1       6.895  12.426   5.439  1.00 20.00           H
ATOM      0  H62  DA A   1       5.991  11.386   6.005  1.00 20.00           H
ATOM      0  H2   DA A   1      10.965  11.240   6.422  1.00 20.00           H
ATOM     22  P    DC A   2       8.850   7.580  13.608  1.00 76.88           P
ATOM     23  OP1  DC A   2       9.291   6.918  14.856  1.00 78.59           O
ATOM     24  OP2  DC A   2       7.419   7.903  13.413  1.00 75.20           O
ATOM     25  O5'  DC A   2       9.687   8.930  13.415  1.00 78.76           O
ATOM     26  C5'  DC A   2      11.109   8.896  13.450  1.00 77.93           C
ATOM     27  C4'  DC A   2      11.692  10.226  13.009  1.00 78.66           C
ATOM     28  O4'  DC A   2      11.337  10.474  11.627  1.00 79.49           O
ATOM     29  C3'  DC A   2      11.176  11.435  13.770  1.00 78.23           C
ATOM     30  O3'  DC A   2      11.953  11.639  14.943  1.00 77.17           O
ATOM     31  C2'  DC A   2      11.380  12.559  12.760  1.00 79.34           C
ATOM     32  C1'  DC A   2      11.087  11.856  11.433  1.00 78.22           C
ATOM     33  N1   DC A   2       9.679  12.026  10.961  1.00 20.00           N
ATOM     34  C2   DC A   2       9.303  13.211  10.318  1.00 20.00           C
ATOM     35  O2   DC A   2      10.145  14.102  10.153  1.00 20.00           O
ATOM     36  N3   DC A   2       8.021  13.349   9.895  1.00 20.00           N
ATOM     37  C4   DC A   2       7.143  12.364  10.092  1.00 20.00           C
ATOM     38  N4   DC A   2       5.892  12.546   9.657  1.00 20.00           N
ATOM     39  C5   DC A   2       7.508  11.150  10.743  1.00 20.00           C
ATOM     40  C6   DC A   2       8.774  11.024  11.156  1.00 20.00           C
ATOM      0  H5'  DC A   2      11.434   8.188  12.872  1.00 77.93           H
ATOM      0 H5''  DC A   2      11.410   8.689  14.349  1.00 77.93           H
ATOM      0  H4'  DC A   2      12.645  10.137  13.169  1.00 78.66           H
ATOM      0  H3'  DC A   2      10.258  11.359  14.074  1.00 78.23           H
ATOM      0  H2'  DC A   2      10.777  13.302  12.918  1.00 79.34           H
ATOM      0 H2''  DC A   2      12.282  12.915  12.791  1.00 79.34           H
ATOM      0  H1'  DC A   2      11.654  12.253  10.753  1.00 78.22           H
ATOM      0  H41  DC A   2       5.304  11.929   9.770  1.00 20.00           H
ATOM      0  H42  DC A   2       5.676  13.280   9.265  1.00 20.00           H
ATOM      0  H5   DC A   2       6.888  10.470  10.876  1.00 20.00           H
ATOM      0  H6   DC A   2       9.041  10.242  11.582  1.00 20.00           H
ATOM     41  P    DG A   3      11.256  12.163  16.293  1.00 76.88           P
ATOM     42  OP1  DG A   3      12.229  11.992  17.395  1.00 78.59           O
ATOM     43  OP2  DG A   3       9.922  11.529  16.384  1.00 75.20           O
ATOM     44  O5'  DG A   3      11.040  13.726  16.034  1.00 78.76           O
ATOM     45  C5'  DG A   3      12.162  14.572  15.808  1.00 77.93           C
ATOM     46  C4'  DG A   3      11.716  15.951  15.356  1.00 78.66           C
ATOM     47  O4'  DG A   3      11.063  15.849  14.066  1.00 79.49           O
ATOM     48  C3'  DG A   3      10.703  16.628  16.265  1.00 78.23           C
ATOM     49  O3'  DG A   3      11.375  17.365  17.278  1.00 77.17           O
ATOM     50  C2'  DG A   3       9.970  17.548  15.297  1.00 79.34           C
ATOM     51  C1'  DG A   3       9.929  16.697  14.029  1.00 78.22           C
ATOM     52  N9   DG A   3       8.732  15.867  13.922  1.00 20.00           N
ATOM     53  C8   DG A   3       8.613  14.540  14.258  1.00 20.00           C
ATOM     54  N7   DG A   3       7.419  14.058  14.055  1.00 20.00           N
ATOM     55  C5   DG A   3       6.699  15.133  13.551  1.00 20.00           C
ATOM     56  C6   DG A   3       5.345  15.214  13.145  1.00 20.00           C
ATOM     57  O6   DG A   3       4.487  14.320  13.152  1.00 20.00           O
ATOM     58  N1   DG A   3       5.018  16.491  12.696  1.00 20.00           N
ATOM     59  C2   DG A   3       5.888  17.554  12.646  1.00 20.00           C
ATOM     60  N2   DG A   3       5.388  18.709  12.183  1.00 20.00           N
ATOM     61  N3   DG A   3       7.158  17.492  13.023  1.00 20.00           N
ATOM     62  C4   DG A   3       7.493  16.254  13.463  1.00 20.00           C
ATOM      0  H5'  DG A   3      12.739  14.178  15.136  1.00 77.93           H
ATOM      0 H5''  DG A   3      12.685  14.647  16.622  1.00 77.93           H
ATOM      0  H4'  DG A   3      12.530  16.478  15.350  1.00 78.66           H
ATOM      0  H3'  DG A   3      10.111  16.019  16.734  1.00 78.23           H
ATOM      0  H2'  DG A   3       9.081  17.774  15.612  1.00 79.34           H
ATOM      0 H2''  DG A   3      10.443  18.383  15.160  1.00 79.34           H
ATOM      0  H1'  DG A   3       9.921  17.293  13.264  1.00 78.22           H
ATOM      0  H8   DG A   3       9.316  14.035  14.597  1.00 20.00           H
ATOM      0  H1   DG A   3       4.211  16.625  12.430  1.00 20.00           H
ATOM      0  H21  DG A   3       5.890  19.405  12.131  1.00 20.00           H
ATOM      0  H22  DG A   3       4.565  18.754  11.938  1.00 20.00           H
ATOM     63  P    DT A   4      10.798  17.392  18.778  1.00 76.88           P
ATOM     64  OP1  DT A   4      11.827  18.009  19.643  1.00 78.59           O
ATOM     65  OP2  DT A   4      10.292  16.037  19.088  1.00 75.20           O
ATOM     66  O5'  DT A   4       9.540  18.377  18.689  1.00 78.76           O
ATOM     67  C5'  DT A   4       9.726  19.734  18.305  1.00 77.93           C
ATOM     68  C4'  DT A   4       8.390  20.427  18.106  1.00 78.66           C
ATOM     69  O4'  DT A   4       7.675  19.792  17.018  1.00 79.49           O
ATOM     70  C3'  DT A   4       7.446  20.351  19.292  1.00 78.23           C
ATOM     71  O3'  DT A   4       7.714  21.415  20.196  1.00 77.17           O
ATOM     72  C2'  DT A   4       6.084  20.512  18.626  1.00 79.34           C
ATOM     73  C1'  DT A   4       6.286  19.777  17.299  1.00 78.22           C
ATOM     74  N1   DT A   4       5.819  18.358  17.323  1.00 20.00           N
ATOM     75  C2   DT A   4       4.475  18.090  17.202  1.00 20.00           C
ATOM     76  O2   DT A   4       3.629  18.958  17.078  1.00 20.00           O
ATOM     77  N3   DT A   4       4.151  16.760  17.234  1.00 20.00           N
ATOM     78  C4   DT A   4       5.019  15.692  17.371  1.00 20.00           C
ATOM     79  O4   DT A   4       4.630  14.528  17.388  1.00 20.00           O
ATOM     80  C5   DT A   4       6.415  16.042  17.491  1.00 20.00           C
ATOM     81  C7   DT A   4       7.451  14.968  17.644  1.00 20.00           C
ATOM     82  C6   DT A   4       6.743  17.342  17.460  1.00 20.00           C
ATOM      0  H5'  DT A   4      10.241  19.774  17.484  1.00 77.93           H
ATOM      0 H5''  DT A   4      10.238  20.199  18.985  1.00 77.93           H
ATOM      0  H4'  DT A   4       8.621  21.356  17.947  1.00 78.66           H
ATOM      0  H3'  DT A   4       7.520  19.538  19.816  1.00 78.23           H
ATOM      0  H2'  DT A   4       5.373  20.117  19.155  1.00 79.34           H
ATOM      0 H2''  DT A   4       5.854  21.445  18.492  1.00 79.34           H
ATOM      0 HO3'  DT A   4       7.187  21.368  20.849  1.00 77.17           H
ATOM      0  H1'  DT A   4       5.758  20.230  16.623  1.00 78.22           H
ATOM      0  H3   DT A   4       3.316  16.568  17.161  1.00 20.00           H
ATOM      0  H71  DT A   4       8.287  15.266  17.252  1.00 20.00           H
ATOM      0  H72  DT A   4       7.586  14.780  18.586  1.00 20.00           H
ATOM      0  H73  DT A   4       7.152  14.163  17.193  1.00 20.00           H
ATOM      0  H6   DT A   4       7.641  17.571  17.535  1.00 20.00           H
TER
ATOM     83  P     A B   1      19.152   8.087  13.325  1.00 76.88           P
ATOM     84  OP1   A B   1      19.178   6.878  14.186  1.00 78.59           O
ATOM     85  OP2   A B   1      18.253   9.216  13.672  1.00 75.20           O
ATOM     86  O5'   A B   1      20.638   8.654  13.226  1.00 78.76           O
ATOM     87  C5'   A B   1      21.757   7.782  13.299  1.00 77.93           C
ATOM     88  C4'   A B   1      23.035   8.494  12.934  1.00 78.66           C
ATOM     89  O4'   A B   1      22.965   8.956  11.560  1.00 79.49           O
ATOM     90  C3'   A B   1      23.349   9.751  13.730  1.00 78.23           C
ATOM     91  O3'   A B   1      23.937   9.475  14.990  1.00 77.17           O
ATOM     92  C2'   A B   1      24.260  10.524  12.786  1.00 79.34           C
ATOM     93  O2'   A B   1      25.579   9.998  12.818  1.00 79.34           O
ATOM     94  C1'   A B   1      23.640  10.189  11.430  1.00 78.22           C
ATOM     95  N9    A B   1      22.667  11.212  10.999  1.00 20.00           N
ATOM     96  C8    A B   1      21.297  11.130  11.018  1.00 20.00           C
ATOM     97  N7    A B   1      20.697  12.206  10.570  1.00 20.00           N
ATOM     98  C5    A B   1      21.743  13.054  10.232  1.00 20.00           C
ATOM     99  C6    A B   1      21.774  14.354   9.698  1.00 20.00           C
ATOM    100  N6    A B   1      20.680  15.057   9.397  1.00 20.00           N
ATOM    101  N1    A B   1      22.984  14.915   9.480  1.00 20.00           N
ATOM    102  C2    A B   1      24.082  14.210   9.781  1.00 20.00           C
ATOM    103  N3    A B   1      24.181  12.983  10.287  1.00 20.00           N
ATOM    104  C4    A B   1      22.963  12.454  10.492  1.00 20.00           C
ATOM      0  H5'   A B   1      21.622   7.030  12.702  1.00 77.93           H
ATOM      0 H5''   A B   1      21.830   7.422  14.197  1.00 77.93           H
ATOM      0  H4'   A B   1      23.717   7.828  13.113  1.00 78.66           H
ATOM      0  H3'   A B   1      22.558  10.255  13.980  1.00 78.23           H
ATOM      0  H2'   A B   1      24.325  11.469  12.997  1.00 79.34           H
ATOM      0 HO2'   A B   1      25.699   9.585  13.540  1.00 79.34           H
ATOM      0  H1'   A B   1      24.346  10.152  10.766  1.00 78.22           H
ATOM      0  H8    A B   1      20.837  10.381  11.321  1.00 20.00           H
ATOM      0  H61   A B   1      20.754  15.850   9.072  1.00 20.00           H
ATOM      0  H62   A B   1      19.901  14.717   9.528  1.00 20.00           H
ATOM      0  H2    A B   1      24.890  14.639   9.613  1.00 20.00           H
ATOM    105  P     C B   2      23.517  10.332  16.283  1.00 76.88           P
ATOM    106  OP1   C B   2      24.086   9.664  17.481  1.00 78.59           O
ATOM    107  OP2   C B   2      22.056  10.588  16.210  1.00 75.20           O
ATOM    108  O5'   C B   2      24.271  11.721  16.088  1.00 78.76           O
ATOM    109  C5'   C B   2      25.689  11.779  16.026  1.00 77.93           C
ATOM    110  C4'   C B   2      26.168  13.144  15.602  1.00 78.66           C
ATOM    111  O4'   C B   2      25.690  13.438  14.264  1.00 79.49           O
ATOM    112  C3'   C B   2      25.672  14.314  16.437  1.00 78.23           C
ATOM    113  O3'   C B   2      26.410  14.492  17.633  1.00 77.17           O
ATOM    114  C2'   C B   2      25.783  15.485  15.469  1.00 79.34           C
ATOM    115  O2'   C B   2      27.125  15.940  15.384  1.00 79.34           O
ATOM    116  C1'   C B   2      25.413  14.818  14.143  1.00 78.22           C
ATOM    117  N1    C B   2      23.980  14.988  13.809  1.00 20.00           N
ATOM    118  C2    C B   2      23.562  16.185  13.220  1.00 20.00           C
ATOM    119  O2    C B   2      24.401  17.070  12.995  1.00 20.00           O
ATOM    120  N3    C B   2      22.255  16.350  12.910  1.00 20.00           N
ATOM    121  C4    C B   2      21.380  15.377  13.166  1.00 20.00           C
ATOM    122  N4    C B   2      20.101  15.584  12.843  1.00 20.00           N
ATOM    123  C5    C B   2      21.776  14.146  13.764  1.00 20.00           C
ATOM    124  C6    C B   2      23.072  13.996  14.065  1.00 20.00           C
ATOM      0  H5'   C B   2      26.014  11.112  15.401  1.00 77.93           H
ATOM      0 H5''   C B   2      26.061  11.559  16.894  1.00 77.93           H
ATOM      0  H4'   C B   2      27.131  13.078  15.693  1.00 78.66           H
ATOM      0  H3'   C B   2      24.770  14.191  16.771  1.00 78.23           H
ATOM      0  H2'   C B   2      25.239  16.248  15.720  1.00 79.34           H
ATOM      0 HO2'   C B   2      27.522  15.780  16.107  1.00 79.34           H
ATOM      0  H1'   C B   2      25.929  15.233  13.434  1.00 78.22           H
ATOM      0  H41   C B   2      19.515  14.974  12.997  1.00 20.00           H
ATOM      0  H42   C B   2      19.864  16.327  12.481  1.00 20.00           H
ATOM      0  H5    C B   2      21.158  13.473  13.938  1.00 20.00           H
ATOM      0  H6    C B   2      23.361  13.203  14.455  1.00 20.00           H
ATOM    125  P     G B   3      25.668  14.962  18.979  1.00 76.88           P
ATOM    126  OP1   G B   3      26.614  14.789  20.110  1.00 78.59           O
ATOM    127  OP2   G B   3      24.339  14.301  19.026  1.00 75.20           O
ATOM    128  O5'   G B   3      25.435  16.524  18.767  1.00 78.76           O
ATOM    129  C5'   G B   3      26.534  17.409  18.606  1.00 77.93           C
ATOM    130  C4'   G B   3      26.078  18.786  18.191  1.00 78.66           C
ATOM    131  O4'   G B   3      25.414  18.715  16.903  1.00 79.49           O
ATOM    132  C3'   G B   3      25.053  19.449  19.099  1.00 78.23           C
ATOM    133  O3'   G B   3      25.635  20.055  20.240  1.00 77.17           O
ATOM    134  C2'   G B   3      24.371  20.437  18.163  1.00 79.34           C
ATOM    135  O2'   G B   3      25.175  21.593  17.980  1.00 79.34           O
ATOM    136  C1'   G B   3      24.355  19.649  16.853  1.00 78.22           C
ATOM    137  N9    G B   3      23.091  18.912  16.667  1.00 20.00           N
ATOM    138  C8    G B   3      22.884  17.565  16.841  1.00 20.00           C
ATOM    139  N7    G B   3      21.654  17.200  16.601  1.00 20.00           N
ATOM    140  C5    G B   3      21.009  18.377  16.246  1.00 20.00           C
ATOM    141  C6    G B   3      19.659  18.610  15.877  1.00 20.00           C
ATOM    142  O6    G B   3      18.733  17.795  15.786  1.00 20.00           O
ATOM    143  N1    G B   3      19.431  19.953  15.597  1.00 20.00           N
ATOM    144  C2    G B   3      20.377  20.946  15.663  1.00 20.00           C
ATOM    145  N2    G B   3      19.960  22.182  15.354  1.00 20.00           N
ATOM    146  N3    G B   3      21.639  20.742  16.005  1.00 20.00           N
ATOM    147  C4    G B   3      21.883  19.444  16.282  1.00 20.00           C
ATOM      0  H5'   G B   3      27.142  17.054  17.939  1.00 77.93           H
ATOM      0 H5''   G B   3      27.028  17.468  19.438  1.00 77.93           H
ATOM      0  H4'   G B   3      26.897  19.305  18.207  1.00 78.66           H
ATOM      0  H3'   G B   3      24.430  18.821  19.497  1.00 78.23           H
ATOM      0  H2'   G B   3      23.506  20.739  18.483  1.00 79.34           H
ATOM      0 HO2'   G B   3      25.658  21.709  18.658  1.00 79.34           H
ATOM      0  H1'   G B   3      24.448  20.271  16.115  1.00 78.22           H
ATOM      0  H8    G B   3      23.552  16.974  17.103  1.00 20.00           H
ATOM      0  H1    G B   3      18.635  20.179  15.364  1.00 20.00           H
ATOM      0  H21   G B   3      20.512  22.841  15.379  1.00 20.00           H
ATOM      0  H22   G B   3      19.140  22.316  15.131  1.00 20.00           H
ATOM    148  P     U B   4      24.901  19.965  21.666  1.00 76.88           P
ATOM    149  OP1   U B   4      25.844  20.460  22.701  1.00 78.59           O
ATOM    150  OP2   U B   4      24.313  18.608  21.797  1.00 75.20           O
ATOM    151  O5'   U B   4      23.704  21.011  21.551  1.00 78.76           O
ATOM    152  C5'   U B   4      23.961  22.392  21.346  1.00 77.93           C
ATOM    153  C4'   U B   4      22.696  23.148  21.023  1.00 78.66           C
ATOM    154  O4'   U B   4      22.134  22.654  19.779  1.00 79.49           O
ATOM    155  C3'   U B   4      21.560  23.004  22.024  1.00 78.23           C
ATOM    156  O3'   U B   4      21.704  23.850  23.151  1.00 77.17           O
ATOM    157  C2'   U B   4      20.331  23.315  21.180  1.00 79.34           C
ATOM    158  O2'   U B   4      20.183  24.716  21.005  1.00 79.34           O
ATOM    159  C1'   U B   4      20.724  22.700  19.835  1.00 78.22           C
ATOM    160  N1    U B   4      20.198  21.324  19.673  1.00 20.00           N
ATOM    161  C2    U B   4      18.881  21.183  19.278  1.00 20.00           C
ATOM    162  O2    U B   4      18.147  22.133  19.064  1.00 20.00           O
ATOM    163  N3    U B   4      18.452  19.886  19.143  1.00 20.00           N
ATOM    164  C4    U B   4      19.189  18.740  19.358  1.00 20.00           C
ATOM    165  O4    U B   4      18.660  17.639  19.194  1.00 20.00           O
ATOM    166  C5    U B   4      20.542  18.969  19.762  1.00 20.00           C
ATOM    167  C6    U B   4      20.988  20.222  19.901  1.00 20.00           C
ATOM      0  H5'   U B   4      24.598  22.500  20.622  1.00 77.93           H
ATOM      0 H5''   U B   4      24.370  22.768  22.141  1.00 77.93           H
ATOM      0  H4'   U B   4      22.988  24.073  21.008  1.00 78.66           H
ATOM      0  H3'   U B   4      21.519  22.125  22.433  1.00 78.23           H
ATOM      0  H2'   U B   4      19.504  22.987  21.566  1.00 79.34           H
ATOM      0 HO2'   U B   4      20.622  25.120  21.596  1.00 79.34           H
ATOM      0 HO3'   U B   4      21.055  23.736  23.672  1.00 77.17           H
ATOM      0  H1'   U B   4      20.348  23.244  19.125  1.00 78.22           H
ATOM      0  H3    U B   4      17.635  19.776  18.898  1.00 20.00           H
ATOM      0  H5    U B   4      21.109  18.250  19.927  1.00 20.00           H
ATOM      0  H6    U B   4      21.871  20.354  20.163  1.00 20.00           H
TER
END
"""

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "OK. Time: %8.3f"%(time.time()-t0)

