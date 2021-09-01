from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from six.moves import zip

#-----------------------------------------------------------------------------
# This test checks the parameterization of H/D atoms in models containing
# both H and D atoms
#-----------------------------------------------------------------------------

def exercise(pdb_str):
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.use_neutron_distances = True

  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(pdb_interpretation_params=params, make_restraints=True)

  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_parameterization = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart         = sites_cart,
    threshold          = 0.05)
  h_distances        = diagnostics.h_distances

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_parameterization) - h_parameterization.count(None)

  if (pdb_str != pdb_str_02):
    assert (number_h_para == number_h), 'Not all H atoms are parameterized'

# For each H atom, check if distance between computed H and that in input model is
# not too large
  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.1), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_parameterization[ih].htype, atoms[ih].name, ih, labels.resseq.strip())

  #for type1, type2 in zip(type_list, type_list_known):
  #  assert (type1 == type2)


pdb_str_00 = """\
CRYST1   14.815   23.523   21.068  90.00  90.00  90.00 P 1
SCALE1      0.067499  0.000000  0.000000        0.00000
SCALE2      0.000000  0.042512  0.000000        0.00000
SCALE3      0.000000  0.000000  0.047465        0.00000
ATOM      1  N   ARG A 161       5.235  16.231  12.444  1.00  7.68           N
ATOM      2  CA  ARG A 161       6.309  16.550  13.364  1.00  9.47           C
ATOM      3  C   ARG A 161       7.140  17.751  12.856  1.00 10.84           C
ATOM      4  O   ARG A 161       7.669  18.523  13.658  1.00 11.74           O
ATOM      5  CB  ARG A 161       7.186  15.310  13.584  1.00 13.82           C
ATOM      6  CG  ARG A 161       6.357  14.097  14.068  1.00 16.56           C
ATOM      7  CD  ARG A 161       7.166  13.049  14.764  1.00 11.41           C
ATOM      8  NE  ARG A 161       6.350  11.919  15.227  1.00 18.74           N
ATOM      9  CZ  ARG A 161       6.279  10.732  14.648  1.00 12.24           C
ATOM     10  NH1 ARG A 161       6.967  10.428  13.537  1.00 11.51           N
ATOM     11  NH2 ARG A 161       5.503   9.814  15.194  1.00 16.77           N
ATOM     12  DE  ARG A 161       5.791  12.059  16.068  1.00 22.49           D
ATOM     13 DH11 ARG A 161       7.580  11.112  13.087  1.00 13.81           D
ATOM     14 DH21 ARG A 161       5.000  10.019  16.043  1.00 20.13           D
ATOM     15  HA  ARG A 161       5.875  16.826  14.327  1.00 11.37           H
ATOM     16  HB2 ARG A 161       7.666  15.039  12.643  1.00 16.58           H
ATOM     17  HB3 ARG A 161       7.940  15.534  14.339  1.00 16.58           H
ATOM     18  HG2 ARG A 161       5.598  14.449  14.766  1.00 19.87           H
ATOM     19  HG3 ARG A 161       5.878  13.631  13.205  1.00 19.87           H
ATOM     20  HD2 ARG A 161       7.917  12.663  14.076  1.00 13.70           H
ATOM     21  HD3 ARG A 161       7.651  13.492  15.635  1.00 13.70           H
ATOM     22 HH12AARG A 161       6.870   9.515  13.114  0.26 13.81           H
ATOM     23 HH22AARG A 161       5.423   8.890  14.766  0.54 20.13           H
ATOM     24 DH12BARG A 161       6.870   9.515  13.114  0.74 13.81           D
ATOM     25 DH22BARG A 161       5.423   8.890  14.766  0.46 20.13           D
TER
ATOM     26  N   GLY B   2      12.057   9.922   8.317  1.00 27.73           N
ATOM     27  CA  GLY B   2      11.194  10.314   7.218  1.00 26.03           C
ATOM     28  C   GLY B   2       9.719  10.240   7.559  1.00 24.25           C
ATOM     29  O   GLY B   2       9.041  11.263   7.653  1.00 20.41           O
ATOM     30  N   ILE B   3       9.220   9.022   7.744  1.00 27.73           N
ATOM     31  CA  ILE B   3       7.818   8.811   8.076  1.00 26.03           C
ATOM     32  C   ILE B   3       7.565   9.227   9.520  1.00 24.25           C
ATOM     33  O   ILE B   3       8.276   8.806  10.431  1.00 20.41           O
ATOM     34  CB  ILE B   3       7.404   7.345   7.845  1.00 27.04           C
ATOM     35  CG1 ILE B   3       7.452   7.006   6.354  1.00 28.51           C
ATOM     36  CG2 ILE B   3       6.009   7.089   8.398  1.00 26.88           C
ATOM     37  CD1 ILE B   3       7.628   5.530   6.070  1.00 28.93           C
ATOM     38  DA  ILE B   3       7.203   9.439   7.432  1.00 26.09           D
ATOM     39  DB  ILE B   3       8.109   6.700   8.370  1.00 25.53           D
ATOM     40 DD11 ILE B   3       8.556   5.190   6.530  1.00 27.61           D
ATOM     41 DG12 ILE B   3       6.516   7.323   5.893  1.00 27.88           D
ATOM     42  D  AILE B   3       9.763   8.161   7.671  0.96 27.94           D
ATOM     43 DD12AILE B   3       7.670   5.378   4.992  0.43 28.45           D
ATOM     44 DD13AILE B   3       6.783   4.984   6.490  0.86 28.24           D
ATOM     45 DG13AILE B   3       8.288   7.536   5.896  0.69 28.21           D
ATOM     46 DG21AILE B   3       5.354   7.909   8.103  0.96 27.70           D
ATOM     47 DG22AILE B   3       6.062   7.027   9.485  0.78 27.54           D
ATOM     48 DG23AILE B   3       5.635   6.149   7.993  0.81 27.54           D
ATOM     49  H  BILE B   3       9.763   8.161   7.671  0.04 27.94           H
ATOM     50 HG13BILE B   3       8.288   7.536   5.896  0.31 28.21           H
ATOM     51 HG21BILE B   3       5.354   7.909   8.103  0.04 27.70           H
ATOM     52 HG22BILE B   3       6.062   7.027   9.485  0.22 27.54           H
ATOM     53 HG23BILE B   3       5.635   6.149   7.993  0.19 27.54           H
ATOM     54 HD12BILE B   3       7.670   5.378   4.992  0.57 28.45           H
ATOM     55 HD13BILE B   3       6.783   4.984   6.490  0.14 28.24           H
TER
"""

pdb_str_01 = """
CRYST1   14.630   14.210   14.547  90.00  90.00  90.00 P 1
SCALE1      0.068353  0.000000  0.000000        0.00000
SCALE2      0.000000  0.070373  0.000000        0.00000
SCALE3      0.000000  0.000000  0.068743        0.00000
ATOM   1387  N   ILE A  63      21.154 -27.868 -37.969  1.00 22.34           N
ATOM   1388  CA  ILE A  63      21.399 -26.570 -38.585  1.00 21.95           C
ATOM   1389  C   ILE A  63      20.231 -25.648 -38.320  1.00 23.15           C
ATOM   1390  O   ILE A  63      19.079 -25.992 -38.576  1.00 24.15           O
ATOM   1391  CB  ILE A  63      21.593 -26.710 -40.099  1.00 23.08           C
ATOM   1392  CG1 ILE A  63      22.822 -27.578 -40.395  1.00 23.79           C
ATOM   1393  CG2 ILE A  63      21.745 -25.345 -40.744  1.00 23.03           C
ATOM   1394  CD1 ILE A  63      22.871 -28.137 -41.811  1.00 25.99           C
ATOM   1397  DA AILE A  63      22.287 -26.154 -38.157  0.44 22.47           D
ATOM   1398  HA BILE A  63      22.287 -26.154 -38.157  0.56 22.47           H
ATOM   1399  DB AILE A  63      20.719 -27.170 -40.493  0.99 23.45           D
ATOM   1400  HB BILE A  63      20.719 -27.170 -40.493  0.01 23.45           H
ATOM   1401 DG12AILE A  63      23.709 -26.976 -40.253  0.54 23.93           D
ATOM   1402 HG12BILE A  63      23.709 -26.976 -40.253  0.46 23.93           H
ATOM   1403 DG13AILE A  63      22.850 -28.402 -39.700  0.65 23.97           D
ATOM   1404 HG13BILE A  63      22.850 -28.402 -39.700  0.35 23.97           H
ATOM   1405 HG21 ILE A  63      21.066 -24.629 -40.306  1.00 22.80           H
ATOM   1406 DG22AILE A  63      21.541 -25.426 -41.799  0.51 22.86           D
ATOM   1407 HG22BILE A  63      21.541 -25.426 -41.799  0.49 22.86           H
ATOM   1408 DG23 ILE A  63      22.758 -25.015 -40.597  1.00 23.40           D
ATOM   1409 DD11AILE A  63      21.935 -28.641 -42.033  0.67 24.71           D
ATOM   1410 HD11BILE A  63      21.935 -28.641 -42.033  0.33 24.71           H
ATOM   1411 DD12AILE A  63      23.680 -28.839 -41.891  0.44 24.81           D
ATOM   1412 HD12BILE A  63      23.680 -28.839 -41.891  0.56 24.81           H
ATOM   1413 DD13AILE A  63      23.024 -27.344 -42.516  0.90 24.27           D
ATOM   1414 HD13BILE A  63      23.024 -27.344 -42.516  0.10 24.27           H
"""

pdb_str_02 = """
CRYST1   15.636   16.098   18.562  90.00  90.00  90.00 P 1
SCALE1      0.063955  0.000000  0.000000        0.00000
SCALE2      0.000000  0.062120  0.000000        0.00000
SCALE3      0.000000  0.000000  0.053874        0.00000
ATOM    517  N   GLU A  30       4.726  39.084   1.924  1.00 18.34           N
ATOM    519  C   GLU A  30       6.556  40.664   2.421  1.00 20.67           C
ATOM    533  D   GLU A  31       5.082  40.344   3.856  1.00 21.59           D
ATOM    534  N  AGLU A  31       5.986  40.754   3.619  0.29 20.27           N
ATOM    535  CA AGLU A  31       6.623  41.467   4.722  0.29 21.76           C
ATOM    536  C  AGLU A  31       7.429  40.530   5.617  0.29 20.64           C
ATOM    548  N  CGLU A  31       5.986  40.753   3.620  0.71 20.17           N
ATOM    549  CA CGLU A  31       6.621  41.462   4.726  0.71 21.77           C
ATOM    550  C  CGLU A  31       7.488  40.540   5.574  0.71 20.61           C
"""

pdb_list = [pdb_str_00, pdb_str_01, pdb_str_02]

pdb_list_name = ['pdb_str_00', 'pdb_str_01', 'pdb_str_02']

#type_list_known = ['2tetra', '2tetra', 'alg1b', '3neigbs', '3neigbs',
#  '3neigbs', 'alg1b', '2tetra', '2tetra', '2tetra', '2tetra', 'alg1b',
#  '3neigbs', '2tetra', '2tetra', '3neigbs', '2tetra', '2tetra',
#  'flat_2neigbs', 'alg1b', '3neigbs', '2tetra', '2tetra', 'alg1b',
#  '2tetra', '2tetra', 'alg1b', '2tetra', '2tetra', 'alg1b', '3neigbs',
#  'prop', 'prop', 'prop', 'prop', 'prop', 'prop', '3neigbs', 'prop',
#  'prop', 'prop', 'prop', 'prop', 'prop', '3neigbs', '2tetra', '2tetra',
#  'alg1b', '2tetra', '2tetra', 'alg1b', 'flat_2neigbs']

def run():
  for pdb_str, str_name in zip(pdb_list,pdb_list_name):
    exercise(pdb_str=pdb_str)


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
