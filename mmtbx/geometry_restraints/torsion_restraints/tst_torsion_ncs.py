from __future__ import absolute_import, division, print_function

import iotbx.pdb
from cctbx.array_family import flex
from mmtbx.monomer_library import pdb_interpretation
from cctbx import adp_restraints # import dependency
from mmtbx.geometry_restraints.torsion_restraints import torsion_ncs, utils
from six.moves import cStringIO as StringIO
import mmtbx
import mmtbx.model

pdb_str_1 = """
CRYST1   46.053    9.561   20.871  90.00  97.43  90.00 C 1 2 1       8
ATOM      1  N   LYS A   1       8.669   0.060   3.795  1.00 12.85           N
ATOM      2  CA  LYS A   1       7.501  -0.641   4.336  1.00 11.37           C
ATOM      3  C   LYS A   1       6.186   0.030   3.984  1.00 11.88           C
ATOM      4  O   LYS A   1       6.176   1.208   3.641  1.00 11.02           O
ATOM      5  CB  LYS A   1       7.606  -0.792   5.871  1.00 15.61           C
ATOM      6  CG  LYS A   1       7.543   0.531   6.652  1.00 23.32           C
ATOM      7  CD  LYS A   1       8.894   0.917   7.175  1.00 22.74           C
ATOM      8  CE  LYS A   1       8.851   2.260   7.839  1.00 22.96           C
ATOM      9  NZ  LYS A   1      10.091   2.510   8.611  1.00 18.89           N
ATOM     10  N   LEU A   2       5.063  -0.702   4.139  1.00  8.33           N
ATOM     11  CA  LEU A   2       3.720  -0.164   3.917  1.00  8.69           C
ATOM     12  C   LEU A   2       2.701  -0.803   4.885  1.00  9.37           C
ATOM     13  O   LEU A   2       2.675  -2.017   5.019  1.00  8.16           O
ATOM     14  CB  LEU A   2       3.331  -0.407   2.441  1.00  9.53           C
ATOM     15  CG  LEU A   2       1.966   0.033   1.872  1.00 16.06           C
ATOM     16  CD1 LEU A   2       0.773  -0.824   2.351  1.00 15.86           C
ATOM     17  CD2 LEU A   2       1.747   1.479   1.927  1.00 17.67           C
ATOM     18  N   VAL A   3       1.856   0.019   5.527  1.00  6.29           N
ATOM     19  CA  VAL A   3       0.748  -0.404   6.386  1.00  6.45           C
ATOM     20  C   VAL A   3      -0.501   0.203   5.783  1.00 11.86           C
ATOM     21  O   VAL A   3      -0.569   1.425   5.609  1.00 10.67           O
ATOM     22  CB  VAL A   3       0.865   0.064   7.854  1.00 10.59           C
ATOM     23  CG1 VAL A   3      -0.285  -0.493   8.710  1.00 10.71           C
ATOM     24  CG2 VAL A   3       2.219  -0.288   8.450  1.00 10.64           C
ATOM     25  N   PHE A   4      -1.510  -0.637   5.537  1.00  8.41           N
ATOM     26  CA  PHE A   4      -2.797  -0.223   4.998  1.00  8.51           C
ATOM     27  C   PHE A   4      -3.900  -0.807   5.880  1.00 11.98           C
ATOM     28  O   PHE A   4      -3.871  -2.003   6.216  1.00 10.68           O
ATOM     29  CB  PHE A   4      -2.954  -0.769   3.564  1.00  9.14           C
ATOM     30  CG  PHE A   4      -4.317  -0.526   2.953  1.00  9.69           C
ATOM     31  CD1 PHE A   4      -4.536   0.563   2.127  1.00 10.70           C
ATOM     32  CD2 PHE A   4      -5.362  -1.422   3.159  1.00 12.10           C
ATOM     33  CE1 PHE A   4      -5.787   0.791   1.562  1.00 11.47           C
ATOM     34  CE2 PHE A   4      -6.618  -1.195   2.592  1.00 15.42           C
ATOM     35  CZ  PHE A   4      -6.815  -0.096   1.778  1.00 13.27           C
ATOM     36  N   PHE A   5      -4.902   0.009   6.183  1.00  8.96           N
ATOM     37  CA  PHE A   5      -6.081  -0.452   6.908  1.00  9.67           C
ATOM     38  C   PHE A   5      -7.304   0.208   6.324  1.00 13.65           C
ATOM     39  O   PHE A   5      -7.349   1.428   6.235  1.00 10.13           O
ATOM     40  CB  PHE A   5      -5.986  -0.239   8.445  1.00 11.19           C
ATOM     41  CG  PHE A   5      -7.306  -0.481   9.152  1.00 12.79           C
ATOM     42  CD1 PHE A   5      -7.716  -1.766   9.473  1.00 16.12           C
ATOM     43  CD2 PHE A   5      -8.156   0.580   9.457  1.00 15.48           C
ATOM     44  CE1 PHE A   5      -8.953  -1.987  10.088  1.00 17.28           C
ATOM     45  CE2 PHE A   5      -9.379   0.359  10.098  1.00 18.56           C
ATOM     46  CZ  PHE A   5      -9.771  -0.922  10.402  1.00 15.89           C
ATOM     47  N   ALA A   6      -8.282  -0.608   5.899  1.00 12.21           N
ATOM     48  CA  ALA A   6      -9.570  -0.134   5.389  1.00 15.49           C
ATOM     49  C   ALA A   6     -10.647  -0.909   6.127  1.00 31.54           C
ATOM     50  O   ALA A   6     -10.637  -2.159   6.058  1.00 33.30           O
ATOM     51  CB  ALA A   6      -9.677  -0.381   3.898  1.00 16.16           C
ATOM     52  OXT ALA A   6     -11.417  -0.275   6.874  1.00 51.39           O
TER      53      ALA A   6
ATOM     54  N   LYS B   1     -10.065  -4.813   6.193  1.00 12.06           N
ATOM     55  CA  LYS B   1      -8.858  -5.207   5.451  1.00 12.29           C
ATOM     56  C   LYS B   1      -7.574  -4.609   6.062  1.00 12.37           C
ATOM     57  O   LYS B   1      -7.486  -3.387   6.205  1.00  7.97           O
ATOM     58  CB  LYS B   1      -8.985  -4.776   3.968  1.00 15.40           C
ATOM     59  CG  LYS B   1      -7.989  -5.451   3.027  1.00 27.98           C
ATOM     60  CD  LYS B   1      -8.475  -5.495   1.572  1.00 36.65           C
ATOM     61  CE  LYS B   1      -8.321  -4.183   0.829  1.00 51.92           C
ATOM     62  NZ  LYS B   1      -9.598  -3.419   0.783  1.00 65.39           N
ATOM     63  N   LEU B   2      -6.567  -5.463   6.352  1.00  8.53           N
ATOM     64  CA  LEU B   2      -5.275  -5.044   6.904  1.00  8.19           C
ATOM     65  C   LEU B   2      -4.174  -5.652   6.041  1.00  9.32           C
ATOM     66  O   LEU B   2      -4.161  -6.876   5.835  1.00  8.36           O
ATOM     67  CB  LEU B   2      -5.179  -5.546   8.364  1.00  9.24           C
ATOM     68  CG  LEU B   2      -4.127  -5.040   9.352  1.00 14.56           C
ATOM     69  CD1 LEU B   2      -2.893  -5.867   9.314  1.00 16.64           C
ATOM     70  CD2 LEU B   2      -3.873  -3.540   9.265  1.00 13.86           C
ATOM     71  N   VAL B   3      -3.293  -4.804   5.481  1.00  5.98           N
ATOM     72  CA  VAL B   3      -2.164  -5.264   4.671  1.00  4.31           C
ATOM     73  C   VAL B   3      -0.910  -4.581   5.178  1.00  8.02           C
ATOM     74  O   VAL B   3      -0.872  -3.349   5.309  1.00  7.46           O
ATOM     75  CB  VAL B   3      -2.312  -5.033   3.142  1.00  6.99           C
ATOM     76  CG1 VAL B   3      -1.206  -5.781   2.368  1.00  6.50           C
ATOM     77  CG2 VAL B   3      -3.696  -5.419   2.642  1.00  7.17           C
ATOM     78  N   PHE B   4       0.131  -5.371   5.389  1.00  4.71           N
ATOM     79  CA  PHE B   4       1.445  -4.868   5.777  1.00  5.68           C
ATOM     80  C   PHE B   4       2.526  -5.586   4.957  1.00  9.59           C
ATOM     81  O   PHE B   4       2.440  -6.805   4.790  1.00  9.76           O
ATOM     82  CB  PHE B   4       1.661  -5.087   7.288  1.00  7.35           C
ATOM     83  CG  PHE B   4       3.098  -5.169   7.735  1.00 10.11           C
ATOM     84  CD1 PHE B   4       3.900  -4.033   7.765  1.00 14.40           C
ATOM     85  CD2 PHE B   4       3.655  -6.388   8.121  1.00 13.08           C
ATOM     86  CE1 PHE B   4       5.239  -4.115   8.153  1.00 16.51           C
ATOM     87  CE2 PHE B   4       4.988  -6.466   8.532  1.00 17.07           C
ATOM     88  CZ  PHE B   4       5.770  -5.324   8.557  1.00 16.28           C
ATOM     89  N   PHE B   5       3.546  -4.840   4.448  1.00  7.84           N
ATOM     90  CA  PHE B   5       4.703  -5.447   3.782  1.00  6.46           C
ATOM     91  C   PHE B   5       5.938  -4.609   3.941  1.00 10.69           C
ATOM     92  O   PHE B   5       5.830  -3.420   4.235  1.00  8.77           O
ATOM     93  CB  PHE B   5       4.457  -5.846   2.320  1.00  8.36           C
ATOM     94  CG  PHE B   5       3.975  -4.763   1.395  1.00  8.12           C
ATOM     95  CD1 PHE B   5       4.873  -3.877   0.812  1.00 10.77           C
ATOM     96  CD2 PHE B   5       2.629  -4.640   1.092  1.00  8.53           C
ATOM     97  CE1 PHE B   5       4.421  -2.862  -0.036  1.00 11.77           C
ATOM     98  CE2 PHE B   5       2.170  -3.607   0.275  1.00 11.24           C
ATOM     99  CZ  PHE B   5       3.073  -2.744  -0.307  1.00 10.01           C
ATOM    100  N   ALA B   6       7.120  -5.247   3.770  1.00 10.09           N
ATOM    101  CA  ALA B   6       8.448  -4.651   3.900  1.00 12.61           C
ATOM    102  C   ALA B   6       9.483  -5.594   3.319  1.00 36.73           C
ATOM    103  O   ALA B   6       9.457  -6.801   3.661  1.00 40.97           O
ATOM    104  CB  ALA B   6       8.760  -4.382   5.366  1.00 13.92           C
ATOM    105  OXT ALA B   6      10.304  -5.128   2.503  1.00 60.32           O
TER     106      ALA B   6
"""

def exercise_1(mon_lib_srv, ener_lib):
  pdb_h = iotbx.pdb.input(lines=pdb_str_1, source_info=None).construct_hierarchy()
  log = StringIO()
  dihedral_proxies = utils.get_complete_dihedral_proxies(
                       pdb_hierarchy=pdb_h)
  assert len(dihedral_proxies) == 54, \
      "Expected 54, got %d" % len(dihedral_proxies)

  # default run (1 residue is out of NCS)
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.ncs_search.enabled=True
  pdb_inp = iotbx.pdb.input(source_info=None, lines=flex.split_lines(pdb_str_1))
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.process(pdb_interpretation_params=params,
    make_restraints=True)
  ncs_manager = torsion_ncs.torsion_ncs(
                  model = model,
                  log=log)
  nprox = ncs_manager.get_n_proxies()
  assert nprox == 28, "got %d instead of 28" % nprox


  # supply full NCS
  cuspars = iotbx.phil.parse("""
pdb_interpretation.ncs_search.enabled=True
pdb_interpretation.ncs_group {
  reference        = (chain A )
  selection        = (chain B )
}
""")
  params = mmtbx.model.manager.get_default_pdb_interpretation_scope()
  p = params.fetch(cuspars).extract()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=flex.split_lines(pdb_str_1))
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.process(pdb_interpretation_params=p, make_restraints=True)
  ncs_manager = torsion_ncs.torsion_ncs(
                  model = model,
                  log=log)
  nprox = ncs_manager.get_n_proxies()
  assert nprox == 40, "got %d instead of 40" % nprox

if (__name__ == "__main__"):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  exercise_1(mon_lib_srv, ener_lib)
  print("OK")
