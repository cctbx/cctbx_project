from __future__ import absolute_import, division, print_function
import mmtbx.monomer_library
from scitbx.array_family import flex
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
from libtbx.test_utils import approx_equal

pdb_str1 = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
SCALE1      0.058962  0.000000  0.000000        0.00000
SCALE2      0.000000  0.051401  0.000000        0.00000
SCALE3      0.000000  0.000000  0.050401        0.00000
ATOM      1  N   TYR A  20      10.418  10.903   8.137  1.00 23.81           N
ATOM      2  CA  TYR A  20      10.810   9.981   7.034  1.00 35.69           C
ATOM      3  C   TYR A  20       9.703   8.988   6.760  1.00 28.94           C
ATOM      4  O   TYR A  20       9.364   8.158   7.667  1.00 25.26           O
ATOM      5  CB  TYR A  20      12.115   9.280   7.395  1.00 33.37           C
ATOM      6  CG  TYR A  20      13.316  10.201   7.509  1.00 31.48           C
ATOM      7  CD1 TYR A  20      14.071  10.517   6.393  1.00 36.91           C
ATOM      8  CD2 TYR A  20      13.704  10.721   8.703  1.00 21.36           C
ATOM      9  CE1 TYR A  20      15.188  11.317   6.471  1.00 36.56           C
ATOM     10  CE2 TYR A  20      14.835  11.519   8.817  1.00 29.56           C
ATOM     11  CZ  TYR A  20      15.540  11.837   7.694  1.00 35.08           C
ATOM     12  OH  TYR A  20      16.672  12.664   7.781  1.00 38.73           O
ATOM     13  N  AARG A  21       9.160   8.989   5.577  0.50 38.95           N
ATOM     14  CA AARG A  21       8.073   8.140   5.126  0.50 38.77           C
ATOM     15  C  AARG A  21       6.863   8.252   6.048  0.50 27.69           C
ATOM     16  O  AARG A  21       6.208   7.218   6.352  0.50 22.82           O
ATOM     17  CB AARG A  21       8.542   6.674   5.053  0.50 20.00           C
ATOM     18  CG AARG A  21       7.471   5.695   4.562  0.50 20.00           C
ATOM     19  CD AARG A  21       8.021   4.282   4.444  0.50 20.00           C
ATOM     20  NE AARG A  21       6.983   3.348   3.946  0.50 20.00           N
ATOM     21  CZ AARG A  21       5.878   3.748   3.337  0.50 20.00           C
ATOM     22  NH1AARG A  21       5.643   5.047   3.127  0.50 20.00           N
ATOM     23  NH2AARG A  21       4.997   2.854   2.938  0.50 20.00           N
ATOM     13  N  BTRP A  21       9.160   8.989   5.577  0.50 38.95           N
ATOM     14  CA BTRP A  21       8.073   8.140   5.126  0.50 38.77           C
ATOM     15  C  BTRP A  21       6.863   8.252   6.048  0.50 27.69           C
ATOM     16  O  BTRP A  21       6.208   7.218   6.352  0.50 22.82           O
ATOM     17  CB BTRP A  21       8.532   6.683   5.040  0.50 20.00           C
ATOM     18  CG BTRP A  21       9.522   6.431   3.944  0.50 20.00           C
ATOM     19  CD1BTRP A  21      10.883   6.417   4.051  0.50 20.00           C
ATOM     20  CD2BTRP A  21       9.228   6.154   2.569  0.50 20.00           C
ATOM     21  NE1BTRP A  21      11.453   6.152   2.830  0.50 20.00           N
ATOM     22  CE2BTRP A  21      10.459   5.985   1.902  0.50 20.00           C
ATOM     23  CE3BTRP A  21       8.047   6.033   1.829  0.50 20.00           C
ATOM     24  CZ2BTRP A  21      10.541   5.702   0.540  0.50 20.00           C
ATOM     25  CZ3BTRP A  21       8.131   5.752   0.478  0.50 20.00           C
ATOM     26  CH2BTRP A  21       9.368   5.590  -0.152  0.50 20.00           C
ATOM     29  N  CHIS A  21       9.160   8.989   5.577  1.00 38.95           N
ATOM     30  CA CHIS A  21       8.073   8.140   5.126  1.00 38.77           C
ATOM     31  C  CHIS A  21       6.863   8.252   6.048  1.00 27.69           C
ATOM     32  O  CHIS A  21       6.208   7.218   6.352  1.00 22.82           O
ATOM     33  CB CHIS A  21       8.533   6.684   5.037  1.00 20.00           C
ATOM     34  CG CHIS A  21       9.008   6.118   6.338  1.00 20.00           C
ATOM     35  ND1CHIS A  21       8.412   6.420   7.544  1.00 20.00           N
ATOM     36  CD2CHIS A  21      10.019   5.265   6.623  1.00 20.00           C
ATOM     37  CE1CHIS A  21       9.038   5.781   8.515  1.00 20.00           C
ATOM     38  NE2CHIS A  21      10.017   5.071   7.985  1.00 20.00           N
ATOM     32  N   GLY A  22       6.554   9.436   6.463  1.00 24.85           N
ATOM     33  CA  GLY A  22       5.410   9.707   7.327  1.00 29.53           C
ATOM     34  C   GLY A  22       5.748   9.622   8.821  1.00 33.22           C
ATOM     35  O   GLY A  22       5.479  10.518   9.577  1.00 30.06           O
ATOM     36  N   TYR A  23       6.330   8.487   9.201  1.00 27.25           N
ATOM     37  CA  TYR A  23       6.695   8.251  10.596  1.00 34.16           C
ATOM     38  C   TYR A  23       7.852   9.136  11.044  1.00 23.48           C
ATOM     39  O   TYR A  23       8.641   9.570  10.203  1.00 39.30           O
ATOM     40  CB  TYR A  23       7.052   6.786  10.831  1.00 29.65           C
ATOM     41  CG  TYR A  23       5.903   5.834  10.636  1.00 34.88           C
ATOM     42  CD1 TYR A  23       5.022   5.540  11.687  1.00 30.77           C
ATOM     43  CD2 TYR A  23       5.649   5.271   9.396  1.00 32.29           C
ATOM     44  CE1 TYR A  23       3.990   4.674  11.505  1.00 39.91           C
ATOM     45  CE2 TYR A  23       4.589   4.378   9.190  1.00 30.45           C
ATOM     46  CZ  TYR A  23       3.770   4.107  10.257  1.00 37.13           C
ATOM     47  OH  TYR A  23       2.693   3.226  10.097  1.00 38.59           O
ATOM     48  N   SER A  24       7.948   9.408  12.310  1.00 36.25           N
ATOM     49  CA  SER A  24       9.013  10.264  12.855  1.00 26.44           C
ATOM     50  C   SER A  24       9.607   9.636  14.109  1.00 20.40           C
ATOM     51  O   SER A  24       8.916   9.382  15.096  1.00 39.27           O
ATOM     52  CB  SER A  24       8.482  11.655  13.176  1.00 21.07           C
ATOM     53  OG  SER A  24       7.498  11.601  14.206  1.00 34.95           O
TER      54      SER A  24
END
"""

pdb_str2 = """\
ATOM     28  N   PRO C  29      61.293  16.365  38.366  1.00 45.42           N
ATOM     29  CA  PRO C  29      61.610  17.144  39.554  1.00 51.27           C
ATOM     30  C   PRO C  29      60.399  17.813  40.109  1.00 50.43           C
ATOM     31  O   PRO C  29      59.324  17.664  39.560  1.00 47.67           O
ATOM     32  CB  PRO C  29      62.495  18.235  38.993  1.00 49.45           C
ATOM     33  CG  PRO C  29      61.908  18.507  37.705  1.00 55.39           C
ATOM     34  CD  PRO C  29      61.575  17.152  37.157  1.00 56.41           C
"""

pdb_str3 = """\
ATOM    127  N   PRO L  18      18.596 -78.410  66.845  1.00 33.30           N
ATOM    128  CA  PRO L  18      19.467 -77.234  66.790  1.00 32.71           C
ATOM    129  C   PRO L  18      20.924 -77.585  66.494  1.00 30.96           C
ATOM    130  O   PRO L  18      21.369 -78.705  66.735  1.00 30.74           O
ATOM    131  CB  PRO L  18      19.289 -76.622  68.173  1.00 31.82           C
ATOM    132  CG  PRO L  18      19.175 -77.838  69.037  1.00 33.36           C
ATOM    133  CD  PRO L  18      18.236 -78.728  68.240  1.00 34.01           C
ATOM    134  N   ALA L  19      21.657 -76.614  65.963  1.00 28.50           N
ATOM    135  CA  ALA L  19      23.064 -76.802  65.641  1.00 27.23           C
ATOM    136  C   ALA L  19      23.816 -75.537  66.027  1.00 27.07           C
ATOM    137  O   ALA L  19      23.265 -74.434  65.976  1.00 24.32           O
ATOM    138  CB  ALA L  19      23.231 -77.079  64.157  1.00 30.04           C
ATOM    139  N   SER L  20      25.075 -75.694  66.411  1.00 24.81           N
ATOM    140  CA  SER L  20      25.874 -74.549  66.816  1.00 27.97           C
ATOM    141  C   SER L  20      27.304 -74.675  66.310  1.00 26.96           C
ATOM    142  O   SER L  20      27.912 -75.739  66.408  1.00 27.07           O
ATOM    143  CB  SER L  20      25.863 -74.436  68.344  1.00 29.31           C
ATOM    144  OG  SER L  20      26.505 -73.254  68.779  1.00 39.96           O
"""

def run(pdb_str, expected_ids):
  get_class = iotbx.pdb.common_residue_names_get_class
  mon_lib_srv = mmtbx.monomer_library.server.server()
  rotamer_manager = RotamerEval()
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  result_ids = []
  for residue_group in pdb_hierarchy.residue_groups():
    for conformer in residue_group.conformers():
      for residue in conformer.residues():
        sites_cart = residue.atoms().extract_xyz()
        rotamer_name = rotamer_manager.evaluate_residue(residue=residue)
        print(residue.resname, residue.resseq, rotamer_name)
        result_ids.append(rotamer_name)
        if(get_class(residue.resname) == "common_amino_acid"):
          rotamer_iterator = mon_lib_srv.rotamer_iterator(
              fine_sampling = True,
              comp_id       = residue.resname,
              atom_names    = residue.atoms().extract_name(),
              sites_cart    = sites_cart)
          if(rotamer_iterator is None or
             rotamer_iterator.problem_message is not None or
             rotamer_iterator.rotamer_info is None):
            rotamer_iterator = None
          if(rotamer_iterator is not None):
            d1_min, d2_min = 1.e+9, 1.e+9
            for r, rotamer_sites_cart in rotamer_iterator:
              sites_cart_rot = rotamer_manager.nearest_rotamer_sites_cart(
                residue=residue)
              d1= flex.mean(flex.sqrt((sites_cart - sites_cart_rot).dot()))
              d2= flex.mean(flex.sqrt((sites_cart - rotamer_sites_cart).dot()))
              if(d1 < d1_min):
                d1_min = d1
              if(d2 < d2_min):
                d2_min = d2
            assert approx_equal(d1_min, d2_min)
  assert result_ids == expected_ids

if(__name__ == "__main__"):
  run(pdb_str = pdb_str1, expected_ids=['m-80', 'OUTLIER', 'm100', 'OUTLIER', None, 'm-80', 'p'])
  run(pdb_str = pdb_str2, expected_ids=['Cg_endo'])
  run(pdb_str = pdb_str3, expected_ids=['Cg_endo', 'EXCEPTION', 't'])
