from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   23.341   28.568   19.164  90.00  90.00  90.00 P 1
ATOM      1  N   PHE A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  PHE A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   PHE A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   PHE A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  PHE A  58      10.754  19.416  10.061  1.00  8.60           C
ATOM      6  CG  PHE A  58      11.066  20.536   9.101  1.00 10.13           C
ATOM      7  CD1 PHE A  58      11.673  21.701   9.544  1.00  8.55           C
ATOM      9  CE1 PHE A  58      11.950  22.730   8.669  1.00  9.55           C
ATOM     10  CE2 PHE A  58      11.009  21.455   6.871  1.00 10.11           C
ATOM     11  CZ  PHE A  58      11.621  22.610   7.327  1.00 10.71           C
ATOM     12  N   HIS A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     13  CA  HIS A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     14  C   HIS A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     15  O   HIS A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     16  CB  HIS A  59       7.250  14.977   9.971  1.00  7.59           C
ATOM     17  CG  HIS A  59       6.333  15.874  10.738  1.00  9.83           C
ATOM     19  CD2 HIS A  59       5.689  17.009  10.376  1.00  8.78           C
ATOM     20  CE1 HIS A  59       5.211  16.616  12.481  1.00  9.20           C
ATOM     21  NE2 HIS A  59       5.000  17.452  11.479  1.00  9.98           N
ATOM     22  N   TRP A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     23  CA  TRP A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     24  C   TRP A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     25  O   TRP A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     26  CB  TRP A  60      12.197  12.509  10.046  1.00 11.46           C
ATOM     27  CG  TRP A  60      12.933  13.782  10.392  1.00 14.68           C
ATOM     28  CD1 TRP A  60      12.677  14.629  11.446  1.00 12.97           C
ATOM     29  CD2 TRP A  60      13.987  14.394   9.638  1.00 14.91           C
ATOM     30  NE1 TRP A  60      13.524  15.715  11.398  1.00  9.65           N
ATOM     32  CE3 TRP A  60      14.673  14.043   8.472  1.00  8.58           C
ATOM     33  CZ2 TRP A  60      15.350  16.433   9.839  1.00 12.03           C
ATOM     34  CZ3 TRP A  60      15.670  14.879   8.017  1.00 14.50           C
ATOM     35  CH2 TRP A  60      16.002  16.057   8.697  1.00 11.88           C
ATOM     36  N   ARG A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     37  CA  ARG A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     38  C   ARG A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     39  O   ARG A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     40  CB  ARG A  61       9.318   9.570   5.914  1.00 20.00           C
ATOM     41  CG  ARG A  61       9.425   8.639   4.717  1.00 20.00           C
ATOM     42  CD  ARG A  61       8.264   8.745   3.741  1.00 20.00           C
ATOM     43  NE  ARG A  61       7.574  10.028   3.848  1.00 20.00           N
ATOM     44  CZ  ARG A  61       6.414  10.300   3.269  1.00 20.00           C
ATOM     47  N   PRO A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     48  CA  PRO A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     49  C   PRO A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     50  O   PRO A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     51  CB  PRO A  62      14.025   6.257   7.281  1.00  6.15           C
ATOM     52  CG  PRO A  62      13.701   6.929   8.583  1.00  9.40           C
ATOM     53  CD  PRO A  62      12.359   7.621   8.338  1.00  8.83           C
TER      54      PRO A  62
HETATM   55  O   HOH S  30       9.529  12.770  12.418  1.00  9.26           O
HETATM   56  O   HOH S  63      17.561   6.671   6.582  1.00  8.34           O
HETATM   57  O   HOH S 102      16.677   8.520  10.394  1.00  6.42           O
HETATM   58  O   HOH S 116       7.555  10.547  10.276  1.00 11.56           O
HETATM   59  O   HOH S 167       8.683  23.568  14.164  1.00 17.13           O
HETATM   60  O   HOH S 171      15.615  21.403  10.635  1.00  9.12           O
HETATM   61  O   HOH S 176      10.243  21.293  13.219  1.00 18.50           O
HETATM   62  O   HOH S 192       9.980   5.000   7.291  1.00 24.87           O
HETATM   63  O   HOH S 277      18.341  19.685   8.800  1.00 14.61           O
TER      64      HOH S 277
END
"""

pdb_poor = """\
CRYST1   23.341   28.568   19.164  90.00  90.00  90.00 P 1
ATOM      1  N   PHE A  58       8.659  20.073  11.185  1.00  7.73           N
ATOM      2  CA  PHE A  58       9.250  19.144  10.233  1.00  8.65           C
ATOM      3  C   PHE A  58       9.039  17.721  10.706  1.00  9.84           C
ATOM      4  O   PHE A  58       9.023  17.464  11.919  1.00  8.58           O
ATOM      5  CB  PHE A  58      10.754  19.416  10.061  1.00  8.60           C
ATOM      6  CG  PHE A  58      11.545  19.305  11.340  1.00 10.13           C
ATOM      7  CD1 PHE A  58      12.066  18.088  11.748  1.00  8.55           C
ATOM      9  CE1 PHE A  58      12.781  17.983  12.923  1.00  9.55           C
ATOM     10  CE2 PHE A  58      12.466  20.325  13.323  1.00 10.11           C
ATOM     11  CZ  PHE A  58      12.986  19.103  13.714  1.00 10.71           C
ATOM     12  N   HIS A  59       8.887  16.797   9.756  1.00  9.65           N
ATOM     13  CA  HIS A  59       8.720  15.382  10.080  1.00  5.80           C
ATOM     14  C   HIS A  59       9.606  14.539   9.169  1.00 10.35           C
ATOM     15  O   HIS A  59       9.972  14.972   8.075  1.00 10.56           O
ATOM     16  CB  HIS A  59       7.250  14.977   9.971  1.00  7.59           C
ATOM     17  CG  HIS A  59       6.357  15.744  10.892  1.00  9.83           C
ATOM     19  CD2 HIS A  59       6.246  15.728  12.242  1.00  8.78           C
ATOM     20  CE1 HIS A  59       4.758  17.143  11.468  1.00  9.20           C
ATOM     21  NE2 HIS A  59       5.242  16.606  12.575  1.00  9.98           N
ATOM     22  N   TRP A  60       9.964  13.344   9.625  1.00  9.39           N
ATOM     23  CA  TRP A  60      11.067  12.594   9.017  1.00 11.89           C
ATOM     24  C   TRP A  60      10.635  11.189   8.608  1.00  9.81           C
ATOM     25  O   TRP A  60      10.120  10.434   9.430  1.00  8.97           O
ATOM     26  CB  TRP A  60      12.197  12.509  10.046  1.00 11.46           C
ATOM     27  CG  TRP A  60      13.439  11.737   9.668  1.00 14.68           C
ATOM     28  CD1 TRP A  60      13.807  10.487  10.111  1.00 12.97           C
ATOM     29  CD2 TRP A  60      14.513  12.195   8.837  1.00 14.91           C
ATOM     30  NE1 TRP A  60      15.033  10.140   9.587  1.00  9.65           N
ATOM     32  CE3 TRP A  60      14.745  13.369   8.116  1.00  8.58           C
ATOM     33  CZ2 TRP A  60      16.663  11.282   8.064  1.00 12.03           C
ATOM     34  CZ3 TRP A  60      15.913  13.478   7.392  1.00 14.50           C
ATOM     35  CH2 TRP A  60      16.856  12.443   7.368  1.00 11.88           C
ATOM     36  N   ARG A  61      10.858  10.832   7.348  1.00  7.72           N
ATOM     37  CA  ARG A  61      10.510   9.497   6.870  1.00  9.11           C
ATOM     38  C   ARG A  61      11.692   8.812   6.178  1.00 10.61           C
ATOM     39  O   ARG A  61      11.963   9.081   5.006  1.00 11.05           O
ATOM     40  CB  ARG A  61       9.318   9.570   5.914  1.00 20.00           C
ATOM     41  CG  ARG A  61       8.017   9.993   6.576  1.00 20.00           C
ATOM     42  CD  ARG A  61       6.821  10.020   5.637  1.00 20.00           C
ATOM     43  NE  ARG A  61       5.613  10.499   6.303  1.00 20.00           N
ATOM     44  CZ  ARG A  61       4.463  10.727   5.686  1.00 20.00           C
ATOM     47  N   PRO A  62      12.408   7.927   6.895  1.00 10.62           N
ATOM     48  CA  PRO A  62      13.520   7.240   6.220  1.00  7.26           C
ATOM     49  C   PRO A  62      13.019   6.486   5.000  1.00 10.75           C
ATOM     50  O   PRO A  62      11.946   5.890   5.048  1.00 11.44           O
ATOM     51  CB  PRO A  62      14.025   6.257   7.281  1.00  6.15           C
ATOM     52  CG  PRO A  62      12.802   5.932   8.087  1.00  9.40           C
ATOM     53  CD  PRO A  62      12.029   7.250   8.167  1.00  8.83           C
TER      54      PRO A  62
HETATM   55  O   HOH S  30       9.529  12.770  12.418  1.00  9.26           O
HETATM   56  O   HOH S  63      17.561   6.671   6.582  1.00  8.34           O
HETATM   57  O   HOH S 102      16.677   8.520  10.394  1.00  6.42           O
HETATM   58  O   HOH S 116       7.555  10.547  10.276  1.00 11.56           O
HETATM   59  O   HOH S 167       8.683  23.568  14.164  1.00 17.13           O
HETATM   60  O   HOH S 171      15.615  21.403  10.635  1.00  9.12           O
HETATM   61  O   HOH S 176      10.243  21.293  13.219  1.00 18.50           O
HETATM   62  O   HOH S 192       9.980   5.000   7.291  1.00 24.87           O
HETATM   63  O   HOH S 277      18.341  19.685   8.800  1.00 14.61           O
TER      64      HOH S 277
END
"""

pdb_for_map = pdb_answer

def exercise(i_pdb=0, d_min=1.5, resolution_factor=0.1):
  """
  Partial (incomplete) residues. Just run to make sure it does not crash.
  It will not fix incomplete residues.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["ARG","PHE","TRP","PRO","HIS"],
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)
  #
  ph = t.ph_poor
  for i in [1,]:
    result = mmtbx.refinement.real_space.fit_residues.run(
      pdb_hierarchy     = ph,
      vdw_radii         = t.vdw,
      crystal_symmetry  = t.crystal_symmetry,
      map_data          = t.target_map,
      backbone_sample   = True,
      rotatable_hd      = t.rotatable_hd,
      rotamer_manager   = t.rotamer_manager,
      sin_cos_table     = t.sin_cos_table,
      mon_lib_srv       = t.mon_lib_srv)
    ph = result.pdb_hierarchy
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
