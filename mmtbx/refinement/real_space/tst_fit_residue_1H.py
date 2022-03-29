from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.70 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.70 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.70 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.70 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.70 20.00           N
ATOM      0  HA  ARG A  21       6.805  11.715   8.552  1.00 10.00           H   new
ATOM      0  HB2 ARG A  21       6.647   9.291   8.480  1.00 10.00           H   new
ATOM      0  HB3 ARG A  21       8.315   9.780   8.259  1.00 10.00           H   new
ATOM      0  HG2 ARG A  21       8.764   9.264  10.662  0.70 20.00           H   new
ATOM      0  HG3 ARG A  21       7.104   8.729  10.843  0.70 20.00           H   new
ATOM      0  HD2 ARG A  21       7.476   7.008   9.071  0.70 20.00           H   new
ATOM      0  HD3 ARG A  21       9.141   7.534   8.919  0.70 20.00           H   new
ATOM      0  HE  ARG A  21       9.680   6.499  10.965  0.70 20.00           H   new
ATOM      0 HH11 ARG A  21       6.207   6.415  10.339  0.70 20.00           H   new
ATOM      0 HH12 ARG A  21       5.891   5.227  11.609  0.70 20.00           H   new
ATOM      0 HH21 ARG A  21       9.268   4.978  12.592  0.70 20.00           H   new
ATOM      0 HH22 ARG A  21       7.617   4.418  12.878  0.70 20.00           H   new
TER
HETATM    1  U   U   B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor0 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       6.612   8.656   9.584  0.70 20.00           C
ATOM      7  CD  ARG A  21       6.020   7.767   8.502  0.70 20.00           C
ATOM      8  NE  ARG A  21       4.569   7.657   8.617  0.70 20.00           N
ATOM      9  CZ  ARG A  21       3.771   7.248   7.637  0.70 20.00           C
ATOM     10  NH1 ARG A  21       4.282   6.909   6.460  0.70 20.00           N
ATOM     11  NH2 ARG A  21       2.461   7.180   7.830  0.70 20.00           N
ATOM      0  HA  ARG A  21       6.805  11.715   8.552  1.00 10.00           H   new
ATOM      0  HB2 ARG A  21       7.451   9.672   7.900  1.00 10.00           H   new
ATOM      0  HB3 ARG A  21       8.538   9.528   9.267  1.00 10.00           H   new
ATOM      0  HG2 ARG A  21       7.188   8.048  10.282  0.70 20.00           H   new
ATOM      0  HG3 ARG A  21       5.809   9.122  10.155  0.70 20.00           H   new
ATOM      0  HD2 ARG A  21       6.276   8.169   7.522  0.70 20.00           H   new
ATOM      0  HD3 ARG A  21       6.465   6.774   8.565  0.70 20.00           H   new
ATOM      0  HE  ARG A  21       4.140   7.911   9.507  0.70 20.00           H   new
ATOM      0 HH11 ARG A  21       5.289   6.962   6.307  0.70 20.00           H   new
ATOM      0 HH12 ARG A  21       3.668   6.595   5.709  0.70 20.00           H   new
ATOM      0 HH21 ARG A  21       2.064   7.442   8.732  0.70 20.00           H   new
ATOM      0 HH22 ARG A  21       1.850   6.866   7.076  0.70 20.00           H   new
TER
HETATM   12  U   U   B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor1 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.248  10.924   9.570  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.064  10.337  10.375  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.724   9.472   9.652  1.00 10.00           C
ATOM      6  CG  ARG A  21       8.797   9.112   8.637  1.00 10.00           C
ATOM      7  CD  ARG A  21       9.187   7.647   8.741  1.00 10.00           C
ATOM      8  NE  ARG A  21      10.266   7.301   7.820  1.00 10.00           N
ATOM      9  CZ  ARG A  21      10.871   6.118   7.790  1.00 10.00           C
ATOM     10  NH1 ARG A  21      10.505   5.162   8.634  1.00 10.00           N
ATOM     11  NH2 ARG A  21      11.844   5.891   6.920  1.00 10.00           N
ATOM      0  HA  ARG A  21       6.977  11.151   8.539  1.00 10.00           H   new
ATOM      0  HB2 ARG A  21       8.109   9.283  10.654  1.00 10.00           H   new
ATOM      0  HB3 ARG A  21       6.868   8.812   9.510  1.00 10.00           H   new
ATOM      0  HG2 ARG A  21       8.434   9.323   7.631  1.00 10.00           H   new
ATOM      0  HG3 ARG A  21       9.676   9.737   8.797  1.00 10.00           H   new
ATOM      0  HD2 ARG A  21       9.497   7.426   9.762  1.00 10.00           H   new
ATOM      0  HD3 ARG A  21       8.317   7.025   8.531  1.00 10.00           H   new
ATOM      0  HE  ARG A  21      10.575   8.012   7.158  1.00 10.00           H   new
ATOM      0 HH11 ARG A  21       9.758   5.334   9.307  1.00 10.00           H   new
ATOM      0 HH12 ARG A  21      10.970   4.255   8.610  1.00 10.00           H   new
ATOM      0 HH21 ARG A  21      12.130   6.625   6.272  1.00 10.00           H   new
ATOM      0 HH22 ARG A  21      12.307   4.982   6.898  1.00 10.00           H   new
TER
HETATM   12  U   U   B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor2 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.70 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.70 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.70 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.70 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.70 20.00           N
ATOM      0  HA  ARG A  21       6.805  11.715   8.552  1.00 10.00           H   new
ATOM      0  HB2 ARG A  21       6.647   9.291   8.480  1.00 10.00           H   new
ATOM      0  HB3 ARG A  21       8.315   9.780   8.259  1.00 10.00           H   new
ATOM      0  HG2 ARG A  21       8.764   9.264  10.662  0.70 20.00           H   new
ATOM      0  HG3 ARG A  21       7.104   8.729  10.843  0.70 20.00           H   new
ATOM      0  HD2 ARG A  21       7.476   7.008   9.071  0.70 20.00           H   new
ATOM      0  HD3 ARG A  21       9.141   7.534   8.919  0.70 20.00           H   new
ATOM      0  HE  ARG A  21       9.680   6.499  10.965  0.70 20.00           H   new
ATOM      0 HH11 ARG A  21       6.207   6.415  10.339  0.70 20.00           H   new
ATOM      0 HH12 ARG A  21       5.891   5.227  11.609  0.70 20.00           H   new
ATOM      0 HH21 ARG A  21       9.268   4.978  12.592  0.70 20.00           H   new
ATOM      0 HH22 ARG A  21       7.617   4.418  12.878  0.70 20.00           H   new
TER
HETATM   12  U   U   B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor3 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       8.469   9.666   7.814  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.724   8.229   7.388  0.70 20.00           C
ATOM      8  NE  ARG A  21       9.667   8.147   6.277  0.70 20.00           N
ATOM      9  CZ  ARG A  21       9.331   8.301   5.000  0.70 20.00           C
ATOM     10  NH1 ARG A  21       8.070   8.545   4.668  0.70 20.00           N
ATOM     11  NH2 ARG A  21      10.255   8.208   4.054  0.70 20.00           N
ATOM      0  HA  ARG A  21       6.805  11.715   8.552  1.00 10.00           H   new
ATOM      0  HB2 ARG A  21       7.944   9.204   9.834  1.00 10.00           H   new
ATOM      0  HB3 ARG A  21       6.590   9.203   8.721  1.00 10.00           H   new
ATOM      0  HG2 ARG A  21       8.062  10.230   6.974  0.70 20.00           H   new
ATOM      0  HG3 ARG A  21       9.412  10.139   8.089  0.70 20.00           H   new
ATOM      0  HD2 ARG A  21       9.112   7.663   8.235  0.70 20.00           H   new
ATOM      0  HD3 ARG A  21       7.782   7.764   7.099  0.70 20.00           H   new
ATOM      0  HE  ARG A  21      10.646   7.960   6.494  0.70 20.00           H   new
ATOM      0 HH11 ARG A  21       7.356   8.615   5.393  0.70 20.00           H   new
ATOM      0 HH12 ARG A  21       7.814   8.663   3.688  0.70 20.00           H   new
ATOM      0 HH21 ARG A  21      11.225   8.018   4.306  0.70 20.00           H   new
ATOM      0 HH22 ARG A  21       9.996   8.326   3.075  0.70 20.00           H   new
TER
HETATM   12  U   U   B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

def exercise(pdb_poor_str, i_pdb, d_min = 1.0, resolution_factor = 0.1):
  """
  Fit one residue. There is a huge heavy atom nearby that overlaps with a
  plausible rotamer. With H.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor_str,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["ARG"],
    resolution_factor = resolution_factor)
  #
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = t.vdw,
    crystal_symmetry  = t.crystal_symmetry,
    map_data          = t.target_map,
    backbone_sample   = True,
    rotatable_hd      = t.rotatable_hd,
    rotamer_manager   = t.rotamer_manager,
    sin_cos_table     = t.sin_cos_table,
    mon_lib_srv       = t.mon_lib_srv)
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb))
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer  = t.ph_answer,
    ph_refined = result.pdb_hierarchy,
    tol        = 0.15)

if(__name__ == "__main__"):
  t0 = time.time()
  for i_pdb, pdb_poor_str in enumerate(
                               [pdb_poor0, pdb_poor1, pdb_poor2, pdb_poor3]):
    exercise(
      pdb_poor_str = pdb_poor_str,
      i_pdb        = i_pdb)
  print("Time: %6.4f"%(time.time()-t0))
