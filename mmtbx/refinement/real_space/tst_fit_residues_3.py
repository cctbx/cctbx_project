from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   16.402   13.234   20.833  90.00  90.00  90.00 P 1
ATOM      1  N   PHE A 364       3.630   3.286   4.228  1.00  9.46           N
ATOM      2  CA  PHE A 364       4.866   4.027   4.477  1.00 10.15           C
ATOM      3  C   PHE A 364       6.107   3.219   4.071  1.00 11.38           C
ATOM      4  O   PHE A 364       6.958   3.695   3.331  1.00 11.71           O
ATOM      5  CB  PHE A 364       4.987   4.446   5.947  1.00  9.57           C
ATOM      6  CG  PHE A 364       6.322   5.060   6.261  1.00 10.71           C
ATOM      7  CD1 PHE A 364       6.708   6.251   5.660  1.00 11.42           C
ATOM      8  CD2 PHE A 364       7.223   4.410   7.089  1.00  9.88           C
ATOM      9  CE1 PHE A 364       7.965   6.803   5.906  1.00 12.47           C
ATOM     10  CE2 PHE A 364       8.485   4.967   7.348  1.00 12.45           C
ATOM     11  CZ  PHE A 364       8.844   6.165   6.756  1.00 12.01           C
ATOM     12  CA  THR A 401       5.525   3.293  16.950  1.00  7.97           C
ATOM     13  C   THR A 401       6.000   3.929  15.643  1.00  9.81           C
ATOM     14  O   THR A 401       6.802   3.348  14.900  1.00  9.23           O
ATOM     15  N   ALA A 402       5.502   5.133  15.362  1.00  8.84           N
ATOM     16  CA  ALA A 402       5.936   5.844  14.166  1.00  8.49           C
ATOM     17  C   ALA A 402       7.443   6.099  14.165  1.00  8.21           C
ATOM     18  O   ALA A 402       8.116   5.862  13.150  1.00  8.01           O
ATOM     19  CB  ALA A 402       5.174   7.167  14.035  1.00  7.75           C
ATOM     20  N   PHE A 403       7.977   6.586  15.283  1.00  9.12           N
ATOM     21  CA  PHE A 403       9.405   6.873  15.345  1.00  8.30           C
ATOM     22  C   PHE A 403      10.246   5.627  15.062  1.00  9.33           C
ATOM     23  O   PHE A 403      11.192   5.667  14.279  1.00 11.43           O
ATOM     24  CB  PHE A 403       9.775   7.458  16.711  1.00  9.28           C
ATOM     25  CG  PHE A 403      11.202   7.901  16.805  1.00 12.22           C
ATOM     26  CD1 PHE A 403      11.566   9.176  16.412  1.00 12.95           C
ATOM     27  CD2 PHE A 403      12.176   7.050  17.291  1.00 14.12           C
ATOM     28  CE1 PHE A 403      12.887   9.593  16.493  1.00 14.50           C
ATOM     29  N   PHE A 407       9.045   6.777  10.558  1.00 10.36           N
ATOM     30  CA  PHE A 407       7.962   7.610  10.037  1.00 10.78           C
ATOM     31  C   PHE A 407       7.986   8.907  10.843  1.00  8.79           C
ATOM     32  O   PHE A 407       7.084   9.187  11.629  1.00 10.53           O
ATOM     33  CB  PHE A 407       6.620   6.885  10.202  1.00  7.96           C
ATOM     34  CG  PHE A 407       5.960   6.454   8.919  1.00  9.05           C
ATOM     35  CD1 PHE A 407       6.506   5.396   8.202  1.00  8.53           C
ATOM     36  CD2 PHE A 407       4.826   7.076   8.421  1.00  9.48           C
ATOM     37  CE1 PHE A 407       5.954   4.976   7.010  1.00 10.10           C
ATOM     38  CE2 PHE A 407       4.261   6.653   7.220  1.00  9.34           C
ATOM     39  CZ  PHE A 407       4.841   5.602   6.515  1.00  9.45           C
TER
"""


pdb_poor = pdb_answer

pdb_for_map_1 = pdb_answer

pdb_for_map_2 = """\
CRYST1   16.402   13.234   20.833  90.00  90.00  90.00 P 1
ATOM      1  N   PHE A 364       3.630   3.286   4.228  1.00  9.46           N
ATOM      2  CA  PHE A 364       4.866   4.027   4.477  1.00 10.15           C
ATOM      3  C   PHE A 364       6.107   3.219   4.071  1.00 11.38           C
ATOM      4  O   PHE A 364       6.958   3.695   3.331  1.00 11.71           O
ATOM      5  CB  PHE A 364       4.987   4.446   5.947  1.00  9.57           C
ATOM      6  CG  PHE A 364       6.322   5.060   6.261  1.00 10.71           C
ATOM      7  CD1 PHE A 364       6.708   6.251   5.660  1.00 11.42           C
ATOM      8  CD2 PHE A 364       7.223   4.410   7.089  1.00  9.88           C
ATOM      9  CE1 PHE A 364       7.965   6.803   5.906  1.00 12.47           C
ATOM     10  CE2 PHE A 364       8.485   4.967   7.348  1.00 12.45           C
ATOM     11  CZ  PHE A 364       8.844   6.165   6.756  1.00 12.01           C
ATOM     12  CA  THR A 401       5.525   3.293  16.950  1.00  7.97           C
ATOM     13  C   THR A 401       6.000   3.929  15.643  1.00  9.81           C
ATOM     14  O   THR A 401       6.802   3.348  14.900  1.00  9.23           O
ATOM     15  N   ALA A 402       5.502   5.133  15.362  1.00  8.84           N
ATOM     16  CA  ALA A 402       5.936   5.844  14.166  1.00  8.49           C
ATOM     17  C   ALA A 402       7.443   6.099  14.165  1.00  8.21           C
ATOM     18  O   ALA A 402       8.116   5.862  13.150  1.00  8.01           O
ATOM     19  CB  ALA A 402       5.174   7.167  14.035  1.00  7.75           C
ATOM     20  N   PHE A 403       7.977   6.586  15.283  1.00  9.12           N
ATOM     21  CA  PHE A 403       9.405   6.873  15.345  1.00  8.30           C
ATOM     22  C   PHE A 403      10.246   5.627  15.062  1.00  9.33           C
ATOM     23  O   PHE A 403      11.192   5.667  14.279  1.00 11.43           O
ATOM     24  CB  PHE A 403       9.775   7.458  16.711  1.00  9.28           C
ATOM     25  CG  PHE A 403      11.202   7.901  16.805  1.00 12.22           C
ATOM     26  CD1 PHE A 403      11.566   9.176  16.412  1.00 12.95           C
ATOM     27  CD2 PHE A 403      12.176   7.050  17.291  1.00 14.12           C
ATOM     28  CE1 PHE A 403      12.887   9.593  16.493  1.00 14.50           C
ATOM     29  N   PHE A 407       9.045   6.777  10.558  1.00 10.36           N
ATOM     30  CA  PHE A 407       7.962   7.610  10.037  1.00 10.78           C
ATOM     31  C   PHE A 407       7.986   8.907  10.843  1.00  8.79           C
ATOM     32  O   PHE A 407       7.084   9.187  11.629  1.00 10.53           O
ATOM     33  CB  PHE A 407       6.620   6.885  10.202  1.00  7.96           C
ATOM     34  CG  PHE A 407       5.960   6.454   8.919  1.00  9.05           C
ATOM     35  CD1 PHE A 407       6.506   5.396   8.202  1.00  8.53           C
ATOM     36  CD2 PHE A 407       4.826   7.076   8.421  1.00  9.48           C
ATOM     37  CE1 PHE A 407       5.954   4.976   7.010  1.00 10.10           C
ATOM     38  CE2 PHE A 407       4.261   6.653   7.220  1.00  9.34           C
ATOM     39  CZ  PHE A 407       4.841   5.602   6.515  1.00  9.45           C
TER
"""


def exercise(i_pdb, pdb_for_map, d_min = 1.5, resolution_factor = 0.1):
  """
  Best fitting residue is a rotamer outlier (PHE 407), two scenarious:
    - outlier fits density perfectly
    - outlier fits not so good.
  No better options to fit other than keep the outlier unchanged.
  """
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["PHE","THR","ALA"],
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)
  #
  ph = t.ph_poor
  for i in [1,2,3]:
    result = mmtbx.refinement.real_space.fit_residues.run(
      pdb_hierarchy     = ph,
      vdw_radii         = t.vdw,
      crystal_symmetry  = t.crystal_symmetry,
      map_data          = t.target_map,
      backbone_sample   = False,
      rotatable_hd      = t.rotatable_hd,
      rotamer_manager   = t.rotamer_manager,
      sin_cos_table     = t.sin_cos_table,
      mon_lib_srv       = t.mon_lib_srv)
    ph = result.pdb_hierarchy
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer  = t.ph_answer,
    ph_refined = result.pdb_hierarchy,
    tol        = 0.16)

if(__name__ == "__main__"):
  t0 = time.time()
  for i_pdb, pdb_for_map in enumerate([pdb_for_map_1, pdb_for_map_2]):
    exercise(
      i_pdb       = i_pdb,
      pdb_for_map = pdb_for_map)
  print("Time: %6.4f"%(time.time()-t0))
