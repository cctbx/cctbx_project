from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1
ATOM    369  N   MSE A  67       6.949   5.583   8.194  1.00  9.89           N
ATOM    370  CA  MSE A  67       6.974   5.900   6.771  1.00 11.78           C
ATOM    375  C   MSE A  67       5.598   5.712   6.141  1.00  9.34           C
ATOM    376  O   MSE A  67       4.987   6.669   5.665  1.00 10.34           O
ATOM    371  CB  MSE A  67       8.007   5.034   6.048  1.00 20.39           C
ATOM    372  CG  MSE A  67       9.434   5.367   6.373  1.00 24.20           C
ATOM    373 SE   MSE A  67       9.847   7.212   6.188  1.00 45.90          SE
ATOM    374  CE  MSE A  67      11.413   7.254   7.196  1.00 13.67           C
TER
"""

pdb_poor0 = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1
ATOM    369  N   MSE A  67       6.949   5.583   8.194  1.00  9.89           N
ATOM    375  C   MSE A  67       5.598   5.712   6.141  1.00  9.34           C
ATOM    376  O   MSE A  67       5.000   6.673   5.657  1.00 10.34           O
ATOM    370  CA  MSE A  67       6.977   5.887   6.768  1.00 11.78           C
ATOM    371  CB  MSE A  67       7.997   5.000   6.052  1.00 20.39           C
ATOM    372  CG  MSE A  67       8.820   5.724   5.000  1.00 24.20           C
ATOM    373 SE   MSE A  67      10.538   6.343   5.685  1.00 45.90          Se
ATOM    374  CE  MSE A  67       9.911   7.841   6.765  1.00 13.67           C
TER
"""

def exercise(pdb_poor_str, i_pdb, d_min = 1.0, resolution_factor = 0.1):
  """
  Fit one non-standard residue: MSE.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor_str,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["MSE"],
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
    tol        = 0.44)

if(__name__ == "__main__"):
  t0 = time.time()
  for i_pdb, pdb_poor_str in enumerate([pdb_poor0,]):
    exercise(
      pdb_poor_str = pdb_poor_str,
      i_pdb        = i_pdb)
  print("Time: %6.4f"%(time.time()-t0))
