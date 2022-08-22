from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO H 308       8.327   5.395   8.395  1.00 13.20           N
ATOM      2  CA  PRO H 308       9.300   4.740   7.514  1.00  8.82           C
ATOM      3  C   PRO H 308       8.901   4.816   6.043  1.00  9.64           C
ATOM      4  O   PRO H 308       7.821   4.340   5.693  1.00 11.48           O
ATOM      5  CB  PRO H 308       9.296   3.291   8.005  1.00  6.98           C
ATOM      6  CG  PRO H 308       8.840   3.381   9.418  1.00 10.47           C
ATOM      7  CD  PRO H 308       7.827   4.490   9.456  1.00  9.83           C
TER
"""

pdb_poor = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO H 308       8.327   5.395   8.395  1.00 13.20           N
ATOM      2  CA  PRO H 308       9.300   4.740   7.514  1.00  8.82           C
ATOM      3  C   PRO H 308       8.901   4.816   6.043  1.00  9.64           C
ATOM      4  O   PRO H 308       7.821   4.340   5.693  1.00 11.48           O
ATOM      5  CB  PRO H 308       9.296   3.291   8.005  1.00  6.98           C
ATOM      6  CG  PRO H 308       7.942   3.104   8.592  1.00 10.47           C
ATOM      7  CD  PRO H 308       7.572   4.422   9.209  1.00  9.83           C
TER
"""

def exercise(pdb_poor_str, i_pdb, d_min = 1.0, resolution_factor = 0.1):
  """
  Fit PRO. PRO is a special case. No H.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor_str,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["PRO"],
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
    tol        = 0.002)

def exercise_2(pdb_poor_str, i_pdb, d_min = 1.0, resolution_factor = 0.1):
  """
  Fit PRO. PRO is a special case. No H.
  Same as above, but without a map.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor_str,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["PRO"],
    resolution_factor = resolution_factor)
  #
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = t.vdw,
    crystal_symmetry  = t.crystal_symmetry,
    # map_data          = t.target_map,
    # backbone_sample   = True,
    rotatable_hd      = t.rotatable_hd,
    rotamer_manager   = t.rotamer_manager,
    sin_cos_table     = t.sin_cos_table,
    mon_lib_srv       = t.mon_lib_srv)
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb))

if(__name__ == "__main__"):
  t0 = time.time()
  exercise(
    pdb_poor_str = pdb_poor,
    i_pdb        = 0)
  exercise_2(
    pdb_poor_str = pdb_poor,
    i_pdb        = 1)
  print("Time: %6.4f"%(time.time()-t0))
