from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1           0
SCALE1      0.064358  0.000000  0.000000        0.00000
SCALE2      0.000000  0.077876  0.000000        0.00000
SCALE3      0.000000  0.000000  0.075792        0.00000
ATOM   2006  N   MSE B  37       8.282   9.046   6.190  1.00 10.00           N
ATOM   2007  CA  MSE B  37       6.863   8.719   6.123  1.00 10.00           C
ATOM   2008  C   MSE B  37       6.057   9.906   5.608  1.00 10.00           C
ATOM   2009  O   MSE B  37       5.000   9.734   5.000  1.00 10.00           O
ATOM   2010  CB  MSE B  37       6.347   8.291   7.498  1.00 10.00           C
ATOM   2011  CG  MSE B  37       7.044   7.066   8.067  1.00 10.00           C
ATOM   2012 SE   MSE B  37       6.355   6.551   9.817  1.00 10.00          Se
ATOM   2013  CE  MSE B  37       7.491   5.000  10.145  1.00 10.00           C
ATOM      0  HA  MSE B  37       6.740   7.889   5.427  1.00 10.00           H
ATOM      0  HB2 MSE B  37       6.468   9.121   8.194  1.00 10.00           H
ATOM      0  HB3 MSE B  37       5.279   8.088   7.426  1.00 10.00           H
ATOM      0  HG2 MSE B  37       6.926   6.232   7.375  1.00 10.00           H
ATOM      0  HG3 MSE B  37       8.113   7.265   8.146  1.00 10.00           H
ATOM      0  HE1 MSE B  37       7.241   4.564  11.112  1.00 10.00           H
ATOM      0  HE2 MSE B  37       7.331   4.259   9.361  1.00 10.00           H
ATOM      0  HE3 MSE B  37       8.537   5.308  10.145  1.00 10.00           H
TER
"""

pdb_poor0 = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1
ATOM   2006  N   MSE B  37       8.282   9.046   6.190  1.00 10.00           N
ATOM   2007  CA  MSE B  37       6.863   8.719   6.123  1.00 10.00           C
ATOM   2008  C   MSE B  37       6.057   9.906   5.608  1.00 10.00           C
ATOM   2009  O   MSE B  37       5.000   9.734   5.000  1.00 10.00           O
ATOM   2010  CB  MSE B  37       6.347   8.291   7.498  1.00 10.00           C
ATOM   2011  CG  MSE B  37       4.889   7.862   7.506  1.00 10.00           C
ATOM   2012 SE   MSE B  37       4.305   7.207   9.247  1.00 10.00          Se
ATOM   2013  CE  MSE B  37       5.093   5.424   9.180  1.00 10.00           C
ATOM      9  HA  MSE B  37       6.740   7.889   5.427  1.00 10.00           H
ATOM     10  HB2 MSE B  37       6.960   7.467   7.864  1.00 10.00           H
ATOM     11  HB3 MSE B  37       6.475   9.118   8.196  1.00 10.00           H
ATOM     12  HG2 MSE B  37       4.265   8.705   7.212  1.00 10.00           H
ATOM     13  HG3 MSE B  37       4.740   7.081   6.760  1.00 10.00           H
ATOM     14  HE1 MSE B  37       4.863   4.890  10.101  1.00 10.00           H
ATOM     15  HE2 MSE B  37       4.682   4.876   8.331  1.00 10.00           H
ATOM     16  HE3 MSE B  37       6.174   5.507   9.069  1.00 10.00           H
TER
"""

def exercise(pdb_poor_str, i_pdb, d_min = 1.0, resolution_factor = 0.1):
  """
  Fit one non-standard residue: MSE with H.
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
    tol        = 0.5)

if(__name__ == "__main__"):
  t0 = time.time()
  for i_pdb, pdb_poor_str in enumerate([pdb_poor0,]):
    exercise(
      pdb_poor_str = pdb_poor_str,
      i_pdb        = i_pdb)
  print("Time: %6.4f"%(time.time()-t0))
