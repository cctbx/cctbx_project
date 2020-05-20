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
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.30 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.30 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.30 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.30 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.30 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.30 20.00           N
TER
END
"""

pdb_poor = """\
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
TER
END
"""

pdb_poor_for_map = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.30 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.30 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.30 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.30 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.30 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.30 20.00           N
TER
ATOM      9  O   HOH B  21       8.776  10.791   4.311  1.00  5.00           O
ATOM     10  O   HOH B  22       7.708  11.548   4.090  1.00  5.00           O
ATOM     11  O   HOH B  23       9.698  10.663   3.367  1.00  5.00           O
TER
END
"""

def exercise(d_min = 1.0, resolution_factor = 0.1):
  """
  Fit one residue having weak side chain density. There is a blob nearby that
  overlaps with a plausible rotamer.
  """
 #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = 0,
    d_min             = d_min,
    resolution_factor = resolution_factor,
    residues          = ["ARG"],
    pdb_for_map       = pdb_poor_for_map)
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
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(0))
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer  = t.ph_answer,
    ph_refined = result.pdb_hierarchy,
    tol        = 0.15)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
