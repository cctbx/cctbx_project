from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.171  11.092   9.451  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.017  10.416  10.266  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.564   9.645   9.143  1.00 10.00           C
ATOM      6  CG  ARG A  21       8.560   9.501   8.003  1.00 10.00           C
ATOM      7  CD  ARG A  21       8.125  10.300   6.785  1.00 10.00           C
ATOM      8  NE  ARG A  21       8.926   9.982   5.607  1.00 10.00           N
ATOM      9  CZ  ARG A  21       8.893  10.674   4.473  1.00 10.00           C
ATOM     10  NH1 ARG A  21       8.093  11.727   4.359  1.00 10.00           N
ATOM     11  NH2 ARG A  21       9.655  10.313   3.451  1.00 10.00           N
END
"""

pdb_for_map = """
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  O   HOH S   0       8.318  11.834   9.960  1.00 10.00           O
ATOM      2  O   HOH S   1       7.171  11.092   9.451  1.00 10.00           O
ATOM      3  O   HOH S   2       6.012  11.120  10.440  1.00 10.00           O
ATOM      4  O   HOH S   3       5.017  10.416  10.266  1.00 10.00           O
ATOM      5  O   HOH S   4       7.564   9.645   9.143  1.00 10.00           O
ATOM      6  O   HOH S   5       8.560   9.501   8.003  1.00 10.00           O
ATOM      7  O   HOH S   6       8.125  10.300   6.785  1.00 10.00           O
ATOM      8  O   HOH S   7       8.926   9.982   5.607  1.00 10.00           O
ATOM      9  O   HOH S   8       8.893  10.674   4.473  1.00 10.00           O
ATOM     10  O   HOH S   9       8.093  11.727   4.359  1.00 10.00           O
ATOM     11  O   HOH S  10       9.655  10.313   3.451  1.00 10.00           O
TER
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
ATOM      1  O   HOH S   0       8.318  11.834   9.960  1.00 10.00           O
ATOM      2  O   HOH S   1       7.171  11.092   9.451  1.00 10.00           O
ATOM      3  O   HOH S   2       6.012  11.120  10.440  1.00 10.00           O
ATOM      4  O   HOH S   3       5.017  10.416  10.266  1.00 10.00           O
ATOM      5  O   HOH S   4       7.564   9.645   9.143  1.00 10.00           O
ATOM      6  O   HOH S   5       8.560   9.501   8.003  1.00 10.00           O
ATOM      7  O   HOH S   6       8.125  10.300   6.785  1.00 10.00           O
ATOM      8  O   HOH S   7       8.926   9.982   5.607  1.00 10.00           O
ATOM      9  O   HOH S   8       8.893  10.674   4.473  1.00 10.00           O
ATOM     10  O   HOH S   9       8.093  11.727   4.359  1.00 10.00           O
ATOM     11  O   HOH S  10       9.655  10.313   3.451  1.00 10.00           O
TER
END
"""

def exercise(d_min = 1.0, resolution_factor = 0.1):
  """
  Make sure it kicks off existing water. Simple case: no alternatives.
  """
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = 0,
    d_min             = d_min,
    resolution_factor = resolution_factor,
    residues          = ["ARG"],
    pdb_for_map       = pdb_for_map)
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
  sel = result.pdb_hierarchy.atom_selection_cache().selection("not water")
  result_hierarchy = result.pdb_hierarchy.select(sel)
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer  = t.ph_answer,
    ph_refined = result_hierarchy,
    tol        = 0.4)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
