from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space

pdb_answer = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
ATOM     17  N   TYR A  20      12.397  11.404   6.631  1.00 20.00           N
ATOM     18  CA  TYR A  20      12.988  10.470   5.681  1.00 20.00           C
ATOM     19  C   TYR A  20      11.909   9.888   4.773  1.00 20.00           C
ATOM     20  O   TYR A  20      11.317  10.600   3.962  1.00 20.00           O
ATOM     21  CB  TYR A  20      13.716   9.345   6.422  1.00 20.00           C
ATOM     22  CG  TYR A  20      12.912   8.484   7.253  1.00 20.00           C
ATOM     23  CD1 TYR A  20      12.307   7.345   6.738  1.00 20.00           C
ATOM     24  CD2 TYR A  20      12.712   8.780   8.596  1.00 20.00           C
ATOM     25  CE1 TYR A  20      11.531   6.525   7.534  1.00 20.00           C
ATOM     26  CE2 TYR A  20      11.937   7.967   9.400  1.00 20.00           C
ATOM     27  CZ  TYR A  20      11.349   6.841   8.863  1.00 20.00           C
ATOM     28  OH  TYR A  20      10.575   6.027   9.660  1.00 20.00           O
TER
ATOM     28  O   HOH B   1      10.575   6.027   9.660  1.00 20.00           O
TER
ATOM     28  O   HOH C   2      15.926   5.414   2.947  1.00 20.00           O
TER
ATOM     28  O   HOH E   4      18.225  11.177   9.392  1.00 20.00           O
TER
ATOM     28  O   HOH F   5      19.097  10.731   7.772  1.00 20.00           O
TER
END
"""


pdb_poor = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
ATOM     17  N   TYR A  20      12.397  11.404   6.631  1.00 20.00           N
ATOM     18  CA  TYR A  20      12.988  10.470   5.681  1.00 20.00           C
ATOM     19  C   TYR A  20      11.909   9.888   4.773  1.00 20.00           C
ATOM     20  O   TYR A  20      11.317  10.600   3.962  1.00 20.00           O
ATOM     21  CB  TYR A  20      13.716   9.345   6.422  1.00 20.00           C
ATOM     22  CG  TYR A  20      13.115   8.984   7.766  1.00 20.00           C
ATOM     23  CD1 TYR A  20      12.032   8.119   7.860  1.00 20.00           C
ATOM     24  CD2 TYR A  20      13.632   9.514   8.942  1.00 20.00           C
ATOM     25  CE1 TYR A  20      11.485   7.789   9.085  1.00 20.00           C
ATOM     26  CE2 TYR A  20      13.092   9.190  10.171  1.00 20.00           C
ATOM     27  CZ  TYR A  20      12.018   8.327  10.236  1.00 20.00           C
ATOM     28  OH  TYR A  20      11.474   8.001  11.459  1.00 20.00           O
TER
ATOM     28  O   HOH B   1      10.575   6.027   9.660  1.00 20.00           O
TER
ATOM     28  O   HOH C   2      15.926   5.414   2.947  1.00 20.00           O
TER
ATOM     28  O   HOH E   4      18.225  11.177   9.392  1.00 20.00           O
TER
ATOM     28  O   HOH F   5      19.097  10.731   7.772  1.00 20.00           O
TER
END
"""

pdb_for_map = """\
CRYST1   16.960   19.455   19.841  90.00  90.00  90.00 P 1
ATOM     17  N   TYR A  20      12.397  11.404   6.631  1.00 20.00           N
ATOM     18  CA  TYR A  20      12.988  10.470   5.681  1.00 20.00           C
ATOM     19  C   TYR A  20      11.909   9.888   4.773  1.00 20.00           C
ATOM     20  O   TYR A  20      11.317  10.600   3.962  1.00 20.00           O
ATOM     21  CB  TYR A  20      13.716   9.345   6.422  1.00 20.00           C
ATOM     22  CG  TYR A  20      12.912   8.484   7.253  1.00 20.00           C
ATOM     23  CD1 TYR A  20      12.307   7.345   6.738  1.00 20.00           C
ATOM     24  CD2 TYR A  20      12.712   8.780   8.596  1.00 20.00           C
ATOM     25  CE1 TYR A  20      11.531   6.525   7.534  1.00 20.00           C
ATOM     26  CE2 TYR A  20      11.937   7.967   9.400  1.00 20.00           C
ATOM     27  CZ  TYR A  20      11.349   6.841   8.863  1.00 20.00           C
ATOM     28  OH  TYR A  20      10.575   6.027   9.660  1.00 20.00           O
TER
ATOM     28  O   HOH C   2      15.926   5.414   2.947  1.00 20.00           O
TER
ATOM     28  O   HOH E   4      18.225  11.177   9.392  1.00 20.00           O
TER
ATOM     28  O   HOH F   5      19.097  10.731   7.772  1.00 20.00           O
TER
END
"""

def exercise(d_min = 1.0, resolution_factor = 0.1, i_pdb = 0):
  """
  Run into a water clash if needed: water is considered as just a map peak.
  """
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["TYR"],
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
      backbone_sample   = True,
      rotatable_hd      = t.rotatable_hd,
      rotamer_manager   = t.rotamer_manager,
      sin_cos_table     = t.sin_cos_table,
      mon_lib_srv       = t.mon_lib_srv)
    ph = result.pdb_hierarchy
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer          = t.ph_answer,
    ph_refined         = result.pdb_hierarchy,
    exclude_atom_names = ["CE1","CE2","CD1","CD2"],
    tol                = 0.35)

if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
