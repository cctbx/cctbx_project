from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
from libtbx import group_args
import mmtbx.restraints
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residues
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager

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

def exercise(rotamer_manager, sin_cos_table,
             d_min = 1.0, resolution_factor = 0.1):
  # Run into a water clash if needed: water is considered as just a map peak.
  #
  # answer PDB
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  # answer map
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_for_map)
  pdb_inp.write_pdb_file(file_name = "for_map.pdb")
  xrs_map = pdb_inp.xray_structure_simple()
  f_calc = xrs_map.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  # poor
  mon_lib_srv = monomer_library.server.server()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = flex.std_string(pdb_poor.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor.pdb")
  #
  target_map_object = group_args(
    data             = target_map,
    f_map_diff       = None,
    miller_array     = f_calc,
    crystal_gridding = fft_map)
  grm = mmtbx.restraints.manager(
    geometry=processed_pdb_file.geometry_restraints_manager(show_energies=False),
    normalization = True)
  sm = mmtbx.refinement.real_space.structure_monitor(
    pdb_hierarchy               = pdb_hierarchy_poor,
    xray_structure              = xrs_poor,
    target_map_object           = target_map_object,
    geometry_restraints_manager = grm.geometry)
  sm.show(prefix="start")
  sm.show_residues(map_cc_all=1.0)
  for i in [1]:
    result = mmtbx.refinement.real_space.fit_residues.manager(
      structure_monitor = sm,
      rotamer_manager   = rotamer_manager,
      sin_cos_table     = sin_cos_table,
      mon_lib_srv       = mon_lib_srv)
    print
    sm.show(prefix="MC %s"%str(i))
    print
    sm.show_residues(map_cc_all=1.0)
  #
  sm.pdb_hierarchy.write_pdb_file(file_name = "refined.pdb")
  dist = xrs_answer.mean_distance(other = sm.xray_structure)
  # assert dist < 0.07, dist # this works on chevy
  assert dist < 0.75, dist # to make it work on marbles

if(__name__ == "__main__"):
  t0 = time.time()
  # load rotamer manager
  rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
  # pre-compute sin and cos tables
  sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  exercise(
    rotamer_manager = rotamer_manager,
    sin_cos_table   = sin_cos_table)
  print "Time: %6.4f"%(time.time()-t0)
