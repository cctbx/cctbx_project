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
CRYST1   16.402   13.234   20.833  90.00  90.00  90.00 P 1
ATOM   2685  N   PHE A 364       3.630   3.286   4.228  1.00  9.46           N
ATOM   2686  CA  PHE A 364       4.866   4.027   4.477  1.00 10.15           C
ATOM   2687  CB  PHE A 364       4.987   4.446   5.947  1.00  9.57           C
ATOM   2688  CG  PHE A 364       6.322   5.060   6.261  1.00 10.71           C
ATOM   2689  CD1 PHE A 364       6.708   6.251   5.660  1.00 11.42           C
ATOM   2690  CE1 PHE A 364       7.965   6.803   5.906  1.00 12.47           C
ATOM   2691  CZ  PHE A 364       8.844   6.165   6.756  1.00 12.01           C
ATOM   2692  CE2 PHE A 364       8.485   4.967   7.348  1.00 12.45           C
ATOM   2693  CD2 PHE A 364       7.223   4.410   7.089  1.00  9.88           C
ATOM   2694  C   PHE A 364       6.107   3.219   4.071  1.00 11.38           C
ATOM   2695  O   PHE A 364       6.958   3.695   3.331  1.00 11.71           O
ATOM   2983  CA  THR A 401       5.525   3.293  16.950  1.00  7.97           C
ATOM   2987  C   THR A 401       6.000   3.929  15.643  1.00  9.81           C
ATOM   2988  O   THR A 401       6.802   3.348  14.900  1.00  9.23           O
ATOM   2989  N   ALA A 402       5.502   5.133  15.362  1.00  8.84           N
ATOM   2990  CA  ALA A 402       5.936   5.844  14.166  1.00  8.49           C
ATOM   2991  CB  ALA A 402       5.174   7.167  14.035  1.00  7.75           C
ATOM   2992  C   ALA A 402       7.443   6.099  14.165  1.00  8.21           C
ATOM   2993  O   ALA A 402       8.116   5.862  13.150  1.00  8.01           O
ATOM   2994  N   PHE A 403       7.977   6.586  15.283  1.00  9.12           N
ATOM   2995  CA  PHE A 403       9.405   6.873  15.345  1.00  8.30           C
ATOM   2996  CB  PHE A 403       9.775   7.458  16.711  1.00  9.28           C
ATOM   2997  CG  PHE A 403      11.202   7.901  16.805  1.00 12.22           C
ATOM   2998  CD1 PHE A 403      11.566   9.176  16.412  1.00 12.95           C
ATOM   2999  CE1 PHE A 403      12.887   9.593  16.493  1.00 14.50           C
ATOM   3002  CD2 PHE A 403      12.176   7.050  17.291  1.00 14.12           C
ATOM   3003  C   PHE A 403      10.246   5.627  15.062  1.00  9.33           C
ATOM   3004  O   PHE A 403      11.192   5.667  14.279  1.00 11.43           O
ATOM   3027  N   PHE A 407       9.045   6.777  10.558  1.00 10.36           N
ATOM   3028  CA  PHE A 407       7.962   7.610  10.037  1.00 10.78           C
ATOM   3029  CB  PHE A 407       6.620   6.885  10.202  1.00  7.96           C
ATOM   3030  CG  PHE A 407       5.563   7.251   9.194  1.00  9.05           C
ATOM   3031  CD1 PHE A 407       4.361   6.552   9.194  1.00  8.53           C
ATOM   3032  CE1 PHE A 407       3.369   6.832   8.278  1.00 10.10           C
ATOM   3033  CZ  PHE A 407       3.556   7.820   7.348  1.00  9.45           C
ATOM   3034  CE2 PHE A 407       4.742   8.549   7.336  1.00  9.34           C
ATOM   3035  CD2 PHE A 407       5.743   8.259   8.260  1.00  9.48           C
ATOM   3036  C   PHE A 407       7.986   8.907  10.843  1.00  8.79           C
ATOM   3037  O   PHE A 407       7.084   9.187  11.629  1.00 10.53           O
TER
"""


pdb_poor = pdb_answer

pdb_for_map_1 = pdb_answer

pdb_for_map_2 = """\
CRYST1   16.402   13.234   20.833  90.00  90.00  90.00 P 1
ATOM   2685  N   PHE A 364       3.630   3.286   4.228  1.00  9.46           N
ATOM   2686  CA  PHE A 364       4.866   4.027   4.477  1.00 10.15           C
ATOM   2687  CB  PHE A 364       4.987   4.446   5.947  1.00  9.57           C
ATOM   2688  CG  PHE A 364       6.322   5.060   6.261  1.00 10.71           C
ATOM   2689  CD1 PHE A 364       6.708   6.251   5.660  1.00 11.42           C
ATOM   2690  CE1 PHE A 364       7.965   6.803   5.906  1.00 12.47           C
ATOM   2691  CZ  PHE A 364       8.844   6.165   6.756  1.00 12.01           C
ATOM   2692  CE2 PHE A 364       8.485   4.967   7.348  1.00 12.45           C
ATOM   2693  CD2 PHE A 364       7.223   4.410   7.089  1.00  9.88           C
ATOM   2694  C   PHE A 364       6.107   3.219   4.071  1.00 11.38           C
ATOM   2695  O   PHE A 364       6.958   3.695   3.331  1.00 11.71           O
ATOM   2983  CA  THR A 401       5.525   3.293  16.950  1.00  7.97           C
ATOM   2987  C   THR A 401       6.000   3.929  15.643  1.00  9.81           C
ATOM   2988  O   THR A 401       6.802   3.348  14.900  1.00  9.23           O
ATOM   2989  N   ALA A 402       5.502   5.133  15.362  1.00  8.84           N
ATOM   2990  CA  ALA A 402       5.936   5.844  14.166  1.00  8.49           C
ATOM   2991  CB  ALA A 402       5.174   7.167  14.035  1.00  7.75           C
ATOM   2992  C   ALA A 402       7.443   6.099  14.165  1.00  8.21           C
ATOM   2993  O   ALA A 402       8.116   5.862  13.150  1.00  8.01           O
ATOM   2994  N   PHE A 403       7.977   6.586  15.283  1.00  9.12           N
ATOM   2995  CA  PHE A 403       9.405   6.873  15.345  1.00  8.30           C
ATOM   2996  CB  PHE A 403       9.775   7.458  16.711  1.00  9.28           C
ATOM   2997  CG  PHE A 403      11.202   7.901  16.805  1.00 12.22           C
ATOM   2998  CD1 PHE A 403      11.566   9.176  16.412  1.00 12.95           C
ATOM   2999  CE1 PHE A 403      12.887   9.593  16.493  1.00 14.50           C
ATOM   3002  CD2 PHE A 403      12.176   7.050  17.291  1.00 14.12           C
ATOM   3003  C   PHE A 403      10.246   5.627  15.062  1.00  9.33           C
ATOM   3004  O   PHE A 403      11.192   5.667  14.279  1.00 11.43           O
ATOM   3027  N   PHE A 407       9.045   6.777  10.558  1.00 10.36           N
ATOM   3028  CA  PHE A 407       7.962   7.610  10.037  1.00 10.78           C
ATOM   3029  CB  PHE A 407       6.620   6.885  10.202  1.00  7.96           C
ATOM   3030  CG  PHE A 407       5.563   7.251   9.194  1.00  9.05           C
ATOM   3031  CD1 PHE A 407       4.361   6.552   9.194  1.00  8.53           C
ATOM   3032  CE1 PHE A 407       3.369   6.832   8.278  0.00 10.10           C
ATOM   3033  CZ  PHE A 407       3.556   7.820   7.348  0.00  9.45           C
ATOM   3034  CE2 PHE A 407       4.742   8.549   7.336  0.00  9.34           C
ATOM   3035  CD2 PHE A 407       5.743   8.259   8.260  1.00  9.48           C
ATOM   3036  C   PHE A 407       7.986   8.907  10.843  1.00  8.79           C
ATOM   3037  O   PHE A 407       7.084   9.187  11.629  1.00 10.53           O
TER
"""


def exercise(i_pdb, pdb_for_map, rotamer_manager, sin_cos_table,
             d_min = 1.5, resolution_factor = 0.1):
  # Best fitting residue is a rotamer outlier (PHE 407), two scenarious:
  #   - outlier fits density perfectly
  #   - outlier fits not so good.
  # No better options to fit other than keep the outlier unchanged.
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
  mtz_object.write(file_name = "answer_%s.mtz"%str(i_pdb))
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
  result = mmtbx.refinement.real_space.fit_residues.manager(
    structure_monitor = sm,
    rotamer_manager   = rotamer_manager,
    sin_cos_table     = sin_cos_table,
    mon_lib_srv       = mon_lib_srv)
  #
  sm.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb))
  dist = xrs_answer.mean_distance(other = sm.xray_structure)
  assert dist < 1.e-6, dist

if(__name__ == "__main__"):
  t0 = time.time()
  # load rotamer manager
  rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
  # pre-compute sin and cos tables
  sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  for i_pdb, pdb_for_map in enumerate([pdb_for_map_1, pdb_for_map_2]):
    exercise(
      i_pdb           = i_pdb,
      pdb_for_map     = pdb_for_map,
      rotamer_manager = rotamer_manager,
      sin_cos_table   = sin_cos_table)
  print "Time: %6.4f"%(time.time()-t0)
