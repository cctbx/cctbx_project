from __future__ import division
from __future__ import print_function
from builtins import str
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residues
import iotbx.pdb
import math
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager

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


def exercise(pdb_poor_str, rotamer_manager, sin_cos_table, i_pdb, d_min = 1.0,
             resolution_factor = 0.1):
  # Fit one non-standard residue: MSE .
  #
  # answer
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  pdb_inp.write_pdb_file(file_name = "answer_%s.pdb"%str(i_pdb))
  xrs_answer = pdb_inp.xray_structure_simple()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
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
    raw_records              = flex.std_string(pdb_poor_str.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor_%s.pdb"%str(i_pdb))
  mon_lib_srv = monomer_library.server.server()
  #
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = pdb_hierarchy_poor,
    crystal_symmetry  = xrs_poor.crystal_symmetry(),
    map_data          = target_map,
    do_all            = True,
    rotamer_manager   = rotamer_manager,
    sin_cos_table     = sin_cos_table,
    mon_lib_srv       = mon_lib_srv)
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb))
  dist = flex.max(flex.sqrt((xrs_answer.sites_cart() -
    result.pdb_hierarchy.atoms().extract_xyz()).dot()))
  assert dist < 0.44,  dist

if(__name__ == "__main__"):
  t0 = time.time()
  # load rotamer manager
  rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
  # pre-compute sin and cos tables
  sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  for i_pdb, pdb_poor_str in enumerate([pdb_poor0,]):
    exercise(
      pdb_poor_str    = pdb_poor_str,
      rotamer_manager = rotamer_manager,
      sin_cos_table   = sin_cos_table,
      i_pdb           = i_pdb)
  print("Time: %6.4f"%(time.time()-t0))
