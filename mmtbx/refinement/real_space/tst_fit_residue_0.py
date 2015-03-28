from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
import math
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager

pdb_answer = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1
ATOM    369  N   MSE A  67       6.997   5.582   8.196  1.00  9.89           N
ATOM    370  CA  MSE A  67       6.979   5.904   6.774  1.00 11.78           C
ATOM    371  CB  MSE A  67       7.987   5.037   6.015  1.00 20.39           C
ATOM    372  CG  MSE A  67       9.434   5.204   6.465  1.00 24.20           C
ATOM    373 SE   MSE A  67      10.168   6.982   6.124  1.00 45.90          Se
ATOM    374  CE  MSE A  67       9.820   7.813   7.855  1.00 13.67           C
ATOM    375  C   MSE A  67       5.584   5.722   6.185  1.00  9.34           C
ATOM    376  O   MSE A  67       4.995   6.665   5.656  1.00 10.34           O
TER
"""

pdb_poor0 = """\
CRYST1   15.538   12.841   13.194  90.00  90.00  90.00 P 1
ATOM    369  N   MSE A  67       6.949   5.583   8.194  1.00  9.89           N
ATOM    370  CA  MSE A  67       6.977   5.887   6.768  1.00 11.78           C
ATOM    371  CB  MSE A  67       7.997   5.000   6.052  1.00 20.39           C
ATOM    372  CG  MSE A  67       8.820   5.724   5.000  1.00 24.20           C
ATOM    373 SE   MSE A  67      10.538   6.343   5.685  1.00 45.90          Se
ATOM    374  CE  MSE A  67       9.911   7.841   6.765  1.00 13.67           C
ATOM    375  C   MSE A  67       5.598   5.712   6.141  1.00  9.34           C
ATOM    376  O   MSE A  67       5.000   6.673   5.657  1.00 10.34           O
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
  #
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy_poor.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if(get_class(residue.resname) == "common_amino_acid"):
          # negate map
          t0=time.time() # TIMER START
          negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
            xray_structure          = xrs_poor,
            selection_within_radius = 5,
            iselection              = residue.atoms().extract_i_seq())
          target_map_ = mmtbx.refinement.real_space.\
            negate_map_around_selected_atoms_except_selected_atoms(
              xray_structure   = xrs_poor,
              map_data         = target_map,
              negate_selection = negate_selection,
              atom_radius      = 4)
          print "  time (negate map): %6.4f" % (time.time()-t0)
          # refine
          mmtbx.refinement.real_space.fit_residue.run(
            residue         = residue,
            unit_cell       = xrs_poor.unit_cell(),
            target_map      = target_map_,
            mon_lib_srv     = mon_lib_srv,
            rotamer_manager = rotamer_manager,
            sin_cos_table   = sin_cos_table)
          sites_cart_poor.set_selected(residue.atoms().extract_i_seq(),
            residue.atoms().extract_xyz())
          print "  time (refine): %6.4f" % (time.time()-t0)
  xrs_poor = xrs_poor.replace_sites_cart(sites_cart_poor)
  pdb_hierarchy_poor.adopt_xray_structure(xrs_poor)
  pdb_hierarchy_poor.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb))
  dist = xrs_answer.max_distance(other = xrs_poor)
  assert dist < 0.20,  dist

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
  print "Time: %6.4f"%(time.time()-t0)
