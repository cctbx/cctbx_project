from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager

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
ATOM      9  O   HOH B  21       8.776  10.791   4.311  1.00 55.00           O
ATOM     10  O   HOH B  22       7.708  11.548   4.090  1.00 55.00           O
ATOM     11  O   HOH B  23       9.698  10.663   3.367  1.00 55.00           O
TER
END
"""

def exercise(rotamer_manager, sin_cos_table, d_min = 1.0,
             resolution_factor = 0.1):
  # Fit one residue having weak side chain density. There is a blob nearby that
  # overlaps with a plausible rotamer.
  # Making B of HOH smaller will break the test, indicaing potential problem.
  #
  # answer PDB
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  # answer map
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_poor_for_map)
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
  get_class = iotbx.pdb.common_residue_names_get_class
  residue_poor = None
  for model in pdb_hierarchy_poor.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if(get_class(residue.resname) == "common_amino_acid"):
          t0=time.time() # TIMER START
          # refine
          mmtbx.refinement.real_space.fit_residue.run(
            residue         = residue,
            unit_cell       = xrs_poor.unit_cell(),
            target_map      = target_map,
            mon_lib_srv     = mon_lib_srv,
            rotamer_manager = rotamer_manager,
            sin_cos_table   = sin_cos_table)
          sites_cart_poor.set_selected(residue.atoms().extract_i_seq(),
            residue.atoms().extract_xyz())
          print "time (refine): %6.4f" % (time.time()-t0)
  xrs_poor = xrs_poor.replace_sites_cart(sites_cart_poor)
  pdb_hierarchy_poor.adopt_xray_structure(xrs_poor)
  pdb_hierarchy_poor.write_pdb_file(file_name = "refined.pdb")
  dist = xrs_answer.max_distance(other = xrs_poor)
  assert dist < 0.24, dist

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
