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
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.70 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.70 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.70 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.70 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.70 20.00           N
TER
HETATM    1  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER
END
"""

pdb_poor0 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       6.612   8.656   9.584  0.70 20.00           C
ATOM      7  CD  ARG A  21       6.020   7.767   8.502  0.70 20.00           C
ATOM      8  NE  ARG A  21       4.569   7.657   8.617  0.70 20.00           N
ATOM      9  CZ  ARG A  21       3.771   7.248   7.637  0.70 20.00           C
ATOM     10  NH1 ARG A  21       4.282   6.909   6.460  0.70 20.00           N
ATOM     11  NH2 ARG A  21       2.461   7.180   7.830  0.70 20.00           N
TER
HETATM    1  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor1 = """\
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
HETATM   12  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER
END
"""

pdb_poor2 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       7.923   8.820  10.129  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.312   7.441   9.621  0.70 20.00           C
ATOM      8  NE  ARG A  21       8.694   6.545  10.708  0.70 20.00           N
ATOM      9  CZ  ARG A  21       7.839   5.785  11.385  0.70 20.00           C
ATOM     10  NH1 ARG A  21       6.546   5.811  11.088  0.70 20.00           N
ATOM     11  NH2 ARG A  21       8.275   5.000  12.360  0.70 20.00           N
TER
HETATM    1  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER
END
"""

pdb_poor3 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       8.318  11.834   9.960  1.00 10.00           N
ATOM      2  CA  ARG A  21       7.146  11.154   9.422  1.00 10.00           C
ATOM      3  C   ARG A  21       6.012  11.120  10.440  1.00 10.00           C
ATOM      4  O   ARG A  21       5.000  10.449  10.235  1.00 10.00           O
ATOM      5  CB  ARG A  21       7.505   9.732   8.987  1.00 10.00           C
ATOM      6  CG  ARG A  21       8.469   9.666   7.814  0.70 20.00           C
ATOM      7  CD  ARG A  21       8.724   8.229   7.388  0.70 20.00           C
ATOM      8  NE  ARG A  21       9.667   8.147   6.277  0.70 20.00           N
ATOM      9  CZ  ARG A  21       9.331   8.301   5.000  0.70 20.00           C
ATOM     10  NH1 ARG A  21       8.070   8.545   4.668  0.70 20.00           N
ATOM     11  NH2 ARG A  21      10.255   8.208   4.054  0.70 20.00           N
TER      12      ARG A  21
HETATM    1  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER      14      ION B   1
END
"""

def exercise(pdb_poor_str, rotamer_manager, sin_cos_table, i_pdb, d_min = 1.0,
             resolution_factor = 0.1):
  # Fit one residue. There is a huge heavy atom nearby that overlaps with a
  # plausible rotamer.
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
  if(i_pdb==2): assert dist < 1.e-6, dist
  else:         assert dist < 0.22,  dist

if(__name__ == "__main__"):
  t0 = time.time()
  # load rotamer manager
  rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
  # pre-compute sin and cos tables
  sin_cos_table = scitbx.math.sin_cos_table(n=10000)
  for i_pdb, pdb_poor_str in enumerate(
                               [pdb_poor0, pdb_poor1, pdb_poor2, pdb_poor3]):
    exercise(
      pdb_poor_str    = pdb_poor_str,
      rotamer_manager = rotamer_manager,
      sin_cos_table   = sin_cos_table,
      i_pdb           = i_pdb)
  print "Time: %6.4f"%(time.time()-t0)
