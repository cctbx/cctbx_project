from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.command_line.real_space_refine as rs

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

pdb_poor3 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       9.932  12.353  13.861  1.00 10.00           N
ATOM      2  CA  ARG A  21       9.452  10.985  14.024  1.00 10.00           C
ATOM      3  C   ARG A  21       9.990  10.362  15.308  1.00 10.00           C
ATOM      4  O   ARG A  21       9.369   9.469  15.883  1.00 10.00           O
ATOM      5  CB  ARG A  21       9.849  10.131  12.819  1.00 10.00           C
ATOM      6  CG  ARG A  21       9.292  10.629  11.495  1.00 10.00           C
ATOM      7  CD  ARG A  21       9.724   9.736  10.343  1.00 10.00           C
ATOM      8  NE  ARG A  21       9.197  10.200   9.063  1.00 10.00           N
ATOM      9  CZ  ARG A  21       9.825  11.057   8.264  1.00 10.00           C
ATOM     10  NH1 ARG A  21      11.007  11.547   8.611  1.00 10.00           N
ATOM     11  NH2 ARG A  21       9.270  11.425   7.117  1.00 10.00           N
TER
HETATM   12  U   ION B   1       9.074   7.848   5.000  1.00 10.00           U
TER
"""

pdb_poor4 = """\
CRYST1   14.074   16.834   17.360  90.00  90.00  90.00 P 1
ATOM      1  N   ARG A  21       9.932  12.353  13.861  1.00 10.00           N
ATOM      2  CA  ARG A  21       9.452  10.985  14.024  1.00 10.00           C
ATOM      3  C   ARG A  21       9.990  10.362  15.308  1.00 10.00           C
ATOM      4  O   ARG A  21       9.369   9.469  15.883  1.00 10.00           O
ATOM      5  CB  ARG A  21       9.849  10.131  12.819  1.00 10.00           C
ATOM      6  CG  ARG A  21       9.292  10.629  11.495  1.00 10.00           C
ATOM      7  CD  ARG A  21       9.724   9.736  10.343  1.00 10.00           C
ATOM      8  NE  ARG A  21       9.197  10.200   9.063  1.00 10.00           N
ATOM      9  CZ  ARG A  21       9.825  11.057   8.264  1.00 10.00           C
ATOM     10  NH1 ARG A  21      11.007  11.547   8.611  1.00 10.00           N
ATOM     11  NH2 ARG A  21       9.270  11.425   7.117  1.00 10.00           N
TER
"""

def exercise(pdb_poor_str, i_file, d_min = 1.0, resolution_factor = 0.1):
  # Fit one residue. There is a huge heavy atom nearby that overlaps with a
  # plausible rotamer.
  #
  # answer
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_answer)
  pdb_inp.write_pdb_file(file_name = "answer.pdb")
  xrs_answer = pdb_inp.xray_structure_simple()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer.mtz")
  sites_answer = list(
    pdb_inp.construct_hierarchy().residue_groups())[0].atoms().extract_xyz()
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
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor_%s.pdb"%str(i_file))
  grm = rs.get_geometry_restraints_manager(
    processed_pdb_file = processed_pdb_file,
    xray_structure     = xrs_poor)
  #
  rotamer_manager = RotamerEval()
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy_poor.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if(get_class(residue.resname) == "common_amino_acid"):
          t0 = time.time()
          ro = mmtbx.refinement.real_space.fit_residue.run(
            target_map      = target_map,
            residue         = residue,
            xray_structure  = xrs_poor,
            mon_lib_srv     = mon_lib_srv,
            rotamer_manager = rotamer_manager,
            real_space_gradients_delta  = d_min/4,
            geometry_restraints_manager = grm,
            anchor_site_cart_1=[8.318, 11.834,  9.960],
            anchor_site_cart_2=[6.012, 11.120, 10.440],
            anchor_atom_name_1="N",
            anchor_atom_name_2="C")
          sites_final = residue.atoms().extract_xyz()
          t1 = time.time()-t0
  pdb_hierarchy_poor.adopt_xray_structure(ro.xray_structure)
  pdb_hierarchy_poor.write_pdb_file(file_name = "refined_%s.pdb"%str(i_file))
  dist = flex.mean(flex.sqrt((sites_answer - sites_final).dot()))
  return dist, t1

if(__name__ == "__main__"):
  for i, pdb_poor in enumerate([pdb_poor1, pdb_poor2, pdb_poor3, pdb_poor4]):
    dist, t = exercise(pdb_poor_str = pdb_poor, i_file=i)
    print dist, t
    if(i in [0,1]): assert dist < 0.03
    else:           assert dist < 0.3
