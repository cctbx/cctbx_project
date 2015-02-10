from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval

pdb_answer = """\
CRYST1   14.230   10.991   17.547  90.00  90.00  90.00 P 1
ATOM     43  N   GLY A   8       3.000   7.369   9.591  1.00  6.97           N
ATOM     44  CA  GLY A   8       4.304   7.980   9.410  1.00  7.38           C
ATOM     45  C   GLY A   8       5.375   6.966   9.058  1.00  6.50           C
ATOM     46  O   GLY A   8       5.074   5.814   8.745  1.00  6.23           O
ATOM     47  N   TYR A   9       6.631   7.398   9.110  1.00  6.94           N
ATOM     48  CA  TYR A   9       7.755   6.524   8.796  1.00  7.32           C
ATOM     49  C   TYR A   9       8.787   6.527   9.918  1.00  8.00           C
ATOM     50  O   TYR A   9       9.205   7.585  10.387  1.00  9.95           O
ATOM     51  CB  TYR A   9       8.410   6.945   7.478  1.00  8.02           C
ATOM     52  CG  TYR A   9       7.484   6.880   6.284  1.00  8.02           C
ATOM     53  CD1 TYR A   9       7.345   5.709   5.552  1.00  8.16           C
ATOM     54  CD2 TYR A   9       6.750   7.991   5.889  1.00  8.88           C
ATOM     55  CE1 TYR A   9       6.500   5.645   4.459  1.00  8.74           C
ATOM     56  CE2 TYR A   9       5.903   7.937   4.798  1.00  9.46           C
ATOM     57  CZ  TYR A   9       5.782   6.761   4.087  1.00  9.31           C
ATOM     58  OH  TYR A   9       4.940   6.702   3.000  1.00 11.26           O
ATOM     59  N   ASN A  10       9.193   5.336  10.345  1.00  7.07           N
ATOM     60  CA  ASN A  10      10.177   5.199  11.413  1.00  6.83           C
ATOM     61  C   ASN A  10      11.230   4.143  11.093  1.00  6.42           C
ATOM     62  O   ASN A  10      10.900   3.000  10.778  1.00  6.41           O
ATOM     63  CB  ASN A  10       9.486   4.870  12.738  1.00  7.59           C
ATOM     64  CG  ASN A  10      10.465   4.730  13.887  1.00  8.66           C
ATOM     65  OD1 ASN A  10      10.947   3.635  14.178  1.00  8.74           O
ATOM     66  ND2 ASN A  10      10.766   5.843  14.547  1.00 11.51           N
TER
END
"""

pdb_poor = """\
CRYST1   14.230   10.991   17.547  90.00  90.00  90.00 P 1
ATOM     43  N   GLY A   8       3.000   7.369   9.591  1.00  6.97           N
ATOM     44  CA  GLY A   8       4.304   7.980   9.410  1.00  7.38           C
ATOM     45  C   GLY A   8       5.375   6.966   9.058  1.00  6.50           C
ATOM     46  O   GLY A   8       5.074   5.814   8.745  1.00  6.23           O
ATOM     47  N   TYR A   9      10.261   7.310   2.887  1.00  6.94           N
ATOM     48  CA  TYR A   9      10.362   7.269   4.341  1.00  7.32           C
ATOM     49  C   TYR A   9      11.382   6.231   4.795  1.00  8.00           C
ATOM     50  O   TYR A   9      11.354   5.084   4.350  1.00  9.95           O
ATOM     51  CB  TYR A   9       8.996   6.971   4.965  1.00  8.02           C
ATOM     52  CG  TYR A   9       9.029   6.816   6.469  1.00  8.02           C
ATOM     53  CD1 TYR A   9       9.116   7.925   7.299  1.00  8.16           C
ATOM     54  CD2 TYR A   9       8.973   5.559   7.058  1.00  8.88           C
ATOM     55  CE1 TYR A   9       9.148   7.788   8.675  1.00  8.74           C
ATOM     56  CE2 TYR A   9       9.004   5.412   8.432  1.00  9.46           C
ATOM     57  CZ  TYR A   9       9.091   6.530   9.235  1.00  9.31           C
ATOM     58  OH  TYR A   9       9.122   6.388  10.604  1.00 11.26           O
ATOM     59  N   ASN A  10       9.193   5.336  10.345  1.00  7.07           N
ATOM     60  CA  ASN A  10      10.177   5.199  11.413  1.00  6.83           C
ATOM     61  C   ASN A  10      11.230   4.143  11.093  1.00  6.42           C
ATOM     62  O   ASN A  10      10.900   3.000  10.778  1.00  6.41           O
ATOM     63  CB  ASN A  10       9.486   4.870  12.738  1.00  7.59           C
ATOM     64  CG  ASN A  10      10.465   4.730  13.887  1.00  8.66           C
ATOM     65  OD1 ASN A  10      10.947   3.635  14.178  1.00  8.74           O
ATOM     66  ND2 ASN A  10      10.766   5.843  14.547  1.00 11.51           N
TER      25      ASN A  10
END
"""

def exercise(pdb_poor_str, d_min = 1.0, resolution_factor = 0.25):
  # Fit one residue in many-residues model
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
  # take TYR9
  sites_answer = list(
    pdb_inp.construct_hierarchy().residue_groups())[1].atoms().extract_xyz()
  # poor
  mon_lib_srv = monomer_library.server.server()
  master_params = iotbx.phil.parse(
    input_string=mmtbx.monomer_library.pdb_interpretation.master_params_str,
    process_includes=True).extract()
  master_params.link_distance_cutoff=999
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    params                   = master_params,
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = flex.std_string(pdb_poor_str.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor.pdb")
  grm = mmtbx.restraints.manager(
    geometry=processed_pdb_file.geometry_restraints_manager(show_energies=False),
    normalization = True)
  #
  rotamer_manager = RotamerEval()
  get_class = iotbx.pdb.common_residue_names_get_class
  for model in pdb_hierarchy_poor.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if(get_class(residue.resname) == "common_amino_acid" and
           int(residue.resseq)==9): # take TYR9
          t0 = time.time()
          ro = mmtbx.refinement.real_space.fit_residue.run(
            target_map      = target_map,
            residue         = residue,
            xray_structure  = xrs_poor,
            mon_lib_srv     = mon_lib_srv,
            rotamer_manager = rotamer_manager,
            real_space_gradients_delta  = d_min*resolution_factor,
            geometry_restraints_manager = grm)
          sites_final = residue.atoms().extract_xyz()
          t1 = time.time()-t0
  pdb_hierarchy_poor.adopt_xray_structure(ro.xray_structure)
  pdb_hierarchy_poor.write_pdb_file(file_name = "refined.pdb")
  dist = flex.mean(flex.sqrt((sites_answer - sites_final).dot()))
  return dist, t1

if(__name__ == "__main__"):
  dist, t = exercise(pdb_poor_str = pdb_poor, resolution_factor=0.2)
  print dist, t
  assert dist < 0.128, dist
