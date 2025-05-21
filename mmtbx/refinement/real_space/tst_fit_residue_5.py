from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb

pdb_answer = """\
CRYST1   14.230   10.991   17.547  90.00  90.00  90.00 P 1
ATOM     43  N   GLY A   8       3.000   7.369   9.591  1.00 10.00           N
ATOM     44  CA  GLY A   8       4.304   7.980   9.410  1.00 10.00           C
ATOM     45  C   GLY A   8       5.375   6.966   9.058  1.00 10.00           C
ATOM     46  O   GLY A   8       5.074   5.814   8.745  1.00 10.00           O
ATOM     47  N   TYR A   9       6.631   7.398   9.110  1.00 10.00           N
ATOM     48  CA  TYR A   9       7.755   6.524   8.796  1.00 10.00           C
ATOM     49  C   TYR A   9       8.787   6.527   9.918  1.00 10.00           C
ATOM     50  O   TYR A   9       9.205   7.585  10.387  1.00 10.00           O
ATOM     51  CB  TYR A   9       8.410   6.945   7.478  1.00 10.00           C
ATOM     52  CG  TYR A   9       7.484   6.880   6.284  1.00 10.00           C
ATOM     54  CD1 TYR A   9       6.750   7.991   5.889  1.00 10.00           C
ATOM     53  CD2 TYR A   9       7.345   5.709   5.552  1.00 10.00           C
ATOM     56  CE1 TYR A   9       5.903   7.937   4.798  1.00 10.00           C
ATOM     55  CE2 TYR A   9       6.500   5.645   4.459  1.00 10.00           C
ATOM     57  CZ  TYR A   9       5.782   6.761   4.087  1.00 10.00           C
ATOM     58  OH  TYR A   9       4.940   6.702   3.000  1.00 10.00           O
ATOM     59  N   ASN A  10       9.193   5.336  10.345  1.00 10.00           N
ATOM     60  CA  ASN A  10      10.177   5.199  11.413  1.00 10.00           C
ATOM     61  C   ASN A  10      11.230   4.143  11.093  1.00 10.00           C
ATOM     62  O   ASN A  10      10.900   3.000  10.778  1.00 10.00           O
ATOM     63  CB  ASN A  10       9.486   4.870  12.738  1.00 10.00           C
ATOM     64  CG  ASN A  10      10.465   4.730  13.887  1.00 10.00           C
ATOM     65  OD1 ASN A  10      10.947   3.635  14.178  1.00 10.00           O
ATOM     66  ND2 ASN A  10      10.766   5.843  14.547  1.00 10.00           N
TER
END
"""

pdb_poor = """\
CRYST1   14.230   10.991   17.547  90.00  90.00  90.00 P 1
ATOM     43  N   GLY A   8       3.000   7.369   9.591  1.00 10.00           N
ATOM     44  CA  GLY A   8       4.304   7.980   9.410  1.00 10.00           C
ATOM     45  C   GLY A   8       5.375   6.966   9.058  1.00 10.00           C
ATOM     46  O   GLY A   8       5.074   5.814   8.745  1.00 10.00           O
ATOM     47  N   TYR A   9      10.261   7.310   2.887  1.00 10.00           N
ATOM     48  CA  TYR A   9      10.362   7.269   4.341  1.00 10.00           C
ATOM     49  C   TYR A   9      11.382   6.231   4.795  1.00 10.00           C
ATOM     50  O   TYR A   9      11.354   5.084   4.350  1.00 10.00           O
ATOM     51  CB  TYR A   9       8.996   6.971   4.965  1.00 10.00           C
ATOM     52  CG  TYR A   9       9.029   6.816   6.469  1.00 10.00           C
ATOM     53  CD1 TYR A   9       9.116   7.925   7.299  1.00 10.00           C
ATOM     54  CD2 TYR A   9       8.973   5.559   7.058  1.00 10.00           C
ATOM     55  CE1 TYR A   9       9.148   7.788   8.675  1.00 10.00           C
ATOM     56  CE2 TYR A   9       9.004   5.412   8.432  1.00 10.00           C
ATOM     57  CZ  TYR A   9       9.091   6.530   9.235  1.00 10.00           C
ATOM     58  OH  TYR A   9       9.122   6.388  10.604  1.00 10.00           O
ATOM     59  N   ASN A  10       9.193   5.336  10.345  1.00 10.00           N
ATOM     60  CA  ASN A  10      10.177   5.199  11.413  1.00 10.00           C
ATOM     61  C   ASN A  10      11.230   4.143  11.093  1.00 10.00           C
ATOM     62  O   ASN A  10      10.900   3.000  10.778  1.00 10.00           O
ATOM     63  CB  ASN A  10       9.486   4.870  12.738  1.00 10.00           C
ATOM     64  CG  ASN A  10      10.465   4.730  13.887  1.00 10.00           C
ATOM     65  OD1 ASN A  10      10.947   3.635  14.178  1.00 10.00           O
ATOM     66  ND2 ASN A  10      10.766   5.843  14.547  1.00 10.00           N
TER      25      ASN A  10
END
"""

def exercise(pdb_poor_str, i_pdb = 0, d_min = 1.0, resolution_factor = 0.25):
  # Fit one residue in many-residues model
  #
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor_str,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["GLY","TYR","ASN"],
    resolution_factor = resolution_factor)
  #
  get_class = iotbx.pdb.common_residue_names_get_class
  cmv = mmtbx.refinement.real_space.common_map_values(
    pdb_hierarchy = t.ph_poor,
    unit_cell     = t.xrs_poor.crystal_symmetry().unit_cell(),
    map_data      = t.target_map)
  for model in t.ph_poor.models():
    for chain in model.chains():
      for residue in chain.only_conformer().residues():
        if(get_class(residue.resname) == "common_amino_acid" and
           residue.resseq_as_int()==9): # take TYR9
          t0 = time.time()
          grm = t.model_poor.get_restraints_manager().geometry
          ro = mmtbx.refinement.real_space.fit_residue.run_with_minimization(
            target_map      = t.target_map,
            vdw_radii       = t.vdw,
            residue         = residue,
            xray_structure  = t.xrs_poor,
            mon_lib_srv     = t.mon_lib_srv,
            rotamer_manager = t.rotamer_manager,
            real_space_gradients_delta  = d_min*resolution_factor,
            cmv = cmv,
            geometry_restraints_manager = grm)
          sites_final = residue.atoms().extract_xyz()
          t1 = time.time()-t0
  #
  t.ph_poor.adopt_xray_structure(ro.xray_structure)
  t.ph_poor.write_pdb_file(file_name = "refined.pdb")
  #
  mmtbx.refinement.real_space.check_sites_match(
    ph_answer          = t.ph_answer,
    ph_refined         = t.ph_poor,
    exclude_atom_names = ["CD1","CD2","CE1","CE2"],
    tol                = 0.4)

if(__name__ == "__main__"):
  exercise(pdb_poor_str = pdb_poor, resolution_factor=0.2)
