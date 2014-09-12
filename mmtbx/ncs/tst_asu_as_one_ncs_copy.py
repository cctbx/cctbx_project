from __future__ import division
import mmtbx.ncs
import iotbx.pdb
import mmtbx.utils
from scitbx.array_family import flex

ncs_1_copy="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
END
"""

full_asu="""\
CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
ATOM      1  N   THR B   1       3.097   7.753  15.237  1.00 20.00           N
ATOM      2  CA  THR B   1       3.613   7.314  13.945  1.00 20.00           C
ATOM      3  C   THR B   1       4.967   6.628  14.097  1.00 20.00           C
ATOM      4  O   THR B   1       5.824   7.086  14.852  1.00 20.00           O
ATOM      5  CB  THR B   1       3.752   8.492  12.963  1.00 20.00           C
ATOM      6  OG1 THR B   1       2.473   9.106  12.765  1.00 20.00           O
ATOM      7  CG2 THR B   1       4.291   8.010  11.625  1.00 20.00           C
TER
ATOM      1  N   THR C   1       5.422   0.660  16.493  1.00 20.00           N
ATOM      2  CA  THR C   1       5.208   1.362  15.233  1.00 20.00           C
ATOM      3  C   THR C   1       6.410   2.230  14.875  1.00 20.00           C
ATOM      4  O   THR C   1       6.982   2.901  15.734  1.00 20.00           O
ATOM      5  CB  THR C   1       3.947   2.244  15.284  1.00 20.00           C
ATOM      6  OG1 THR C   1       2.801   1.430  15.559  1.00 20.00           O
ATOM      7  CG2 THR C   1       3.746   2.965  13.960  1.00 20.00           C
TER
"""

def exercise_00():
  #
  of = open("ncs_1.pdb","w")
  print >> of, ncs_1_copy
  of.close()
  of = open("full_asu.pdb","w")
  print >> of, full_asu
  of.close()
  #
  xrs_1_copy = iotbx.pdb.input(source_info=None,
    lines=ncs_1_copy).xray_structure_simple()
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=full_asu)
  xrs_full_asu = pdb_inp.xray_structure_simple()
  cs = pdb_inp.crystal_symmetry_from_cryst1()
  ph = pdb_inp.construct_hierarchy()
  #
  o = mmtbx.ncs.asu_ncs_converter(pdb_hierarchy=ph)
  xrs = o.ph_first_chain.extract_xray_structure(crystal_symmetry=cs)
  mmtbx.utils.assert_xray_structures_equal(x1 = xrs, x2 = xrs_1_copy)
  ###
  o.update_sites_cart(sites_cart_master_ncs_copy = xrs_1_copy.sites_cart())
  x1 = o.pdb_hierarchy.extract_xray_structure(crystal_symmetry=cs)
  mmtbx.utils.assert_xray_structures_equal(
    x1  = xrs_full_asu,
    x2  = o.pdb_hierarchy.extract_xray_structure(crystal_symmetry=cs),
    eps = 1.e-4)
  o.write_pdb_file(file_name="one.pdb", crystal_symmetry=cs, mode="ncs")
  o.write_pdb_file(file_name="all.pdb", crystal_symmetry=cs, mode="asu")
  mmtbx.utils.assert_xray_structures_equal(
    x1  = xrs_full_asu,
    x2  = iotbx.pdb.input(file_name="all.pdb").xray_structure_simple(),
    eps = 1.e-4)
  mmtbx.utils.assert_xray_structures_equal(
    x1  = xrs_1_copy,
    x2  = iotbx.pdb.input(file_name="one.pdb").xray_structure_simple(),
    eps = 1.e-4)
  ###
  sh1 = flex.vec3_double(xrs_1_copy.sites_cart().size(), [1,2,3])
  shf = flex.vec3_double(xrs_full_asu.sites_cart().size(), [1,2,3])
  o.update_sites_cart(sites_cart_master_ncs_copy = xrs_1_copy.sites_cart()+sh1)
  tmp = xrs_full_asu.deep_copy_scatterers()
  tmp.set_sites_cart(sites_cart = xrs_full_asu.sites_cart()+shf)
  mmtbx.utils.assert_xray_structures_equal(
    x1  = tmp,
    x2  = o.pdb_hierarchy.extract_xray_structure(crystal_symmetry=cs),
    eps = 1.e-4)
  o.write_pdb_file(file_name="one.pdb", crystal_symmetry=cs, mode="ncs")
  o.write_pdb_file(file_name="all.pdb", crystal_symmetry=cs, mode="asu")
  #

if __name__ == "__main__":
  exercise_00()
