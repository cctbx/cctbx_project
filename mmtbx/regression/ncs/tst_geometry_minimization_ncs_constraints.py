from __future__ import absolute_import, division, print_function

import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.model
from scitbx.array_family import flex
import iotbx.ncs
import iotbx.pdb
import sys
from libtbx import easy_run
from six.moves import range

pdb_string1 = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      1  N   THR A   1      11.582  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.571  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.546  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.768   9.193   9.237  1.00 10.00           N
ATOM      2  CA  THR B   1       7.284   8.754   7.945  1.00 10.00           C
ATOM      3  C   THR B   1       8.638   7.868   8.097  1.00 10.00           C
ATOM      4  O   THR B   1       9.495   8.426   8.852  1.00 10.00           O
ATOM      5  CB  THR B   1       7.423   9.832   6.963  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.144  10.446   6.765  1.00 10.00           O
ATOM      7  CG2 THR B   1       7.962   9.350   5.625  1.00 10.00           C
TER
ATOM      1  N   THR C   1       9.093   2.000  10.493  1.00 10.00           N
ATOM      2  CA  THR C   1       8.879   2.702   9.233  1.00 10.00           C
ATOM      3  C   THR C   1      10.081   3.570   8.875  1.00 10.00           C
ATOM      4  O   THR C   1      10.652   4.241   9.734  1.00 10.00           O
ATOM      5  CB  THR C   1       7.618   3.584   9.284  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.472   2.770   9.559  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.417   4.305   7.960  1.00 10.00           C
TER
"""

def tst_1(prefix="gm_ncs_constr_tst1"):
  log = sys.stdout
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_string1.split('\n'))
  # print dir(pdb_in)
  pdb_int_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pdb_int_params.pdb_interpretation.ncs_search.enabled = True
  model = mmtbx.model.manager(model_input = pdb_in)
  model.process(pdb_interpretation_params=pdb_int_params,
    make_restraints=True)
  ncs_obj = iotbx.ncs.input(hierarchy=model.get_hierarchy())
  ncs_restraints_group_list = ncs_obj.get_ncs_restraints_group_list()
  ncs_obj.show(format='phil')
  grm = model.get_restraints_manager()
  tmp_xrs = model.get_xray_structure().deep_copy_scatterers()
  refine_selection = flex.size_t(range(model.get_number_of_atoms()))

  # print "refining sites"
  cycle = 0
  tfg_obj = mmtbx.refinement.minimization_ncs_constraints.\
      target_function_and_grads_geometry_minimization(
          xray_structure=tmp_xrs,
          ncs_restraints_group_list=ncs_restraints_group_list,
          refine_selection=refine_selection,
          restraints_manager=grm.geometry,
          refine_sites=True,
          refine_transformations=False,
          )
  minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
    target_and_grads_object      = tfg_obj,
    xray_structure               = tmp_xrs,
    ncs_restraints_group_list    = ncs_restraints_group_list,
    refine_selection             = refine_selection,
    finite_grad_differences_test = False,
    max_iterations               = 100,
    refine_sites                 = True,
    refine_transformations       = False)
  refined_pdb_h = model.get_hierarchy().deep_copy()
  refined_pdb_h.adopt_xray_structure(tmp_xrs)
  refined_pdb_h.write_pdb_file("refined_%d.pdb" % cycle)
  new_ncs_obj = iotbx.ncs.input(hierarchy=refined_pdb_h)
  spec =  new_ncs_obj.get_ncs_info_as_spec()
  new_nrgl = new_ncs_obj.get_ncs_restraints_group_list()
  assert ncs_restraints_group_list == new_nrgl
  overall_rmsd_after = spec.overall_rmsd()
  assert overall_rmsd_after < 1e-6


def tst_2(prefix="gm_ncs_constr_tst2"):
  pdb_str = """\
CRYST1   18.415   14.419   12.493  90.00  90.00  90.00 P 1
ATOM      5  O   HOH S   3       5.648   5.054   4.347  1.00 10.00           O
ATOM      6  O   HOH S   4       2.000   5.459   5.854  1.00 10.00           O
ATOM      7  O   HOH S   5       3.498   4.230   8.986  1.00 10.00           O
TER
ATOM      1  N   THR A   1      11.782  12.419   4.645  1.00 10.00           N
ATOM      2  CA  THR A   1      11.671  11.061   4.125  1.00 10.00           C
ATOM      3  C   THR A   1      11.746  10.033   5.249  1.00 10.00           C
ATOM      4  O   THR A   1      12.561  10.157   6.163  1.00 10.00           O
ATOM      5  CB  THR A   1      12.772  10.760   3.092  1.00 10.00           C
ATOM      6  OG1 THR A   1      12.672  11.682   2.000  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.635   9.340   2.565  1.00 10.00           C
TER
ATOM      1  N   THR D   1      13.010   5.595  10.010  1.00 10.00           N
ATOM      2  CA  THR D   1      14.035   5.945   9.034  1.00 10.00           C
ATOM      3  C   THR D   1      15.310   6.423   9.720  1.00 10.00           C
ATOM      4  O   THR D   1      16.415   6.054   9.323  1.00 10.00           O
ATOM      5  CB  THR D   1      13.543   7.038   8.066  1.00 10.00           C
ATOM      6  OG1 THR D   1      13.237   8.229   8.802  1.00 10.00           O
ATOM      7  CG2 THR D   1      12.300   6.572   7.323  1.00 10.00           C
TER
ATOM      1  N   THR B   1       6.279   9.382   8.893  1.00 10.00           N
ATOM      2  CA  THR B   1       7.095   8.715   7.884  1.00 10.00           C
ATOM      3  C   THR B   1       8.190   7.873   8.531  1.00 10.00           C
ATOM      4  O   THR B   1       8.828   8.299   9.493  1.00 10.00           O
ATOM      5  CB  THR B   1       7.741   9.729   6.922  1.00 10.00           C
ATOM      6  OG1 THR B   1       6.718  10.485   6.262  1.00 10.00           O
ATOM      7  CG2 THR B   1       8.586   9.011   5.880  1.00 10.00           C
TER
ATOM      1  N   THR C   1       8.953   1.712   9.905  1.00 10.00           N
ATOM      2  CA  THR C   1       8.573   2.732   8.935  1.00 10.00           C
ATOM      3  C   THR C   1       9.744   3.658   8.623  1.00 10.00           C
ATOM      4  O   THR C   1      10.481   4.067   9.520  1.00 10.00           O
ATOM      5  CB  THR C   1       7.385   3.573   9.434  1.00 10.00           C
ATOM      6  OG1 THR C   1       6.258   2.721   9.672  1.00 10.00           O
ATOM      7  CG2 THR C   1       7.008   4.629   8.406  1.00 10.00           C
TER
ATOM      3  O   HOH E   1       4.021   6.412   7.953  1.00 10.00           O
ATOM      4  O   HOH E   2       4.186   2.319   6.140  1.00 10.00           O
TER      39      HOH E   2
END
"""
  f = open("%s_start.pdb" % prefix, 'w')
  f.write(pdb_str)
  f.close()
  #
  #  No NCS constraints - big rmsd
  #
  cmd = " ".join([
      "phenix.geometry_minimization",
      "%s_start.pdb" % prefix,
      # "ncs_search.enabled=True",
      "output_file_name_prefix=%s_minimized" % prefix])
  assert not easy_run.call(cmd)
  pdb_inp = iotbx.pdb.input(
      source_info=None, file_name="%s_minimized.pdb" % prefix)
  h = pdb_inp.construct_hierarchy()
  new_ncs_obj = iotbx.ncs.input(hierarchy=h)
  spec =  new_ncs_obj.get_ncs_info_as_spec()
  overall_rmsd_after = spec.overall_rmsd()
  assert overall_rmsd_after > 0.1, overall_rmsd_after # big rmsd
  #
  # NCS constraints - small rmsd
  #
  cmd = " ".join([
      "phenix.geometry_minimization",
      "%s_start.pdb" % prefix,
      "ncs_search.enabled=True",
      "output_file_name_prefix=%s_minimized_ncs" % prefix])
  assert not easy_run.call(cmd)
  pdb_inp = iotbx.pdb.input(
      source_info=None, file_name="%s_minimized_ncs.pdb" % prefix)
  h = pdb_inp.construct_hierarchy()
  new_ncs_obj = iotbx.ncs.input(hierarchy=h)
  spec =  new_ncs_obj.get_ncs_info_as_spec()
  overall_rmsd_after = spec.overall_rmsd()
  assert overall_rmsd_after < 0.001, overall_rmsd_after # small rmsd

def tst_3(prefix='gm_ncs_constr_tst3'):
  "Part of the model not in NCS"
  pdb_str = """\
CRYST1   35.937   25.866   35.477  90.00  90.00  90.00 P 21 21 21
ATOM      2  N   GLY A   1      -9.803   4.708   5.110  1.00 10.00           N
ATOM      3  CA  GLY A   1      -9.830   3.946   4.636  1.00 10.00           C
ATOM      4  C   GLY A   1      -8.132   2.550   4.372  1.00 10.00           C
ATOM      5  O   GLY A   1      -8.236   1.780   5.203  1.00 10.00           O
ATOM      6  N   ASN A   2      -7.622   3.238   3.060  1.00 10.00           N
ATOM      7  CA  ASN A   2      -6.590   1.914   3.237  1.00 10.00           C
ATOM      8  C   ASN A   2      -5.318   2.878   3.858  1.00 10.00           C
ATOM      9  O   ASN A   2      -5.365   3.945   3.936  1.00 10.00           O
ATOM     10  CB  ASN A   2      -5.853   2.413   1.389  1.00 10.00           C
ATOM     11  CG  ASN A   2      -6.894   1.495   0.829  1.00 10.00           C
ATOM     12  OD1 ASN A   2      -6.962   0.414   0.807  1.00 10.00           O
ATOM     13  ND2 ASN A   2      -7.645   2.660   0.235  1.00 10.00           N
ATOM     14  N   ASN A   3      -4.072   1.986   4.052  1.00 10.00           N
ATOM     15  CA  ASN A   3      -3.346   2.405   4.945  1.00 10.00           C
ATOM     16  C   ASN A   3      -1.584   1.504   4.409  1.00 10.00           C
ATOM     17  O   ASN A   3      -1.629   0.775   4.037  1.00 10.00           O
ATOM     18  CB  ASN A   3      -3.488   1.743   6.668  1.00 10.00           C
ATOM     19  CG  ASN A   3      -2.429   2.335   7.530  1.00 10.00           C
ATOM     20  OD1 ASN A   3      -1.942   3.238   7.987  1.00 10.00           O
ATOM     21  ND2 ASN A   3      -1.118   1.421   7.372  1.00 10.00           N
TER      48      ASN A   6
ATOM     47  N   TYR X   0      10.635   1.682   6.545  1.00 10.00           N
ATOM     48  CA  TYR X   0      11.465   1.717   7.567  1.00 10.00           C
ATOM     49  C   TYR X   0      13.059   1.541   6.953  1.00 10.00           C
ATOM     50  O   TYR X   0      13.237   1.134   5.945  1.00 10.00           O
ATOM     51  CB  TYR X   0      11.503   0.538   8.640  1.00 10.00           C
ATOM     52  CG  TYR X   0      10.281   0.695   9.372  1.00 10.00           C
ATOM     53  CD1 TYR X   0       9.090  -0.145   8.934  1.00 10.00           C
ATOM     54  CD2 TYR X   0      10.023   1.488  10.501  1.00 10.00           C
ATOM     55  CE1 TYR X   0       7.879   0.003   9.599  1.00 10.00           C
ATOM     56  CE2 TYR X   0       8.816   1.543  11.173  1.00 10.00           C
ATOM     57  CZ  TYR X   0       7.748   0.799  10.718  1.00 10.00           C
ATOM     58  OH  TYR X   0       6.545   0.850  11.384  1.00 10.00           O
ATOM     59  OXT TYR X   0      14.040   2.285   7.850  1.00 10.00           O
TER      62      TYR X   0
ATOM     61  N   GLY B   1      -9.774  11.509   5.385  1.00 30.00           N
ATOM     62  CA  GLY B   1      -9.683  10.144   4.557  1.00 30.00           C
ATOM     63  C   GLY B   1      -8.474   9.623   4.570  1.00 30.00           C
ATOM     64  O   GLY B   1      -8.306   8.679   5.018  1.00 30.00           O
ATOM     65  N   ASN B   2      -7.698   9.438   3.350  1.00 30.00           N
ATOM     66  CA  ASN B   2      -6.623   8.731   2.952  1.00 30.00           C
ATOM     67  C   ASN B   2      -5.306   9.804   3.583  1.00 30.00           C
ATOM     68  O   ASN B   2      -5.328  10.870   4.049  1.00 30.00           O
ATOM     69  CB  ASN B   2      -6.282   8.832   1.303  1.00 30.00           C
ATOM     70  CG  ASN B   2      -7.529   8.151   0.436  1.00 30.00           C
ATOM     71  OD1 ASN B   2      -7.543   7.027   0.515  1.00 30.00           O
ATOM     72  ND2 ASN B   2      -8.277   9.217  -0.129  1.00 30.00           N
ATOM     73  N   ASN B   3      -4.453   8.737   4.086  1.00 30.00           N
ATOM     74  CA  ASN B   3      -3.378   9.132   5.122  1.00 30.00           C
ATOM     75  C   ASN B   3      -1.816   8.645   4.171  1.00 30.00           C
ATOM     76  O   ASN B   3      -2.134   7.308   4.047  1.00 30.00           O
ATOM     77  CB  ASN B   3      -3.355   8.584   6.179  1.00 30.00           C
ATOM     78  CG  ASN B   3      -2.040   8.729   7.440  1.00 30.00           C
ATOM     79  OD1 ASN B   3      -2.045   9.966   7.506  1.00 30.00           O
ATOM     80  ND2 ASN B   3      -1.044   7.966   7.577  1.00 30.00           N
END
"""
  f = open("%s_start.pdb" % prefix, 'w')
  f.write(pdb_str)
  f.close()
  #
  #  No NCS constraints - big rmsd
  #
  cmd = " ".join([
      "phenix.geometry_minimization",
      "%s_start.pdb" % prefix,
      "output_file_name_prefix=%s_minimized" % prefix])
  assert not easy_run.call(cmd)
  pdb_inp = iotbx.pdb.input(
      source_info=None, file_name="%s_minimized.pdb" % prefix)
  h = pdb_inp.construct_hierarchy()
  new_ncs_obj = iotbx.ncs.input(hierarchy=h)
  spec =  new_ncs_obj.get_ncs_info_as_spec()
  overall_rmsd_after = spec.overall_rmsd()
  assert overall_rmsd_after > 0.01, overall_rmsd_after # big rmsd
  #
  # NCS constraints - small rmsd
  #
  cmd = " ".join([
      "phenix.geometry_minimization",
      "%s_start.pdb" % prefix,
      "ncs_search.enabled=True",
      "output_file_name_prefix=%s_minimized_ncs" % prefix])
  assert not easy_run.call(cmd)
  pdb_inp = iotbx.pdb.input(
      source_info=None, file_name="%s_minimized_ncs.pdb" % prefix)
  h = pdb_inp.construct_hierarchy()
  new_ncs_obj = iotbx.ncs.input(hierarchy=h)
  spec =  new_ncs_obj.get_ncs_info_as_spec()
  overall_rmsd_after = spec.overall_rmsd()
  assert overall_rmsd_after < 0.001, overall_rmsd_after # small rmsd

if __name__ == "__main__":
  tst_1()
  tst_2()
  tst_3()
