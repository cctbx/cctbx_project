from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from libtbx.utils import format_cpu_times, null_out

pdb_str_1="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HELIX    1   1 THR A    1  THR A    2  1                                   6
HELIX    1   1 THR B    1  THR B    2  1                                   6
SHEET    1   A 2 THR A   1  THR A   3  0
SHEET    2   A 2 THR B   4  THR B   5 -1  O  THR B   4   N  THR A   2
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THR A   1       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THR A   1       5.075   5.694   8.249  1.00 20.00           C
ATOM      4  O   THR A   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THR A   5       6.101   6.421   6.092  1.00 20.00           C
ATOM      6  OG1 THR A   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THR A   7       5.964   5.000   5.565  1.00 20.00           C
TER
END
"""

pdb_str_2="""
HELIX    1   1 THR A    1  THR A    2  1                                   6
SHEET    1   A 2 THR A   1  THR A   3  0
SHEET    2   A 2 THR B   4  THR B   5 -1  O  THR B   4   N  THR A   2
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       5.111   8.080   7.645  1.00 20.00           N
ATOM      2  CA  THR A   1       5.000   6.722   7.125  1.00 20.00           C
ATOM      3  C   THR A   1       5.075   5.694   8.249  1.00 20.00           C
ATOM      4  O   THR B   4       5.890   5.818   9.163  1.00 20.00           O
ATOM      5  CB  THR B   5       6.101   6.421   6.092  1.00 20.00           C
ATOM      6  OG1 THR C   6       6.001   7.343   5.000  1.00 20.00           O
ATOM      7  CG2 THR C   7       5.964   5.000   5.565  1.00 20.00           C
TER
END
"""

def exercise_1():
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_1)
  model = mmtbx.model.manager(
      model_input = inp,
      log = null_out())
  model.process(make_restraints=True)
  assert model.get_number_of_atoms() == 21, model.get_number_of_atoms()
  assert model.get_hierarchy().atoms_size() == 21
  assert model.get_xray_structure().scatterers().size() == 21
  ss = model.get_ss_annotation()
  # print ss.as_pdb_str()
  # STOP()
  assert ss.get_n_helices() == 3
  # because the second strand contains chain B which is not in ATOM records
  # whole sheet got discarded.
  assert ss.get_n_sheets() == 0
  rm = model.get_restraints_manager()
  assert rm.geometry.pair_proxies().bond_proxies.simple.size() == 6
  # since No NCS was set, these functions return the whole thing and no
  # master selection
  assert model.get_master_hierarchy().atoms_size() == 21
  assert model.get_master_selection().size() == 0
  # print model.model_as_pdb()
  # print "="*40

  # Here we set NCS constraints
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_1)
  pdb_int_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pdb_int_params.pdb_interpretation.ncs_search.enabled=True
  model = mmtbx.model.manager(
      model_input = inp,
      log = null_out())
  model.process(pdb_interpretation_params = pdb_int_params)
  # model.get_xray_structure()
  assert not model.ncs_constraints_present()
  assert model.get_ncs_obj() is not None
  model.setup_ncs_constraints_groups()
  # print model.get_ncs_obj()
  assert model.ncs_constraints_present()
  assert model.get_master_hierarchy().atoms_size() == 7
  # print model.get_master_hierarchy().as_pdb_string()
  # print list(model.get_master_selection())
  assert list(model.get_master_selection()).count(True) == 7

def exercise_2():
  """
  Same as 1 but automatic NCS search procedure does not match short chains,
  in this case chains B,C, so they left out of NCS.
  Not clear if we should utilize MTRIX instead of searching for NCS
  because currently we don't output them and in consecutive runs NCS
  search would be utilized anyway, potentially yelding different groups.
  """
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str_2)
  pdb_int_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pdb_int_params.pdb_interpretation.ncs_search.enabled=True
  model = mmtbx.model.manager(
      model_input = inp,
      log = null_out())
  model.process(pdb_interpretation_params = pdb_int_params)
  # model.get_xray_structure()
  ss = model.get_ss_annotation()
  assert ss.get_n_helices() == 3
  assert ss.get_n_sheets() == 3

  assert not model.ncs_constraints_present()
  assert model.get_ncs_obj() is not None
  model.setup_ncs_constraints_groups()
  # print model.get_ncs_obj()
  assert model.ncs_constraints_present()
  assert model.get_master_hierarchy().atoms_size() == 15
  # print model.get_master_hierarchy().as_pdb_string()
  assert list(model.get_master_selection()).count(True) == 15

if (__name__ == "__main__"):
  exercise_1()
  exercise_2()
  print(format_cpu_times())
  print("OK")
