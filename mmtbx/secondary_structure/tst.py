
from mmtbx.secondary_structure import manager
from iotbx import file_reader
import libtbx.load_env
import cStringIO
import os

def exercise () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_file_h = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print "Skipping exercise(): input file not available."
    return False
  log = cStringIO.StringIO()
  pdb_in = file_reader.any_file(pdb_file_h, force_type="pdb").file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple()
  sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=sec_str_from_pdb_file)
  m.find_automatically(log=log)
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.bonds.size() == 109
  m.params.h_bond_restraints.substitute_n_for_h = True
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.flag_use_bond.count(True) == 106
  (frac_alpha, frac_beta) = m.calculate_structure_content()
  assert ("%.3f" % frac_alpha) == "0.643"
  assert ("%.3f" % frac_beta) == "0.075"
  del m
  # using KSDSSP
  if (not libtbx.env.has_module(name="ksdssp")):
    print "Skipping KSDSSP tests: ksdssp module not available."
  else:
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=None)
    m.find_automatically(log=log)
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.bonds.size() == 93
    m.params.h_bond_restraints.substitute_n_for_h = True
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.flag_use_bond.count(True) == 86
    (frac_alpha, frac_beta) = m.calculate_structure_content()
    assert ("%.3f" % frac_alpha) == "0.552"
    assert ("%.3f" % frac_beta) == "0.066"
    del m
    del pdb_hierarchy
    del xray_structure
  # without hydrogens
  pdb_in = file_reader.any_file(pdb_file, force_type="pdb").file_object
  pdb_hierarchy = pdb_in.construct_hierarchy()
  xray_structure = pdb_in.xray_structure_simple()
  sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
  m = manager(pdb_hierarchy=pdb_hierarchy,
    xray_structure=xray_structure,
    sec_str_from_pdb_file=sec_str_from_pdb_file)
  m.find_automatically(log=log)
  bonds_table = m.get_bonds_table(log=log)
  assert bonds_table.bonds.size() == 109
  del m
  # using KSDSSP
  if (libtbx.env.has_module(name="ksdssp")):
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=None)
    m.find_automatically(log=log)
    bonds_table = m.get_bonds_table(log=log)
    assert bonds_table.bonds.size() == 93
  # Nucleic acids (requires REDUCE and PROBE)
  if (libtbx.env.has_module(name="reduce") and
      libtbx.env.has_module(name="probe")):
    pass # TODO
  print "OK"

if __name__ == "__main__" :
  exercise()
