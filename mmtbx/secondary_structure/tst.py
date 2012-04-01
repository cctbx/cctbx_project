
from mmtbx.secondary_structure import manager
from mmtbx.secondary_structure.base_pairing import pair_database
from mmtbx.geometry_restraints import hbond
from iotbx import file_reader
import libtbx.load_env
from libtbx.utils import null_out
import os

def exercise_protein () :
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_file_h = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print "Skipping exercise(): input file not available."
    return False
  if pdb_file_h is None :
    print "Skipping exercise(): input file not available."
    return False
  run_ksdssp = True
  if (not libtbx.env.has_module(name="ksdssp")):
    print "Skipping KSDSSP tests: ksdssp module not available."
    run_ksdssp = False
  log = null_out()
  expected_distances = [2.9, 1.975]
  for k, file_name in enumerate([pdb_file, pdb_file_h]) :
    pdb_in = file_reader.any_file(file_name, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    xray_structure = pdb_in.xray_structure_simple()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    m.find_automatically(log=log)
    m.params.h_bond_restraints.remove_outliers = False
    hbond_params = hbond.master_phil.extract()
    build_proxies = m.create_hbond_proxies(hbond_params=hbond_params,
      log=log)
    proxies = build_proxies.proxies
    assert (len(proxies) == len(build_proxies.exclude_nb_list) == 109)
    assert (type(proxies[0]).__name__ == "h_bond_simple_proxy")
    assert (proxies[0].distance_ideal == expected_distances[k])
    (frac_alpha, frac_beta) = m.calculate_structure_content()
    assert ("%.3f" % frac_alpha == "0.643")
    assert ("%.3f" % frac_beta == "0.075")
    # Make sure the hydrogen auto-detection override is working
    m.params.h_bond_restraints.substitute_n_for_h = True
    build_proxies = m.create_hbond_proxies(
      log=log,
      as_python_objects=True)
    proxies = build_proxies.proxies
    atom_ids = []
    for i_seq in proxies[0].i_seqs :
      atom_ids.append(pdb_hierarchy.atoms()[i_seq].id_str())
    assert (atom_ids == ['pdb=" N   ARG A  41 "', 'pdb=" O   ASP A  37 "'])
    assert (proxies[0].distance_ideal == 2.9)
    if (run_ksdssp) :
      m = manager(pdb_hierarchy=pdb_hierarchy,
        xray_structure=xray_structure,
        sec_str_from_pdb_file=None)
      m.find_automatically(log=log)
      m.params.h_bond_restraints.remove_outliers = False
      hbond_params.restraint_type = "simple"
      build_proxies = m.create_hbond_proxies(hbond_params=hbond_params, log=log)
      assert (build_proxies.proxies.size() == 81)

def exercise_nucleic_acids () :
  pdb_file_rna = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1u8d.pdb",
    test=os.path.isfile)
  pdb_file_dna = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/dnatest.pdb",
    test=os.path.isfile)
  if (pdb_file_rna is None) :
    print "Skipping exercise_nucleic_acids(): input file not available."
    return False
  # Nucleic acids (requires REDUCE and PROBE)
  if (libtbx.env.has_module(name="reduce") and
      libtbx.env.has_module(name="probe")):
    pdb_in = file_reader.any_file(pdb_file_rna, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    xray_structure = pdb_in.xray_structure_simple()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    log = null_out()
    m.find_automatically(log=log)
    build_proxies = m.create_hbond_proxies(log=log)
    assert (build_proxies.proxies.size() == 70)
    db = pair_database()
    assert (db.get_atoms("A-U", "WWT", True) == [('H61', 'O4'), ('N1', 'H3')])
    assert (db.get_atoms("DA-DT", "WWT", False) == [('N6', 'O4'), ('N1', 'N3')])
    assert (db.get_atoms("DT-DA", "WWT", False) == [('O4', 'N6'), ('N3', 'N1')])
    assert (db.get_atoms("G-C", "WWT", False) == [('O6', 'N4'), ('N1', 'N3'),
      ('N2', 'O2')])
    assert (db.get_atoms("C-G", "WWT", False) == [('N4', 'O6'), ('N3', 'N1'),
      ('O2', 'N2')])
    assert db.get_pair_type("A-U", [('H61', 'O4'), ('N1', 'H3')], True) == "XX"
    assert db.get_pair_type("A-U", [('N1', 'H3'), ('H61', 'O4')], True) == "XX"
    assert db.get_pair_type("C-G", [('N4', 'O6'), ('N3', 'N1'), ('O2', 'N2')], False) == "XIX"
    # DNA
    pdb_in = file_reader.any_file(pdb_file_dna, force_type="pdb").file_object
    pdb_hierarchy = pdb_in.construct_hierarchy()
    pdb_hierarchy.atoms().reset_i_seq()
    xray_structure = pdb_in.xray_structure_simple()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      sec_str_from_pdb_file=sec_str_from_pdb_file)
    log = null_out()
    m.find_automatically(log=log)
    assert (len(m.params.nucleic_acids.base_pair) == 4)
  else:
    print "Skipping base-pairing tests: reduce or probe module not available."
    pass
  print "OK"

if __name__ == "__main__" :
  exercise_protein()
  exercise_nucleic_acids()
