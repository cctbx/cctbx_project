from libtbx.test_utils import approx_equal
from smtbx.refinement.restraints import adp_restraints
from smtbx import development
from cctbx.array_family import flex
from cctbx import crystal

def get_pair_sym_table(xray_structure):
  asu_mappings = xray_structure.asu_mappings(buffer_thickness=3.5)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  scattering_types = xray_structure.scatterers().extract_scattering_types()
  pair_asu_table.add_covalent_pairs(
    scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
  return pair_asu_table.extract_pair_sym_table()

def exercise_adp_similarity():
  xray_structure = development.sucrose()
  pair_sym_table = get_pair_sym_table(xray_structure)
  for table in (None,pair_sym_table):
    if table is None: xs = xray_structure
    else: xs = None
    restraints = \
      adp_restraints.adp_similarity_restraints(
        xray_structure=xs,
        pair_sym_table=table)
    assert restraints.proxies.size() == 24
    i_seqs = (9,14,28,32,36,38)
    restraints = \
      adp_restraints.adp_similarity_restraints(
        xray_structure=xs,
        pair_sym_table=table,
        i_seqs=i_seqs)
    expected_i_seqs = ((9,32),(14,36),(32,36),(36,38))
    expected_weights = (625,156.25,625,625)
    proxies = restraints.proxies
    assert proxies.size() == len(expected_i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs, expected_i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])
    # add more restraints to same shared proxy
    i_seqs = (3,23,40,42)
    restraints = \
      adp_restraints.adp_similarity_restraints(
        xray_structure=xs,
        pair_sym_table=table,
        proxies=proxies,
        i_seqs=i_seqs)
    expected_i_seqs = (
      (9,32),(14,36),(32,36),(36,38),(3,23),(40,42))
    expected_weights = (625,156.25,625,625,156.25,625)
    proxies = restraints.proxies
    assert proxies.size() == len(expected_i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs, expected_i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])

def exercise_rigid_bond():
  xray_structure = development.sucrose()
  pair_sym_table = get_pair_sym_table(xray_structure)
  for table in (None,pair_sym_table):
    if table is None: xs = xray_structure
    else: xs = None
    restraints = \
      adp_restraints.rigid_bond_restraints(
        xray_structure=xs,
        pair_sym_table=table)
    assert restraints.proxies.size() == 60
    i_seqs = (9,14,28,32,36,38)
    restraints = \
      adp_restraints.rigid_bond_restraints(
        xray_structure=xs,
        pair_sym_table=table,
        i_seqs=i_seqs)
    expected_i_seqs = (
      (9,32),(9,36),(14,36),(14,32),(14,38),(32,36),(32,38),(36,38))
    expected_weights = [10000]*len(expected_i_seqs)
    proxies = restraints.proxies
    assert proxies.size() == len(expected_i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs, expected_i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])
    # add more restraints to same shared proxy
    i_seqs = (10,40,42)
    restraints = \
      adp_restraints.rigid_bond_restraints(
        xray_structure=xs,
        pair_sym_table=table,
        proxies=proxies,
        i_seqs=i_seqs)
    expected_i_seqs = (
      (9,32),(9,36),(14,36),(14,32),(14,38),(32,36),
      (32,38),(36,38),(10,42),(10,40),(40,42))
    expected_weights = [10000]*len(expected_i_seqs)
    proxies = restraints.proxies
    assert proxies.size() == len(expected_i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs, expected_i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])

def exercise_isotropic_adp():
  xray_structure = development.sucrose()
  xray_structure.scatterers()[10].set_use_u_iso_only()
  pair_sym_table = get_pair_sym_table(xray_structure)
  for table in (None,pair_sym_table):
    restraints = \
      adp_restraints.isotropic_adp_restraints(
        xray_structure=xray_structure,
        pair_sym_table=table)
    assert restraints.proxies.size() == 22
    i_seqs = (9,14,28,32,36,38)
    expected_weights = (100,25,100,100,100,100)
    restraints = \
      adp_restraints.isotropic_adp_restraints(
        xray_structure=xray_structure,
        pair_sym_table=table,
        i_seqs=i_seqs)
    proxies = restraints.proxies
    assert proxies.size() == len(i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs[0], i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])
    # add more restraints to same shared proxy
    i_seqs = (3,5,42)
    restraints = \
      adp_restraints.isotropic_adp_restraints(
        xray_structure=xray_structure,
        pair_sym_table=table,
        proxies=proxies,
        i_seqs=i_seqs)
    expected_i_seqs = (9,14,28,32,36,38,3,5,42)
    expected_weights = (100,25,100,100,100,100,25,25,100)
    proxies = restraints.proxies
    assert proxies.size() == len(expected_i_seqs)
    for i in range(proxies.size()):
      assert approx_equal(proxies[i].i_seqs[0], expected_i_seqs[i])
      assert approx_equal(proxies[i].weight, expected_weights[i])

def run():
  exercise_isotropic_adp()
  exercise_rigid_bond()
  exercise_adp_similarity()
  print "OK"

if __name__ == "__main__":
  run()
