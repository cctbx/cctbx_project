from __future__ import absolute_import, division, print_function
import time
from mmtbx.validation import validate_ligands
from mmtbx.validation.validate_ligands import master_params_str
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
import libtbx.load_env
from libtbx.test_utils import approx_equal
import os.path
from six.moves import zip

def tst_get_adps(vl_manager):
  '''
  Test getting ADPs of ligands and surrounding atoms
  '''
  n_iso_answer = (0,0,0,0,0)
  n_aniso_answer = (9,9,5,5,6)
  #print(vl_manager)
  for id_tuple, ligand_dict in vl_manager.items():
    #print(ligand_dict)
    for altloc, lr in ligand_dict.items():
      adps = lr.get_adps()
      id_str = lr.id_str
      if (id_str.strip() == 'A BEN    2 A'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 9)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [4.7, 8.1, 6.0], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.1,12.9,7.8], eps=0.1))
      if (id_str.strip() == 'A BEN    2 B'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 9)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [5.1, 8.2, 6.4], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.1,12.9,7.8], eps=0.1))
      if (id_str.strip() == 'A SO4    3'):
        assert(adps.n_iso == 5)
        assert(adps.n_aniso == 0)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [7.4,13.1,10.2], eps=0.1)
        assert(approx_equal(
          [adps.b_min_within, adps.b_max_within, adps.b_mean_within],
          [6.3, 11.1, 7.8], eps=0.1))
      if (id_str.strip() == 'A SO4    4'):
        assert(adps.n_iso == 0)
        assert(adps.n_aniso == 5)
        assert(adps.b_min_within is None)
        assert(adps.n_above_100 == 0)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [10.3,14.6,12.3], eps=0.1)
      if (id_str.strip() == 'A GOL    5'):
        assert(adps.n_iso == 6)
        assert(adps.n_aniso == 0)
        assert(adps.b_min_within is None)
        assert(adps.n_above_100 == 3)
        assert approx_equal([adps.b_min, adps.b_max, adps.b_mean],
          [58.7,114.9,96.9], eps=0.1)


def run_test1():
  '''
  Test if iselection for ligand PG5 (chain A resseq 201) is correct.
  '''
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/one_chain_ligand_water.pdb",
    test=os.path.isfile)

  pdb_inp = iotbx.pdb.input(pdb_fname)
  model = mmtbx.model.manager(model_input = pdb_inp, stop_for_unknowns=False)
  model.process(make_restraints=True)
  model.set_log(null_out())

  params = iotbx.phil.parse(
    input_string=master_params_str, process_includes=True).extract()
  # do not place H atoms for this test
  #params.validate_ligands.place_hydrogens = False

  vl_manager = validate_ligands.manager(
    model = model,
    fmodel = None,
    params = params.validate_ligands,
    log   = null_out())
  vl_manager.run()

  tst_get_ligands(vl_manager = vl_manager)
  tst_get_overlaps(vl_manager = vl_manager)


def tst_get_ligands(vl_manager):
  '''
  Test finding ligand
  '''
  assert (len(vl_manager) == 1)
  # test iselection
  for id_tuple, ligand_dict in vl_manager.items():
    assert (id_tuple == ('', 'A', ' 201'))
    lr = ligand_dict['']
    assert (list(lr.isel) == [84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95])


def tst_get_overlaps(vl_manager):
  '''
  Test nonbonded overlaps
  '''
  for id_tuple, ligand_dict in vl_manager.items():
    for altloc, lr in ligand_dict.items():
      clashes_result = lr.get_overlaps()
      assert(clashes_result.n_clashes == 5)
      assert approx_equal(clashes_result.clashscore, 31.6, eps=1.0)
# anaconda
#(['pdb=" HE3 MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 17, 54, 2.0370952358689647, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 15, 34, 2.946989989803154, 3.4, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 15, 56, 2.4839921497460486, 2.92, '', None)
#
# MAC
#(['pdb=" CE  MET A 107 "', 'pdb=" C8  PG5 A 201 "'], 16, 35, 2.946989989803154, 3.4, '', None),
#(['pdb=" HE3 MET A 107 "', 'pdb=" H83 PG5 A 201 "'], 18, 57, 2.026073542594147, 2.44, '', None),
#(['pdb=" CE  MET A 107 "', 'pdb=" H81 PG5 A 201 "'], 16, 55, 2.4973179613337146, 2.92, '', None)
#
# Readyset gives different names to H atoms.


def tst_get_occupancies(vl_manager):
  '''
  Test occupancy determination
  '''
  assert (len(vl_manager) == 4)
  id_tuple_answer = [('', 'A', '   2'), ('', 'A', '   3'), ('', 'A', '   4'), ('', 'A', '   5')]
  ligand_dict_length_answer = [2, 1, 1, 1]
  occupancy_answer = []
  # TODO six zip me
  for id_tuple, id_tuple_answer, length_answer in zip(vl_manager.keys(), id_tuple_answer, ligand_dict_length_answer):
    ligand_dict = vl_manager[id_tuple]
    assert (id_tuple == id_tuple_answer)
    assert (len(ligand_dict) == length_answer)
    for altloc, lr in ligand_dict.items():
      occs = lr.get_occupancies()
      id_str = lr.id_str
      if (id_str.strip() == 'A BEN    2 A'):
        assert(occs.occ_mean == 0.56)
      if (id_str.strip() == 'A BEN    2 B'):
        assert(occs.occ_mean == 0.44)
      if (id_str.strip() == 'A SO4    3'):
        assert(occs.occ_mean == 0.65)
      if (id_str.strip() == 'A SO4    4'):
        assert(occs.occ_mean == 0.48)
      if (id_str.strip() == 'A GOL    5'):
        assert(occs.occ_mean == 0.67)


def run_test2():
  '''
  Test
  - occupancy determination for ligands
  - adp determination for ligands and neighbors
  Tests are combined to decrease computing time (restraints manager is slow).
  '''
  pdb_fname = libtbx.env.find_in_repositories(
    relative_path="mmtbx/regression/pdbs/two_chains_ligand_water.pdb",
    test=os.path.isfile)

  pdb_inp = iotbx.pdb.input(pdb_fname)
  model = mmtbx.model.manager(model_input = pdb_inp)
  model.set_log(null_out())
  params = iotbx.phil.parse(
    input_string=master_params_str, process_includes=True).extract()
  vl_manager = validate_ligands.manager(
    model = model,
    fmodel = None,
    params = params.validate_ligands,
    log   = null_out)
  vl_manager.run()

  tst_get_occupancies(vl_manager = vl_manager)
  tst_get_adps(vl_manager = vl_manager)


def run():
  run_test1()
  run_test2()


if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f" % (time.time()-t0))
