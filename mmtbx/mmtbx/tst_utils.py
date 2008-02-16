from cctbx.array_family import flex
import math, time, sys, os
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
import mmtbx.model
from libtbx import introspection
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import model_statistics
from mmtbx import utils

def exercise_00():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  processed_pdb_files_srv = utils.process_pdb_file_srv(log = StringIO())
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  xray_structure = processed_pdb_file.xray_structure()
  aal = processed_pdb_file.all_chain_proxies.stage_1.atom_attributes_list
  #
  base = [ [[2],[3]], [[8,9,10,6,7],[11,12,13,14,15]], [[16],[17]], [[24,25,26,27],[28,29,30,31]], [[21]], [[23]] ]
  # default
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False)
  assert approx_equal(res, base)
  # default + add water
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False)
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  # default + add H occupancies
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    as_flex_arrays    = False)
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[19]], [[20]], [[21]], [[23]]]
  assert approx_equal(res, target)
  # default + add H occupancies + add water
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    add_water         = True,
    as_flex_arrays    = False)
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  # 1
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0]], [[1]], [[4]], [[5]], [[21]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0]], [[1]], [[4]], [[5]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0]], [[1]], [[4]], [[5]], [[19]], [[20]], [[21]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    add_water         = True,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0]], [[1]], [[4]], [[5]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  # 2
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0, 1]], [[4, 5]], [[21]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0, 1]], [[4, 5]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0, 1]], [[4, 5]], [[19]], [[20]], [[21]], [[23]]]
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_hydrogens     = True,
    add_water         = True,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0, 1]], [[4, 5]], [[18]], [[19]], [[20]], [[21]], [[22]], [[23]]]
  assert approx_equal(res, target)
  # 3
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and (name S or name O1)'],
    other_group_selection_strings = ['resseq 0 and (name O3 or name O4)'])
  target = [[[2], [3]], [[8, 9, 10, 6, 7], [11, 12, 13, 14, 15]], [[16], [17]], [[24, 25, 26, 27], [28, 29, 30, 31]], [[0]], [[1]], [[4, 5]], [[21]], [[23]]]
  assert approx_equal(res, target)

def run():
  exercise_00()

if (__name__ == "__main__"):
  run()
