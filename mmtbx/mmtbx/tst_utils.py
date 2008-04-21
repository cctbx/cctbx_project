from cctbx.array_family import flex
import math, time, sys, os
from libtbx import adopt_init_args
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import mmtbx.model
from libtbx import introspection
import libtbx.load_env
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import model_statistics
from mmtbx import utils

def extract_serials(atoms, occ_groups):
  r = []
  for i in occ_groups:
    ri = []
    for j in i:
      ri.append([int(atoms[k].serial) for k in j])
    r.append(ri)
  return r

def exercise_00(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  xray_structure = processed_pdb_file.xray_structure()
  #
  base = [ [[2],[3]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[24,25,26,27],[28,29,30,31]] ]
  # default
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target = base[:]
  target.insert(3, [[21]])
  target.insert(4, [[23]])
  assert approx_equal(res, target)
  # default + add water
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  base_21_23 = target[:]
  target.extend([[[18]], [[19]], [[20]], [[22]]])
  assert approx_equal(res, target)
  # 1
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target = base_21_23[:]
  target.extend([[[0]], [[1]], [[4]], [[5]]])
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target.extend([[[18]], [[19]], [[20]], [[22]]])
  assert approx_equal(res, target)
  # 2
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target = base_21_23[:]
  target.extend([[[0, 1]], [[4, 5]]])
  assert approx_equal(res, target)
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False,
    other_group_selection_strings = ['resseq 0 and (name S or name O1)','resseq 0 and (name O3 or name O4)'])
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target.extend([[[18]], [[19]], [[20]], [[22]]])
  assert approx_equal(res, target)
  # 3
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and (name S or name O1)'],
    other_group_selection_strings = ['resseq 0 and (name O3 or name O4)'])
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target = base_21_23[:]
  target.extend([[[0]], [[1]], [[4, 5]]])
  assert approx_equal(res, target)

def exercise_01(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_h.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[0,1,2,3,4],[5,6,7,8,9]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_00(verbose=verbose)
  exercise_01(verbose=verbose)
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
