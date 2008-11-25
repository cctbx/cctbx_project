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

def make_up_other_constrained_groups_obj(selections):
  result = []
  class foo:
    def __init__(self, selection):
      self.selection=selection
  for sel in selections:
    result.append( foo(selection = sel) )
  return result


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
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[['resseq 0 and (name S or name O1)']], [['resseq 0 and (name O3 or name O4)']]])
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_constrained_groups = other_constrained_groups)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target = base_21_23[:]
  target.extend([[[0, 1]], [[4, 5]]])
  assert approx_equal(res, target)
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[['resseq 0 and (name S or name O1)']], [['resseq 0 and (name O3 or name O4)']]])
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    add_water         = True,
    as_flex_arrays    = False,
    other_constrained_groups = other_constrained_groups)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  target.extend([[[18]], [[19]], [[20]], [[22]]])
  assert approx_equal(res, target)
  # 3
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[['resseq 0 and (name O3 or name O4)']]])
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and (name S or name O1)'],
    other_constrained_groups = other_constrained_groups)
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

def exercise_02(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_h.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[0,1,2,3,4,10,12,14,16,18,20,22], [5,6,7,8,9,11,13,15,17,19,21,23]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    ignore_hydrogens  = False,
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_03(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_hd.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[7]], [[8,9,10,11], [12,13,14]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    ignore_hydrogens  = False,
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_04(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_hd.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[8]], [[9],[12]], [[10],[13]], [[11],[14]], [[7]] ]
  res = utils.occupancy_selections(
    all_chain_proxies   = processed_pdb_file.all_chain_proxies,
    xray_structure      = processed_pdb_file.xray_structure(),
    ignore_hydrogens    = True,
    expect_exangable_hd = True,
    as_flex_arrays      = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_05(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_lys_arg_ser_tyr_neutron_hd.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[9],[12]],  [[10],[13]], [[11],[14]], [[33],[37]], [[34],[38]],
           [[35],[39]], [[36],[40]], [[59],[65]], [[60],[66]], [[61],[67]],
           [[62],[68]], [[63],[69]], [[64],[70]], [[80],[82]], [[81],[83]],
           [[103],[105]], [[104],[106]]]
  res = utils.occupancy_selections(
    all_chain_proxies   = processed_pdb_file.all_chain_proxies,
    xray_structure      = processed_pdb_file.xray_structure(),
    ignore_hydrogens    = True,
    expect_exangable_hd = True,
    as_flex_arrays      = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_06(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/NAD_594_HD.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log = log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[65],[77]], [[66],[78]], [[67],[79]], [[68],[80]], [[69],[81]],
           [[70],[82]], [[71],[83]], [[72],[84]], [[73],[85]], [[74],[86]],
           [[75],[87]], [[76],[88]],
           [[124],[127]], [[125],[128]], [[126],[129]],
           [[62]], [[113]] ]
  res = utils.occupancy_selections(
    all_chain_proxies   = processed_pdb_file.all_chain_proxies,
    xray_structure      = processed_pdb_file.xray_structure(),
    ignore_hydrogens    = True,
    expect_exangable_hd = True,
    as_flex_arrays      = False)
  assert approx_equal(res, base)

def exercise_07(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[0, 1, 2, 3, 4]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[['resseq 0']]])
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_08(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answers = [
    [ [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[0,1,2,3,4,5]] ],
    [ [[4],[5]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[6,7,8,9,10,11,12,13,14,15]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[16,17]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[18,19,20]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[23]], [[24,25,26,27],[28,29,30,31]], [[21]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[22]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23,24,25,26,27,28,29,30,31]] ],
    [ [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]] ]
  ]
  group_selections = ['resseq 0',
                      'resseq 1',
                      'resseq 2',
                      'resseq 3',
                      'resseq 4',
                      'resseq 5',
                      'resseq 6',
                      'resseq 0:6',
                      'not all']
  for group_selection, answer in zip(group_selections, answers):
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [[[group_selection]]])
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      other_constrained_groups = other_constrained_groups,
      as_flex_arrays    = False)
    assert approx_equal(result, answer)

def exercise_09(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answers = [
    [ [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[0]], [[1]], [[2]], [[3]], [[4]], [[5]] ],
    [ [[4],[5]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[6]], [[7]], [[8]], [[9]], [[10]], [[11]], [[12]], [[13]], [[14]], [[15]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[16]], [[17]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[18]], [[19]], [[20]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[23]], [[24,25,26,27],[28,29,30,31]], [[21]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]], [[22]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24]], [[25]], [[26]], [[27]], [[28]], [[29]], [[30]], [[31]] ],
    [ [[0]],[[1]],[[2]],[[3]],[[4]],[[5]],[[6]],[[7]],[[8]],[[9]],[[10]],[[11]],[[12]],[[13]],[[14]],[[15]],[[16]],[[17]],[[18]],[[19]],[[20]],[[21]],[[22]],[[23]],[[24]],[[25]],[[26]],[[27]],[[28]],[[29]],[[30]],[[31]] ],
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24,25,26,27],[28,29,30,31]] ]
  ]
  individual_selections = ['resseq 0',
                           'resseq 1',
                           'resseq 2',
                           'resseq 3',
                           'resseq 4',
                           'resseq 5',
                           'resseq 6',
                           'resseq 0:6',
                           'not all']
  for individual_selection, answer in zip(individual_selections, answers):
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      other_individual_selection_strings = [individual_selection],
      as_flex_arrays    = False)
    assert approx_equal(result, answer)

def exercise_10(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  e = None
  try:
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [[['resseq 0']]])
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      other_constrained_groups = other_constrained_groups,
      other_individual_selection_strings = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception, e: pass
  assert e.__str__() == "Duplicate selection: same atoms selected for individual and group occupancy refinement."

def exercise_11(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  e = None
  try:
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      remove_selection = ['resseq 0'],
      other_individual_selection_strings = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception, e: pass
  assert e.__str__() == "Duplicate selection: occupancies of same atoms selected to be fixed and to be refined."
  e = None
  try:
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [[['resseq 0']]])
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      other_constrained_groups = other_constrained_groups,
      remove_selection = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception, e: pass
  assert e.__str__() == "Duplicate selection: occupancies of same atoms selected to be fixed and to be refined."

def exercise_12(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[4],[5]], [[16],[17]], [[21]], [[23,24,25,26,27,28,29,30,31]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[['resseq 6']]])
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    remove_selection = ['resseq 1'],
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)
  #
  answer = [ [[4],[5]], [[16],[17]], [[21]], [[23]], [[24]], [[25]], [[26]], [[27]], [[28]], [[29]], [[30]], [[31]] ]
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    remove_selection = ['resseq 1'],
    other_individual_selection_strings = ['resseq 6'],
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_00(verbose=verbose)
  exercise_01(verbose=verbose)
  exercise_02(verbose=verbose)
  exercise_03(verbose=verbose)
  exercise_04(verbose=verbose)
  exercise_05(verbose=verbose)
  exercise_06(verbose=verbose)
  exercise_07(verbose=verbose)
  exercise_08(verbose=verbose)
  exercise_09(verbose=verbose)
  exercise_10(verbose=verbose)
  exercise_11(verbose=verbose)
  exercise_12(verbose=verbose)
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
