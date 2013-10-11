from __future__ import division
import sys, os
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.utils import format_cpu_times, null_out, Sorry
import libtbx.load_env
from cStringIO import StringIO
from mmtbx import utils
import iotbx.pdb
from scitbx.array_family import flex

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
    selections = [ ['resseq 0 and (name S or name O1)'], ['resseq 0 and (name O3 or name O4)'] ])
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
    selections = [ ['resseq 0 and (name S or name O1)'], ['resseq 0 and (name O3 or name O4)'] ])
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
    selections = [ ['resseq 0 and (name O3 or name O4)'] ])
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
  base = [ [[0,1,2,3,4,10,12,14,16,18,20,22], [5,6,7,8,9,11,13,15,17,19,21,23]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_02(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/occ_mix1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[0,1,2,3,4,5,6,7,8,9,10,11,12], [14,15,16,17,18,19,20,21,22,23,24,25,26]], [[13],[27]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
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
  base = [ [[7]], [[8]], [[9],[12]], [[10],[13]], [[11],[14]] ]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
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
    as_flex_arrays      = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_06(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/NAD_594_HD.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log = log,
    stop_for_unknowns=False)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [ [[62]], [[113]], [[65],[77]],  [[66],[78]],  [[67],[79]], [[68],[80]],
                            [[69],[81]],  [[70],[82]],  [[71],[83]], [[72],[84]],
                            [[73],[85]],  [[74],[86]],  [[75],[87]], [[76],[88]],
                            [[124],[127]],[[125],[128]],[[126],[129]]]
  res = utils.occupancy_selections(
    all_chain_proxies   = processed_pdb_file.all_chain_proxies,
    xray_structure      = processed_pdb_file.xray_structure(),
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
    selections = [ ['resseq 0'] ])
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
    [ [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]] ]
  ]
  group_selections = ['resseq 0',
                      'resseq 1',
                      'resseq 2',
                      'resseq 3',
                      'resseq 4',
                      'resseq 5',
                      'resseq 6',
                      'resseq 0:6']
  for group_selection, answer in zip(group_selections, answers):
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [ [group_selection] ])
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
    [ [[4],[5]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[21]], [[23]], [[24]], [[25]], [[26]], [[27]], [[28]], [[29]], [[30]], [[31]] ]
  ]
  individual_selections = ['resseq 0',
                           'resseq 1',
                           'resseq 2',
                           'resseq 3',
                           'resseq 4',
                           'resseq 5',
                           'resseq 6',
                           'resseq 0:6']
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
      selections = [ ['resseq 0'] ])
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
      selections = [ ['resseq 0'] ])
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
    selections = [ ['resseq 6'] ])
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

def exercise_13(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9]], [[10]], [[0],[1]], [[2],[3]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and name N','chain A and resseq 1 and name CA'],
                   ['chain A and resseq 1 and name C','chain A and resseq 1 and name O'] ]
    )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_14(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9]], [[10]], [[0,1,2],[3,4]], [[5],[6]], [[7]] ]

  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and (name N or name CA or name C)', 'chain A and resseq 1 and (name O or name CB)'],
                   ['chain A and resseq 1 and name CG','chain A and resseq 1 and name CD'],
                   ['chain A and resseq 1 and name CE'] ]
    )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_15(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9]], [[0,1,2],[10]], [[5,7]] ]

  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and (name N or name CA or name C)', 'chain S and resseq 1'],
                   ['chain A and resseq 1 and name CG or chain A and resseq 1 and name CE'] ]
    )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_16(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9],[10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [
      ['chain A and resseq 1 and name NZ and altloc A', 'chain A and resseq 1 and name NZ and altloc B', 'chain S and resseq 1'] ]
    )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_17(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8,9,10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [
      ['chain A and resseq 1 and name NZ and altloc A or chain A and resseq 1 and name NZ and altloc B or chain S and resseq 1'] ]
    )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_18(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9],[10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
   selections = [
    ['chain A and resseq 1 and name NZ and altloc A','chain A and resseq 1 and name NZ and altloc B','chain S and resseq 1 and altloc C']]
  )
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_19(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[8],[9],[10]] ]
  tmp = "chain A and resseq 1 and name XX and altloc A"
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[
      tmp,
      'chain A and resseq 1 and name NZ and altloc B',
      'chain S and resseq 1']])
  try:
    result = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure,
      other_constrained_groups = other_constrained_groups,
      as_flex_arrays    = False)
  except Exception, e: pass
  assert str(e) == \
    'Selection string results in empty selection (selects no atoms): "%s"' \
      % tmp

def exercise_20(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ile_2conf_h.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  xray_structure = processed_pdb_file.xray_structure()
  answer = [ [[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]] ]
  result = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = xray_structure,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_21(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_3.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [[[2], [3]],
         [[6, 7, 8, 9, 10], [11, 12, 13, 14, 15]],
         [[16], [17]],
         [[21]],
         [[23]],
         [[24, 25, 26, 27], [28, 29, 30, 31]],
         [[36]],
         [[47]],
         [[48]],
         [[49]],
         [[50]],
         [[51]],
         [[53]],
         [[56, 57, 58, 59]],
         [[60, 61, 62, 63]],
         [[64, 65, 66, 67, 68]],
         [[37], [40]],
         [[38], [41]],
         [[39], [42]],
         [[43, 44, 45, 46]]]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_22(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_4.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [[[0, 1, 2, 3, 8, 9, 10, 11, 12], [4, 5, 6, 7, 13, 14, 15, 16, 17]]]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_23(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_5.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    pdb_file_names = [pdb_file])
  #
  base = [[[1, 2, 3, 4, 5, 6]], [[7, 8, 9, 10, 11], [12, 13, 14, 15, 16]]]
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  res = extract_serials(processed_pdb_file.all_chain_proxies.pdb_atoms, res)
  assert approx_equal(res, base)

def exercise_24(verbose):
  pdb_str1="""\
CRYST1   10.707   11.101   13.552  90.00  90.00  90.00 P 1
ATOM      0  N  AALA A   9       3.452   6.807   3.508  0.19  9.33      A    N
ATOM      1  CA AALA A   9       4.572   6.204   4.211  0.19  9.82      A    C
ATOM      2  C  AALA A   9       4.165   5.990   5.664  0.19 10.34      A    C
ATOM      3  O  AALA A   9       3.000   6.165   6.021  0.19 10.96      A    O
ATOM      4  CB AALA A   9       5.792   7.098   4.116  0.19 10.31      A    C
ATOM      5  H  AALA A   9       3.466   7.667   3.487  0.19  8.78      A    H
ATOM      6  HA AALA A   9       4.802   5.351   3.810  0.19  9.23      A    H
ATOM      7  HB1AALA A   9       6.533   6.686   4.588  0.19  9.91      A    H
ATOM      8  HB2AALA A   9       6.031   7.221   3.184  0.19  9.91      A    H
ATOM      9  HB3AALA A   9       5.594   7.960   4.515  0.19  9.91      A    H
ATOM     10  N  BALA A   9       3.348   6.697   3.518  0.28  8.28      A    N
ATOM     11  CA BALA A   9       4.461   6.052   4.195  0.28  9.14      A    C
ATOM     12  C  BALA A   9       4.138   5.964   5.683  0.28  9.84      A    C
ATOM     13  O  BALA A   9       3.003   6.215   6.089  0.28 10.68      A    O
ATOM     14  CB BALA A   9       5.726   6.829   3.952  0.28  9.20      A    C
ATOM     15  H  BALA A   9       3.422   7.551   3.454  0.28  8.78      A    H
ATOM     16  HA BALA A   9       4.597   5.156   3.849  0.28  9.23      A    H
ATOM     17  HB1BALA A   9       6.465   6.395   4.406  0.28  9.91      A    H
ATOM     18  HB2BALA A   9       5.907   6.863   3.000  0.28  9.91      A    H
ATOM     19  HB3BALA A   9       5.623   7.731   4.294  0.28  9.91      A    H
ATOM     20  N  CALA A   9       3.608   6.763   3.402  0.28  8.32      A    N
ATOM     21  CA CALA A   9       4.617   6.060   4.177  0.28  9.56      A    C
ATOM     22  C  CALA A   9       4.219   6.081   5.651  0.28 10.15      A    C
ATOM     23  O  CALA A   9       3.126   6.528   6.006  0.28 10.64      A    O
ATOM     24  CB CALA A   9       5.981   6.684   3.973  0.28 10.39      A    C
ATOM     25  H  CALA A   9       3.801   7.579   3.210  0.28  8.78      A    H
ATOM     26  HA CALA A   9       4.671   5.139   3.876  0.28  9.23      A    H
ATOM     27  HB1CALA A   9       6.639   6.202   4.497  0.28  9.91      A    H
ATOM     28  HB2CALA A   9       6.220   6.639   3.034  0.28  9.91      A    H
ATOM     29  HB3CALA A   9       5.959   7.611   4.257  0.28  9.91      A    H
ATOM     30  N  DALA A   9       3.518   6.930   3.530  0.25  8.78      A    N
ATOM     31  CA DALA A   9       4.639   6.333   4.232  0.25  9.23      A    C
ATOM     32  C  DALA A   9       4.203   6.093   5.674  0.25 10.10      A    C
ATOM     33  O  DALA A   9       3.051   6.346   6.031  0.25 10.72      A    O
ATOM     34  CB DALA A   9       5.837   7.255   4.177  0.25  9.91      A    C
ATOM     35  H  DALA A   9       3.490   7.789   3.568  0.25  8.78      A    H
ATOM     36  HA DALA A   9       4.898   5.494   3.819  0.25  9.23      A    H
ATOM     37  HB1DALA A   9       6.581   6.848   4.648  0.25  9.91      A    H
ATOM     38  HB2DALA A   9       6.086   7.408   3.252  0.25  9.91      A    H
ATOM     39  HB3DALA A   9       5.614   8.101   4.595  0.25  9.91      A    H
ATOM     40  N   VAL A  10       5.119   5.606   6.502  1.00 11.13      A    N
ATOM     41  CA  VAL A  10       4.846   5.470   7.925  1.00 12.50      A    C
ATOM     42  C   VAL A  10       4.347   6.801   8.520  1.00 11.26      A    C
ATOM     43  O   VAL A  10       4.763   7.871   8.095  1.00 11.53      A    O
ATOM     44  HA  VAL A  10       4.118   4.835   8.017  1.00 12.50      A    H
ATOM     45  CB AVAL A  10       5.994   4.806   8.722  0.21 14.17      A    C
ATOM     46  CG1AVAL A  10       6.640   3.699   7.889  0.21 14.17      A    C
ATOM     47  CG2AVAL A  10       7.005   5.815   9.197  0.21 15.20      A    C
ATOM     48  H  AVAL A  10       5.926   5.421   6.269  0.19 11.13      A    H
ATOM     49  HB AVAL A  10       5.616   4.404   9.520  0.21 14.91      A    H
ATOM     50 HG11AVAL A  10       7.358   3.289   8.396  0.21 16.29      A    H
ATOM     51 HG12AVAL A  10       5.975   3.028   7.671  0.21 16.29      A    H
ATOM     52 HG13AVAL A  10       6.998   4.077   7.070  0.21 16.29      A    H
ATOM     53 HG21AVAL A  10       7.707   5.363   9.691  0.21 15.63      A    H
ATOM     54 HG22AVAL A  10       7.391   6.271   8.433  0.21 15.63      A    H
ATOM     55 HG23AVAL A  10       6.570   6.462   9.774  0.21 15.63      A    H
ATOM     56  CB BVAL A  10       6.135   4.987   8.645  0.79 14.91      A    C
ATOM     57  CG1BVAL A  10       6.081   5.228  10.144  0.79 16.28      A    C
ATOM     58  CG2BVAL A  10       6.351   3.507   8.360  0.79 15.63      A    C
ATOM     59  H  BVAL A  10       5.928   5.441   6.263  0.28 11.13      A    H
ATOM     60  HB BVAL A  10       6.879   5.504   8.299  0.79 14.91      A    H
ATOM     61 HG11BVAL A  10       6.902   4.913  10.552  0.79 16.29      A    H
ATOM     62 HG12BVAL A  10       5.978   6.177  10.316  0.79 16.29      A    H
ATOM     63 HG13BVAL A  10       5.328   4.748  10.522  0.79 16.29      A    H
ATOM     64 HG21BVAL A  10       7.156   3.205   8.809  0.79 15.63      A    H
ATOM     65 HG22BVAL A  10       5.590   3.000   8.685  0.79 15.63      A    H
ATOM     66 HG23BVAL A  10       6.445   3.372   7.404  0.79 15.63      A    H
ATOM     67  H  CVAL A  10       5.907   5.353   6.270  0.28 11.13      A    H
ATOM     68  H  DVAL A  10       5.903   5.349   6.260  0.25 11.13      A    H
TER
END
"""
  pdb_str2="""\
CRYST1   10.707   11.101   13.552  90.00  90.00  90.00 P 1
ATOM      0  N  AALA A   9       3.452   6.807   3.508  0.19  9.33      A    N
ATOM      1  CA AALA A   9       4.572   6.204   4.211  0.19  9.82      A    C
ATOM      2  C  AALA A   9       4.165   5.990   5.664  0.19 10.34      A    C
ATOM      3  O  AALA A   9       3.000   6.165   6.021  0.19 10.96      A    O
ATOM      4  CB AALA A   9       5.792   7.098   4.116  0.19 10.31      A    C
ATOM      5  D  AALA A   9       3.466   7.667   3.487  0.19  8.78      A    D
ATOM      6  DA AALA A   9       4.802   5.351   3.810  0.19  9.23      A    D
ATOM      7  DB1AALA A   9       6.533   6.686   4.588  0.19  9.91      A    D
ATOM      8  DB2AALA A   9       6.031   7.221   3.184  0.19  9.91      A    D
ATOM      9  DB3AALA A   9       5.594   7.960   4.515  0.19  9.91      A    D
ATOM     10  N  BALA A   9       3.348   6.697   3.518  0.28  8.28      A    N
ATOM     11  CA BALA A   9       4.461   6.052   4.195  0.28  9.14      A    C
ATOM     12  C  BALA A   9       4.138   5.964   5.683  0.28  9.84      A    C
ATOM     13  O  BALA A   9       3.003   6.215   6.089  0.28 10.68      A    O
ATOM     14  CB BALA A   9       5.726   6.829   3.952  0.28  9.20      A    C
ATOM     15  D  BALA A   9       3.422   7.551   3.454  0.28  8.78      A    D
ATOM     16  DA BALA A   9       4.597   5.156   3.849  0.28  9.23      A    D
ATOM     17  DB1BALA A   9       6.465   6.395   4.406  0.28  9.91      A    D
ATOM     18  DB2BALA A   9       5.907   6.863   3.000  0.28  9.91      A    D
ATOM     19  DB3BALA A   9       5.623   7.731   4.294  0.28  9.91      A    D
ATOM     20  N  CALA A   9       3.608   6.763   3.402  0.28  8.32      A    N
ATOM     21  CA CALA A   9       4.617   6.060   4.177  0.28  9.56      A    C
ATOM     22  C  CALA A   9       4.219   6.081   5.651  0.28 10.15      A    C
ATOM     23  O  CALA A   9       3.126   6.528   6.006  0.28 10.64      A    O
ATOM     24  CB CALA A   9       5.981   6.684   3.973  0.28 10.39      A    C
ATOM     25  D  CALA A   9       3.801   7.579   3.210  0.28  8.78      A    D
ATOM     26  DA CALA A   9       4.671   5.139   3.876  0.28  9.23      A    D
ATOM     27  DB1CALA A   9       6.639   6.202   4.497  0.28  9.91      A    D
ATOM     28  DB2CALA A   9       6.220   6.639   3.034  0.28  9.91      A    D
ATOM     29  DB3CALA A   9       5.959   7.611   4.257  0.28  9.91      A    D
ATOM     30  N  DALA A   9       3.518   6.930   3.530  0.25  8.78      A    N
ATOM     31  CA DALA A   9       4.639   6.333   4.232  0.25  9.23      A    C
ATOM     32  C  DALA A   9       4.203   6.093   5.674  0.25 10.10      A    C
ATOM     33  O  DALA A   9       3.051   6.346   6.031  0.25 10.72      A    O
ATOM     34  CB DALA A   9       5.837   7.255   4.177  0.25  9.91      A    C
ATOM     35  D  DALA A   9       3.490   7.789   3.568  0.25  8.78      A    D
ATOM     36  DA DALA A   9       4.898   5.494   3.819  0.25  9.23      A    D
ATOM     37  DB1DALA A   9       6.581   6.848   4.648  0.25  9.91      A    D
ATOM     38  DB2DALA A   9       6.086   7.408   3.252  0.25  9.91      A    D
ATOM     39  DB3DALA A   9       5.614   8.101   4.595  0.25  9.91      A    D
ATOM     40  N   VAL A  10       5.119   5.606   6.502  1.00 11.13      A    N
ATOM     41  CA  VAL A  10       4.846   5.470   7.925  1.00 12.50      A    C
ATOM     42  C   VAL A  10       4.347   6.801   8.520  1.00 11.26      A    C
ATOM     43  O   VAL A  10       4.763   7.871   8.095  1.00 11.53      A    O
ATOM     44  HA  VAL A  10       4.118   4.835   8.017  1.00 12.50      A    D
ATOM     45  CB AVAL A  10       5.994   4.806   8.722  0.21 14.17      A    C
ATOM     46  CG1AVAL A  10       6.640   3.699   7.889  0.21 14.17      A    C
ATOM     47  CG2AVAL A  10       7.005   5.815   9.197  0.21 15.20      A    C
ATOM     48  D  AVAL A  10       5.926   5.421   6.269  0.19 11.13      A    D
ATOM     49  DB AVAL A  10       5.616   4.404   9.520  0.21 14.91      A    D
ATOM     50 DG11AVAL A  10       7.358   3.289   8.396  0.21 16.29      A    D
ATOM     51 DG12AVAL A  10       5.975   3.028   7.671  0.21 16.29      A    D
ATOM     52 DG13AVAL A  10       6.998   4.077   7.070  0.21 16.29      A    D
ATOM     53 DG21AVAL A  10       7.707   5.363   9.691  0.21 15.63      A    D
ATOM     54 DG22AVAL A  10       7.391   6.271   8.433  0.21 15.63      A    D
ATOM     55 DG23AVAL A  10       6.570   6.462   9.774  0.21 15.63      A    D
ATOM     56  CB BVAL A  10       6.135   4.987   8.645  0.79 14.91      A    C
ATOM     57  CG1BVAL A  10       6.081   5.228  10.144  0.79 16.28      A    C
ATOM     58  CG2BVAL A  10       6.351   3.507   8.360  0.79 15.63      A    C
ATOM     59  D  BVAL A  10       5.928   5.441   6.263  0.28 11.13      A    D
ATOM     60  DB BVAL A  10       6.879   5.504   8.299  0.79 14.91      A    D
ATOM     61 DG11BVAL A  10       6.902   4.913  10.552  0.79 16.29      A    D
ATOM     62 DG12BVAL A  10       5.978   6.177  10.316  0.79 16.29      A    D
ATOM     63 DG13BVAL A  10       5.328   4.748  10.522  0.79 16.29      A    D
ATOM     64 DG21BVAL A  10       7.156   3.205   8.809  0.79 15.63      A    D
ATOM     65 DG22BVAL A  10       5.590   3.000   8.685  0.79 15.63      A    D
ATOM     66 DG23BVAL A  10       6.445   3.372   7.404  0.79 15.63      A    D
ATOM     67  D  CVAL A  10       5.907   5.353   6.270  0.28 11.13      A    D
ATOM     68  D  DVAL A  10       5.903   5.349   6.260  0.25 11.13      A    D
TER
END
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  for pdb_str in [pdb_str1, pdb_str2]:
    processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
    processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
      raw_records = pdb_str.splitlines())
    res = utils.occupancy_selections(
      all_chain_proxies = processed_pdb_file.all_chain_proxies,
      xray_structure    = processed_pdb_file.xray_structure(),
      as_flex_arrays    = False)
    answer = \
      [[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 48],
        [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 59],
        [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 67],
        [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 68]],
       [[45, 46, 47, 49, 50, 51, 52, 53, 54, 55],
       [56, 57, 58, 60, 61, 62, 63, 64, 65, 66]]]
    assert approx_equal(res, answer)

def exercise_25(verbose):
  pdb_str="""\
CRYST1   10.707   11.101   13.552  90.00  90.00  90.00 P 1
ATOM      0  N   ALA A   9       3.452   6.807   3.508  1.00  9.33      A    N
ATOM      1  CA  ALA A   9       4.572   6.204   4.211  1.00  9.82      A    C
ATOM      2  C   ALA A   9       4.165   5.990   5.664  1.00 10.34      A    C
ATOM      3  O   ALA A   9       3.000   6.165   6.021  1.00 10.96      A    O
ATOM      4  CB  ALA A   9       5.792   7.098   4.116  1.00 10.31      A    C
ATOM      5  HA  ALA A   9       4.802   5.351   3.810  1.00  9.23      A    H
ATOM      6  HB1 ALA A   9       6.533   6.686   4.588  1.00  9.91      A    H
ATOM      7  HB2 ALA A   9       6.031   7.221   3.184  1.00  9.91      A    H
ATOM      8  HB3 ALA A   9       5.594   7.960   4.515  1.00  9.91      A    H
ATOM      9  H  AALA A   9       3.466   7.667   3.487  0.40  8.78      A    H
ATOM     10  D  BALA A   9       3.466   7.667   3.487  0.60  8.78      A    D
ATOM     11  N   VAL A  10       5.119   5.606   6.502  1.00 11.13      A    N
ATOM     12  CA  VAL A  10       4.846   5.470   7.925  1.00 12.50      A    C
ATOM     13  C   VAL A  10       4.347   6.801   8.520  1.00 11.26      A    C
ATOM     14  O   VAL A  10       4.763   7.871   8.095  1.00 11.53      A    O
ATOM     15  HA  VAL A  10       4.118   4.835   8.017  1.00 12.50      A    H
ATOM     16  CB  VAL A  10       5.994   4.806   8.722  1.00 14.17      A    C
ATOM     17  CG1 VAL A  10       6.640   3.699   7.889  1.00 14.17      A    C
ATOM     18  CG2 VAL A  10       7.005   5.815   9.197  1.00 15.20      A    C
ATOM     19  HB  VAL A  10       5.616   4.404   9.520  1.00 14.91      A    H
ATOM     20 HG11 VAL A  10       7.358   3.289   8.396  1.00 16.29      A    H
ATOM     21 HG12 VAL A  10       5.975   3.028   7.671  1.00 16.29      A    H
ATOM     22 HG13 VAL A  10       6.998   4.077   7.070  1.00 16.29      A    H
ATOM     23 HG21 VAL A  10       7.707   5.363   9.691  1.00 15.63      A    H
ATOM     24 HG22 VAL A  10       7.391   6.271   8.433  1.00 15.63      A    H
ATOM     25 HG23 VAL A  10       6.570   6.462   9.774  1.00 15.63      A    H
ATOM     26  H  AVAL A  10       5.926   5.421   6.269  0.30 11.13      A    H
ATOM     27  D  BVAL A  10       5.926   5.421   6.269  0.70 11.13      A    D
TER
END
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [ [[9],[10]], [[26],[27]] ]
  assert approx_equal(res, answer)

def exercise_26(verbose):
  pdb_str="""\
CRYST1   71.040   72.017   72.362  90.00 100.48  90.00 C 1 2 1
ATOM     96  N   PRO L   5       2.689  13.877  15.387  1.00 13.65           N
ATOM     97  CA  PRO L   5       1.824  14.762  14.572  1.00 17.31           C
ATOM     98  C   PRO L   5       0.338  14.432  14.641  1.00 20.79           C
ATOM     99  O   PRO L   5      -0.466  15.376  14.642  1.00 20.37           O
ATOM    100  CB  PRO L   5       2.330  14.534  13.143  1.00 20.71           C
ATOM    101  CG  PRO L   5       3.772  14.184  13.326  1.00 20.25           C
ATOM    102  CD  PRO L   5       3.871  13.403  14.633  1.00 16.57           C
ATOM    103  HA  PRO L   5       1.981  15.805  14.846  1.00 17.31           H
ATOM    104  HB2 PRO L   5       1.780  13.709  12.691  1.00 20.71           H
ATOM    105  HB3 PRO L   5       2.220  15.447  12.558  1.00 20.71           H
ATOM    106  HG2 PRO L   5       4.103  13.567  12.492  1.00 20.25           H
ATOM    107  HG3 PRO L   5       4.363  15.098  13.382  1.00 20.25           H
ATOM    108  HD2 PRO L   5       3.805  12.331  14.446  1.00 16.57           H
ATOM    109  HD3 PRO L   5       4.791  13.666  15.154  1.00 16.57           H
ATOM    110  N   LEU L   6      -0.052  13.175  14.677  1.00 13.93           N
ATOM    111  CA  LEU L   6      -1.446  12.769  14.667  1.00 15.53           C
ATOM    112  C   LEU L   6      -2.079  12.634  16.029  1.00 17.57           C
ATOM    113  O   LEU L   6      -3.268  12.311  16.111  1.00 18.17           O
ATOM    114  CB  LEU L   6      -1.648  11.435  13.889  1.00 17.76           C
ATOM    115  CG  LEU L   6      -1.291  11.544  12.396  1.00 18.22           C
ATOM    116  CD1 LEU L   6      -1.474  10.257  11.651  1.00 18.93           C
ATOM    117  CD2 LEU L   6      -2.125  12.629  11.689  1.00 22.55           C
ATOM    118  HA  LEU L   6      -2.017  13.534  14.144  1.00 15.53           H
ATOM    119  HB2 LEU L   6      -1.011  10.669  14.331  1.00 17.76           H
ATOM    120  HB3 LEU L   6      -2.693  11.135  13.959  1.00 17.76           H
ATOM    121  HG  LEU L   6      -0.242  11.827  12.310  1.00 18.22           H
ATOM    122 HD11 LEU L   6      -0.750  10.210  10.838  1.00 18.93           H
ATOM    123 HD12 LEU L   6      -1.319   9.426  12.338  1.00 18.93           H
ATOM    124 HD13 LEU L   6      -2.488  10.221  11.252  1.00 18.93           H
ATOM    125 HD21 LEU L   6      -2.084  12.462  10.613  1.00 22.55           H
ATOM    126 HD22 LEU L   6      -3.156  12.565  12.037  1.00 22.55           H
ATOM    127 HD23 LEU L   6      -1.712  13.609  11.929  1.00 22.55           H
ATOM    128  H  ALEU L   6       0.595  12.387  14.715  0.50 13.93           H
ATOM    129  D  BLEU L   6       0.595  12.387  14.715  0.50 13.93           D
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [ [[32], [33]] ]
  assert approx_equal(res, answer)

def exercise_27(verbose):
  pdb_str="""\
CRYST1   64.714   39.225   38.645  90.00 117.38  90.00 C 1 2 1
ATOM      0  N   SER A  -1      20.605   9.913  24.660  1.00 32.98           N
ATOM      1  CA  SER A  -1      21.415  10.057  23.431  1.00 25.22           C
ATOM      2  C   SER A  -1      20.514  10.247  22.233  1.00 25.05           C
ATOM      3  O   SER A  -1      19.332   9.926  22.266  1.00 28.08           O
ATOM      4  CB  SER A  -1      22.253   8.810  23.194  1.00 28.97           C
ATOM      5  OG  SER A  -1      21.417   7.708  22.900  1.00 37.21           O
ATOM      6  H1  SER A  -1      19.896  10.449  24.612  1.00 38.17           H
ATOM      7  H2  SER A  -1      20.335   9.069  24.737  1.00 27.38           H
ATOM      8  H3  SER A  -1      21.098  10.134  25.368  1.00 38.75           H
ATOM      9  HA  SER A  -1      21.997  10.829  23.514  1.00 12.22           H
ATOM     10  HB2 SER A  -1      22.844   8.970  22.440  1.00 22.78           H
ATOM     11  HB3 SER A  -1      22.771   8.614  23.990  1.00 30.47           H
ATOM     12  HG  SER A  -1      21.872   7.007  22.826  1.00 42.35           H
ATOM     13  N  AMET A   0      21.097  10.723  21.147  0.49 20.67           N
ATOM     14  CA AMET A   0      20.340  10.870  19.929  0.49 21.49           C
ATOM     15  C  AMET A   0      21.236  10.795  18.720  0.49 18.70           C
ATOM     16  O  AMET A   0      22.394  11.216  18.750  0.49 19.47           O
ATOM     17  CB AMET A   0      19.569  12.183  19.945  0.49 22.62           C
ATOM     18  CG AMET A   0      20.423  13.414  20.138  0.49 24.87           C
ATOM     19  SD AMET A   0      19.580  14.932  19.650  0.49 29.00           S
ATOM     20  CE AMET A   0      17.946  14.760  20.377  0.49 36.23           C
ATOM     21  H  AMET A   0      21.920  10.964  21.095  0.49 28.25           H
ATOM     22  HA AMET A   0      19.697  10.146  19.870  0.49  7.25           H
ATOM     23  HB2AMET A   0      19.093  12.280  19.105  0.49 13.51           H
ATOM     24  HB3AMET A   0      18.941  12.141  20.681  0.49  7.62           H
ATOM     25  HG2AMET A   0      20.671  13.490  21.072  0.49 26.02           H
ATOM     26  HG3AMET A   0      21.219  13.333  19.589  0.49 30.87           H
ATOM     27  HE1AMET A   0      17.284  14.819  19.669  0.49 20.79           H
ATOM     28  HE2AMET A   0      17.863  13.908  20.829  0.49  8.45           H
ATOM     29  HE3AMET A   0      17.812  15.481  21.012  0.49 30.25           H
ATOM     30  N  BMET A   0      21.082  10.809  21.171  0.51 21.19           N
ATOM     31  CA BMET A   0      20.368  11.023  19.923  0.51 23.13           C
ATOM     32  C  BMET A   0      21.273  10.654  18.766  0.51 21.10           C
ATOM     33  O  BMET A   0      22.496  10.703  18.893  0.51 19.93           O
ATOM     34  CB BMET A   0      19.961  12.488  19.782  0.51 27.15           C
ATOM     35  CG BMET A   0      19.070  12.993  20.889  0.51 29.67           C
ATOM     36  SD BMET A   0      18.685  14.739  20.684  0.51 41.63           S
ATOM     37  CE BMET A   0      17.734  15.043  22.171  0.51 35.23           C
ATOM     38  HA BMET A   0      19.568  10.476  19.897  0.51 36.28           H
ATOM     39  HB2BMET A   0      20.762  13.035  19.778  0.51  8.59           H
ATOM     40  HB3BMET A   0      19.485  12.602  18.945  0.51 27.25           H
ATOM     41  HG2BMET A   0      18.236  12.497  20.877  0.51 21.33           H
ATOM     42  HG3BMET A   0      19.519  12.877  21.741  0.51 34.36           H
ATOM     43  HE1BMET A   0      17.141  15.795  22.018  0.51 42.08           H
ATOM     44  HE2BMET A   0      17.217  14.249  22.380  0.51 22.21           H
ATOM     45  HE3BMET A   0      18.343  15.241  22.899  0.51 40.99           H
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [[[13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
             [30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45]]]
  assert approx_equal(res, answer)

def exercise_28(verbose):
  pdb_str="""\
CRYST1   64.360   64.360   46.038  90.00  90.00 120.00 P 63
ATOM      0  N   ASP A  48       8.896  25.394  -7.791  1.00  8.05           N
ATOM      1  CA  ASP A  48       8.495  26.452  -6.936  1.00  8.42           C
ATOM      2  C   ASP A  48       8.287  26.047  -5.477  1.00  8.20           C
ATOM      3  O   ASP A  48       8.309  26.881  -4.579  1.00 10.68           O
ATOM      4  CB  ASP A  48       7.216  27.151  -7.426  1.00  9.40           C
ATOM      5  CG  ASP A  48       7.457  27.744  -8.791  1.00 10.91           C
ATOM      6  OD1 ASP A  48       8.234  28.729  -8.836  1.00 16.64           O
ATOM      7  OD2 ASP A  48       6.845  27.293  -9.764  1.00 12.53           O
ATOM      8  HA  ASP A  48       9.193  27.122  -6.935  1.00  8.42           H
ATOM      9  HB2 ASP A  48       6.494  26.507  -7.490  1.00  9.40           H
ATOM     10  HB3 ASP A  48       6.981  27.867  -6.815  1.00  9.40           H
ATOM     11  H  AASP A  48       8.303  25.156  -8.367  0.50  8.04           H
ATOM     12  H  BASP A  48       8.242  25.041  -8.223  0.50  8.04           H
ATOM     13  N  ALEU A  49       8.083  24.740  -5.245  0.79  7.34           N
ATOM     14  CA ALEU A  49       7.817  24.239  -3.906  0.79  6.67           C
ATOM     15  C  ALEU A  49       8.124  22.738  -3.941  0.79  5.81           C
ATOM     16  O  ALEU A  49       7.880  22.074  -4.958  0.79  6.71           O
ATOM     17  CB ALEU A  49       6.385  24.559  -3.494  0.79  7.19           C
ATOM     18  CG ALEU A  49       5.914  24.092  -2.111  0.79  7.07           C
ATOM     19  CD1ALEU A  49       4.885  25.059  -1.536  0.79  8.84           C
ATOM     20  CD2ALEU A  49       5.323  22.713  -2.192  0.79  7.46           C
ATOM     21  H  ALEU A  49       8.095  24.131  -5.852  0.79  7.25           H
ATOM     22  HA ALEU A  49       8.421  24.661  -3.275  0.79  7.14           H
ATOM     23  HB2ALEU A  49       6.277  25.523  -3.518  0.79  9.16           H
ATOM     24  HB3ALEU A  49       5.791  24.158  -4.147  0.79  9.16           H
ATOM     25  HG ALEU A  49       6.673  24.062  -1.508  0.79  6.91           H
ATOM     26 HD11ALEU A  49       4.592  24.730  -0.672  0.79  9.95           H
ATOM     27 HD12ALEU A  49       5.294  25.933  -1.437  0.79  9.95           H
ATOM     28 HD13ALEU A  49       4.130  25.113  -2.143  0.79  9.95           H
ATOM     29 HD21ALEU A  49       4.960  22.476  -1.324  0.79  8.29           H
ATOM     30 HD22ALEU A  49       4.616  22.710  -2.856  0.79  8.29           H
ATOM     31 HD23ALEU A  49       6.015  22.082  -2.442  0.79  8.29           H
ATOM     32  N  BLEU A  49       7.975  24.768  -5.242  0.21  7.25           N
ATOM     33  CA BLEU A  49       7.654  24.205  -3.941  0.21  7.15           C
ATOM     34  C  BLEU A  49       8.003  22.716  -3.887  0.21  7.83           C
ATOM     35  O  BLEU A  49       7.689  22.025  -4.858  0.21  5.06           O
ATOM     36  CB BLEU A  49       6.162  24.365  -3.605  0.21  9.16           C
ATOM     37  CG BLEU A  49       5.681  23.652  -2.331  0.21  6.91           C
ATOM     38  CD1BLEU A  49       6.301  24.276  -1.095  0.21  9.95           C
ATOM     39  CD2BLEU A  49       4.156  23.640  -2.248  0.21  8.29           C
ATOM     40  H  BLEU A  49       7.943  24.178  -5.867  0.21  7.25           H
ATOM     41  HA BLEU A  49       8.173  24.662  -3.262  0.21  7.14           H
ATOM     42  HB2BLEU A  49       5.975  25.310  -3.494  0.21  9.16           H
ATOM     43  HB3BLEU A  49       5.645  24.021  -4.346  0.21  9.16           H
ATOM     44  HG BLEU A  49       5.963  22.725  -2.358  0.21  6.91           H
ATOM     45 HD11BLEU A  49       6.470  23.579  -0.443  0.21  9.95           H
ATOM     46 HD12BLEU A  49       7.132  24.697  -1.346  0.21  9.95           H
ATOM     47 HD13BLEU A  49       5.691  24.937  -0.731  0.21  9.95           H
ATOM     48 HD21BLEU A  49       3.888  23.174  -1.441  0.21  8.29           H
ATOM     49 HD22BLEU A  49       3.834  24.555  -2.225  0.21  8.29           H
ATOM     50 HD23BLEU A  49       3.802  23.184  -3.027  0.21  8.29           H
ATOM     51  N   VAL A  50       8.616  22.239  -2.807  1.00  5.93           N
ATOM     52  CA  VAL A  50       8.845  20.793  -2.609  1.00  5.53           C
ATOM     53  C   VAL A  50       7.981  20.307  -1.457  1.00  5.75           C
ATOM     54  O   VAL A  50       7.971  20.912  -0.389  1.00  6.63           O
ATOM     55  CB  VAL A  50      10.325  20.527  -2.343  1.00  6.31           C
ATOM     56  CG1 VAL A  50      10.556  19.043  -2.072  1.00  7.62           C
ATOM     57  CG2 VAL A  50      11.170  20.998  -3.512  1.00  7.52           C
ATOM     58  HA  VAL A  50       8.593  20.305  -3.404  1.00  5.53           H
ATOM     59  HB  VAL A  50      10.599  21.022  -1.555  1.00  6.31           H
ATOM     60 HG11 VAL A  50      11.507  18.860  -2.118  1.00  7.62           H
ATOM     61 HG12 VAL A  50      10.221  18.824  -1.188  1.00  7.62           H
ATOM     62 HG13 VAL A  50      10.087  18.523  -2.744  1.00  7.62           H
ATOM     63 HG21 VAL A  50      12.097  20.765  -3.345  1.00  7.52           H
ATOM     64 HG22 VAL A  50      10.860  20.562  -4.321  1.00  7.52           H
ATOM     65 HG23 VAL A  50      11.081  21.960  -3.600  1.00  7.52           H
ATOM     66  H  AVAL A  50       8.830  22.718  -2.125  0.79  5.93           H
ATOM     67  H  BVAL A  50       8.914  22.729  -2.166  0.21  5.93           H
TER
END
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [ [[11],[12]],
             [[13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,66],
              [32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,67]]]
  assert approx_equal(res, answer)

def exercise_29(verbose):
  pdb_str="""\
CRYST1  148.270   44.010   47.390  90.00 101.57  90.00 C 1 2 1
ATOM      0  N   GLY A 285     -41.269  16.430  -4.458  1.00 18.77           N
ATOM      1  CA  GLY A 285     -41.021  16.772  -5.854  1.00 20.45           C
ATOM      2  H   GLY A 285     -42.080  16.182  -4.313  1.00 22.53           H
ATOM      3  C  AGLY A 285     -41.133  18.291  -6.119  0.85 20.52           C
ATOM      4  O  AGLY A 285     -41.030  18.770  -7.258  0.85 22.89           O
ATOM      5  HA2AGLY A 285     -40.130  16.482  -6.104  0.85 24.54           H
ATOM      6  HA3AGLY A 285     -41.663  16.314  -6.418  0.85 24.54           H
ATOM      7  C  BGLY A 285     -40.556  18.155  -6.113  0.15 20.45           C
ATOM      8  O  BGLY A 285     -39.925  18.445  -7.127  0.15 21.06           O
ATOM      9  HA2BGLY A 285     -40.352  16.166  -6.208  0.15 24.54           H
ATOM     10  HA3BGLY A 285     -41.839  16.638  -6.357  0.15 24.54           H
ATOM     11  N  AASN A 286     -41.375  19.070  -5.066  0.75 20.63           N
ATOM     12  CA AASN A 286     -41.558  20.524  -5.179  0.75 21.34           C
ATOM     13  C  AASN A 286     -40.921  21.176  -3.941  0.75 19.76           C
ATOM     14  O  AASN A 286     -41.136  20.695  -2.825  0.75 18.94           O
ATOM     15  CB AASN A 286     -43.061  20.822  -5.246  0.75 23.19           C
ATOM     16  CG AASN A 286     -43.390  22.293  -5.087  0.75 24.76           C
ATOM     17  OD1AASN A 286     -43.580  22.784  -3.975  0.75 25.15           O
ATOM     18  ND2AASN A 286     -43.491  22.996  -6.206  0.75 26.38           N
ATOM     19  H  AASN A 286     -41.441  18.778  -4.260  0.75 24.76           H
ATOM     20  HA AASN A 286     -41.121  20.863  -5.988  0.75 25.61           H
ATOM     21  HB2AASN A 286     -43.400  20.532  -6.107  0.75 27.82           H
ATOM     22  HB3AASN A 286     -43.509  20.338  -4.535  0.75 27.82           H
ATOM     23 HD21AASN A 286     -43.371  22.614  -6.967  0.75 31.65           H
ATOM     24 HD22AASN A 286     -43.677  23.835  -6.171  0.75 31.65           H
ATOM     25  N  BASN A 286     -40.878  19.026  -5.184  0.25 20.30           N
ATOM     26  CA BASN A 286     -40.589  20.401  -5.396  0.25 20.20           C
ATOM     27  C  BASN A 286     -40.224  21.016  -4.085  0.25 18.88           C
ATOM     28  O  BASN A 286     -40.136  20.364  -3.047  0.25 18.65           O
ATOM     29  CB BASN A 286     -41.798  21.088  -6.023  0.25 22.27           C
ATOM     30  CG BASN A 286     -42.950  21.238  -5.058  0.25 23.28           C
ATOM     31  OD1BASN A 286     -42.781  21.720  -3.938  0.25 23.18           O
ATOM     32  ND2BASN A 286     -44.137  20.828  -5.491  0.25 24.35           N
ATOM     33  H  BASN A 286     -41.259  18.841  -4.435  0.25 24.36           H
ATOM     34  HA BASN A 286     -39.828  20.488  -6.007  0.25 24.24           H
ATOM     35  HB2BASN A 286     -41.538  21.974  -6.321  0.25 26.72           H
ATOM     36  HB3BASN A 286     -42.105  20.561  -6.777  0.25 26.72           H
ATOM     37 HD21BASN A 286     -44.216  20.499  -6.282  0.25 29.22           H
ATOM     38 HD22BASN A 286     -44.826  20.891  -4.981  0.25 29.22           H
ATOM     39  CA  GLU A 287     -39.388  22.905  -3.000  1.00 16.67           C
ATOM     40  C   GLU A 287     -40.376  23.372  -1.952  1.00 15.65           C
ATOM     41  O   GLU A 287     -40.132  23.201  -0.755  1.00 14.31           O
ATOM     42  CB  GLU A 287     -38.514  24.074  -3.481  1.00 17.80           C
ATOM     43  CG  GLU A 287     -37.273  23.645  -4.302  1.00 19.41           C
ATOM     44  CD  GLU A 287     -36.290  24.789  -4.558  1.00 20.84           C
ATOM     45  OE1 GLU A 287     -36.554  25.925  -4.128  1.00 21.26           O
ATOM     46  OE2 GLU A 287     -35.220  24.552  -5.185  1.00 22.93           O
ATOM     47  HB2 GLU A 287     -39.052  24.654  -4.041  1.00 21.36           H
ATOM     48  HB3 GLU A 287     -38.200  24.566  -2.707  1.00 21.36           H
ATOM     49  HG2 GLU A 287     -36.801  22.949  -3.818  1.00 23.29           H
ATOM     50  HG3 GLU A 287     -37.568  23.308  -5.163  1.00 23.29           H
ATOM     51  N  AGLU A 287     -40.109  22.235  -4.122  0.02 18.26           N
ATOM     52  H  AGLU A 287     -39.954  22.592  -4.889  0.02 21.91           H
ATOM     53  HA AGLU A 287     -38.796  22.250  -2.576  0.02 20.01           H
ATOM     54  N  BGLU A 287     -40.017  22.305  -4.119  0.98 18.44           N
ATOM     55  H  BGLU A 287     -40.228  22.836  -4.762  0.98 22.13           H
ATOM     56  HA BGLU A 287     -38.799  22.245  -2.580  0.98 20.01           H
TER
END
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [ [[3,4,5,6,19],
              [7,8,9,10,33]],
             [[11,12,13,14,15,16,17,18,20,21,22,23,24,52],
              [25,26,27,28,29,30,31,32,34,35,36,37,38,55]],
             [[51,53],
              [54,56]]]
  assert approx_equal(res, answer)

def exercise_30(verbose):
  pdb_str="""\
CRYST1   42.198  121.958   37.277  90.00  90.00  90.00 P 21 21 2
ATOM      0  CG  GLU A 115      30.700  22.521   0.401  0.55 25.56           C
ATOM      1  CD  GLU A 115      31.809  23.320  -0.265  1.00 25.96           C
ATOM      2  OE1 GLU A 115      32.842  22.797  -0.723  1.00 24.92           O
ATOM      3  OE2 GLU A 115      31.621  24.544  -0.376  1.00 27.30           O
ATOM      4  N  AGLU A 115      27.819  20.841  -1.012  0.44 19.61           N
ATOM      5  CA AGLU A 115      28.757  21.222  -0.004  0.44 20.79           C
ATOM      6  C  AGLU A 115      28.192  21.930   1.203  0.44 19.50           C
ATOM      7  O  AGLU A 115      27.475  22.922   1.098  0.44 20.38           O
ATOM      8  CB AGLU A 115      29.799  22.079  -0.601  0.44 23.59           C
ATOM      9  N  BGLU A 115      27.018  20.969  -0.446  0.56 27.49           N
ATOM     10  CA BGLU A 115      28.194  21.387   0.311  0.56 26.06           C
ATOM     11  C  BGLU A 115      27.541  21.859   1.611  0.56 25.00           C
ATOM     12  O  BGLU A 115      26.660  22.715   1.640  0.56 26.43           O
ATOM     13  CB BGLU A 115      29.189  22.459  -0.356  0.56 26.03           C
ATOM     14  N  AVAL A 116      28.585  21.407   2.363  0.53 19.29           N
ATOM     15  CA AVAL A 116      28.181  21.931   3.670  0.53 18.27           C
ATOM     16  C  AVAL A 116      29.427  21.990   4.589  0.53 17.81           C
ATOM     17  O  AVAL A 116      30.464  21.420   4.280  0.53 17.67           O
ATOM     18  CB AVAL A 116      27.090  21.046   4.342  0.53 20.31           C
ATOM     19  CG1AVAL A 116      25.743  21.168   3.633  0.53 22.78           C
ATOM     20  CG2AVAL A 116      27.498  19.598   4.395  0.53 20.85           C
ATOM     21  H  AVAL A 116      29.104  20.724   2.421  0.53 23.15           H
ATOM     22  HA AVAL A 116      27.827  22.838   3.564  0.53 21.92           H
ATOM     23  HB AVAL A 116      26.967  21.353   5.264  0.53 24.37           H
ATOM     24  N  BVAL A 116      27.987  21.231   2.690  0.47 21.87           N
ATOM     25  CA BVAL A 116      27.614  21.560   4.041  0.47 19.86           C
ATOM     26  C  BVAL A 116      28.915  21.857   4.746  0.47 19.34           C
ATOM     27  O  BVAL A 116      29.983  21.603   4.213  0.47 18.81           O
ATOM     28  CB BVAL A 116      26.938  20.336   4.707  0.47 19.81           C
ATOM     29  CG1BVAL A 116      25.591  20.061   4.058  0.47 21.33           C
ATOM     30  CG2BVAL A 116      27.825  19.086   4.627  0.47 19.25           C
ATOM     31  H  BVAL A 116      28.539  20.573   2.651  0.47 26.24           H
ATOM     32  HA BVAL A 116      27.021  22.340   4.070  0.47 23.83           H
ATOM     33  HB BVAL A 116      26.782  20.535   5.654  0.47 23.76           H
TER
END
"""
  if (verbose): log = sys.stdout
  else: log = StringIO()
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records = pdb_str.splitlines())
  res = utils.occupancy_selections(
    all_chain_proxies = processed_pdb_file.all_chain_proxies,
    xray_structure    = processed_pdb_file.xray_structure(),
    as_flex_arrays    = False)
  answer = [ [[0]],
             [[4, 5, 6, 7, 8, 21],
              [9, 10, 11, 12, 13, 31]],
             [[14, 15, 16, 17, 18, 19, 20, 22, 23],
              [24, 25, 26, 27, 28, 29, 30, 32, 33]] ]
  assert approx_equal(res, answer)

def exercise_d_data_target_d_atomic_params():
  import iotbx.pdb
  import mmtbx.f_model
  from scitbx.array_family import flex
  good = """
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  1.00 10.00           O
HETATM  115  O   HOH A  18       5.000   5.000   8.000  1.00 10.00           O
TER
END
  """
  bad = """
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  0.50 10.00           O
HETATM  115  O   HOH A  19       5.000   5.000   8.000  0.90 10.00           O
TER
END
  """
  for update_f_part1_for in [None, "map"]:
    pdb_inp = iotbx.pdb.input(source_info=None, lines=bad)
    xrs = pdb_inp.xray_structure_simple()
    xrs.scattering_type_registry(table = "wk1995")
    #
    xb = iotbx.pdb.input(source_info=None, lines=good).xray_structure_simple()
    xb.scattering_type_registry(table = "wk1995")
    f_obs = abs(xb.structure_factors(d_min=1,
      algorithm="direct").f_calc()).set_observation_type_xray_amplitude()
    #
    sfp = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    sfp.algorithm="direct"
    for target_name in ["ls_wunit_kunit", "ls_wunit_k1", "ml"]:
      print "%s"%target_name, "-"*50
      fmodel = mmtbx.f_model.manager(
        f_obs                        = f_obs,
        xray_structure               = xrs.deep_copy_scatterers(),
        target_name                  = target_name,
        sf_and_grads_accuracy_params = sfp)
      fmodel.update_all_scales(update_f_part1_for=None)
      alpha_beta = fmodel.alpha_beta()
      print "R-work: %6.4f"%fmodel.r_work()
      #
      tg = mmtbx.utils.experimental_data_target_and_gradients(fmodel = fmodel,
        alpha_beta=alpha_beta)
      tg.show()
      eps = 1.e-6
      # occupancies
      go = tg.grad_occ()
      for i in [0,1]:
        xrs1 = fmodel.xray_structure.deep_copy_scatterers()
        xrs2 = fmodel.xray_structure.deep_copy_scatterers()
        #
        xrs1.scatterers()[i].occupancy+=+eps
        tg.update_xray_structure(xray_structure=xrs1, alpha_beta=alpha_beta)
        t1 = tg.target()
        #
        xrs2.scatterers()[i].occupancy+=-eps
        tg.update_xray_structure(xray_structure=xrs2, alpha_beta=alpha_beta)
        t2 = tg.target()
        #
        gfd = (t1-t2)/(2*eps)
        print gfd, go[i]
        assert approx_equal(go[i], gfd, 1.e-5)
      # sites_cart
      gc = tg.grad_sites_cart()
      uc = fmodel.xray_structure.unit_cell()
      for i in [0,1]:
        gfd = []
        for e in [(eps,0,0),(0,eps,0),(0,0,eps)]:
          xrs1 = fmodel.xray_structure.deep_copy_scatterers()
          xrs2 = fmodel.xray_structure.deep_copy_scatterers()
          #
          s1 = flex.vec3_double([uc.orthogonalize(xrs1.scatterers()[i].site)])
          s1 = s1+flex.vec3_double([e])
          xrs1.scatterers()[i].site = uc.fractionalize(s1[0])
          tg.update_xray_structure(xray_structure=xrs1, alpha_beta=alpha_beta)
          t1 = tg.target()
          #
          s2 = flex.vec3_double([uc.orthogonalize(xrs2.scatterers()[i].site)])
          s2 = s2-flex.vec3_double([e])
          xrs2.scatterers()[i].site = uc.fractionalize(s2[0])
          tg.update_xray_structure(xray_structure=xrs2, alpha_beta=alpha_beta)
          t2 = tg.target()
          #
          gfd.append( (t1-t2)/(2*eps) )
        print gfd, list(gc[i])
        assert approx_equal(gc[i], gfd, 1.e-5)

def exercise_d_data_target_d_atomic_params2():
  import iotbx.pdb
  import mmtbx.f_model
  from scitbx.array_family import flex
  good = """
CRYST1   26.771   27.605   16.145  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1      21.771  22.605   5.434  1.00 10.00           N
ATOM      2  CA  ASP A   1      20.399  22.215   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1      19.546  21.744   6.168  1.00 10.00           C
ATOM      4  O   ASP A   1      19.997  20.958   7.001  1.00 10.00           O
ATOM      5  N   VAL A   2      18.313  22.233   6.229  1.00 10.00           N
ATOM      6  CA  VAL A   2      17.402  21.833   7.287  1.00 10.00           C
ATOM      7  C   VAL A   2      16.896  20.436   6.954  1.00 10.00           C
ATOM      8  O   VAL A   2      16.413  20.188   5.850  1.00 10.00           O
ATOM      9  N   GLN A   3      17.009  19.531   7.918  1.00 10.00           N
ATOM     10  CA  GLN A   3      16.590  18.147   7.740  1.00 10.00           C
ATOM     11  C   GLN A   3      15.154  17.904   8.207  1.00 10.00           C
ATOM     12  O   GLN A   3      14.813  18.175   9.356  1.00 10.00           O
ATOM     13  N   MET A   4      14.318  17.383   7.310  1.00 10.00           N
ATOM     14  CA  MET A   4      12.921  17.091   7.630  1.00 10.00           C
ATOM     15  C   MET A   4      12.750  15.585   7.849  1.00 10.00           C
ATOM     16  O   MET A   4      13.057  14.784   6.965  1.00 10.00           O
TER
ATOM     17  N   THR B   5      14.260  17.201  11.025  1.00 10.00           N
ATOM     18  CA  THR B   5      14.076  15.783  11.355  1.00 10.00           C
ATOM     19  C   THR B   5      12.612  15.405  11.550  1.00 10.00           C
ATOM     20  O   THR B   5      11.958  15.902  12.463  1.00 10.00           O
ATOM     21  N   GLN B   6      12.104  14.518  10.697  1.00 10.00           N
ATOM     22  CA  GLN B   6      10.712  14.084  10.797  1.00 10.00           C
ATOM     23  C   GLN B   6      10.571  12.704  11.424  1.00 10.00           C
ATOM     24  O   GLN B   6      11.336  11.787  11.118  1.00 10.00           O
ATOM     25  N   THR B   7       9.565  12.570  12.283  1.00 10.00           N
ATOM     26  CA  THR B   7       9.266  11.322  12.973  1.00 10.00           C
ATOM     27  C   THR B   7       7.756  11.127  12.915  1.00 10.00           C
ATOM     28  O   THR B   7       7.000  12.070  13.145  1.00 10.00           O
ATOM     29  N   PRO B   8       7.291   9.907  12.609  1.00 10.00           N
ATOM     30  CA  PRO B   8       8.060   8.698  12.310  1.00 10.00           C
ATOM     31  C   PRO B   8       8.406   8.619  10.827  1.00 10.00           C
ATOM     32  O   PRO B   8       8.058   9.514  10.054  1.00 10.00           O
ATOM     33  N   LEU B   9       9.093   7.553  10.427  1.00 10.00           N
ATOM     34  CA  LEU B   9       9.459   7.385   9.023  1.00 10.00           C
ATOM     35  C   LEU B   9       8.232   7.000   8.213  1.00 10.00           C
ATOM     36  O   LEU B   9       8.026   7.506   7.113  1.00 10.00           O
TER
  """
  bad = """
CRYST1   26.771   27.605   16.145  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1      21.771  22.605   5.434  1.00 10.00           N
ATOM      2  CA  ASP A   1      20.399  22.215   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1      19.546  21.744   6.168  1.00 10.00           C
ATOM      4  O   ASP A   1      19.997  20.958   7.001  1.00 10.00           O
ATOM      5  N   VAL A   2      18.313  22.233   6.229  1.00 10.00           N
ATOM      6  CA  VAL A   2      17.402  21.833   7.287  1.00 10.00           C
ATOM      7  C   VAL A   2      16.896  20.436   6.954  1.00 10.00           C
ATOM      8  O   VAL A   2      16.413  20.188   5.850  1.00 10.00           O
ATOM      9  N   GLN A   3      17.009  19.531   7.918  1.00 10.00           N
ATOM     10  CA  GLN A   3      16.590  18.147   7.740  1.00 10.00           C
ATOM     11  C   GLN A   3      15.154  17.904   8.207  1.00 10.00           C
ATOM     12  O   GLN A   3      14.813  18.175   9.356  1.00 10.00           O
ATOM     13  N   MET A   4      14.318  17.383   7.310  1.00 10.00           N
ATOM     14  CA  MET A   4      12.921  17.091   7.630  1.00 10.00           C
ATOM     15  C   MET A   4      12.750  15.585   7.849  1.00 10.00           C
ATOM     16  O   MET A   4      13.057  14.784   6.965  1.00 10.00           O
TER
ATOM     17  N   THR B   5      14.260  17.201  11.025  1.00 10.00           N
ATOM     18  CA  THR B   5      14.076  15.783  11.355  1.00 10.00           C
ATOM     19  C   THR B   5      12.612  15.405  11.550  1.00 10.00           C
ATOM     20  O   THR B   5      11.958  15.902  12.463  1.00 10.00           O
ATOM     21  N   GLN B   6      12.104  14.518  10.697  0.90 10.00           N
ATOM     22  CA  GLN B   6      10.712  14.084  10.797  0.90 10.00           C
ATOM     23  C   GLN B   6      10.571  12.704  11.424  0.90 10.00           C
ATOM     24  O   GLN B   6      11.336  11.787  11.118  0.90 10.00           O
ATOM     25  N   THR B   7       9.565  12.570  12.283  1.00 10.00           N
ATOM     26  CA  THR B   7       9.266  11.322  12.973  1.00 10.00           C
ATOM     27  C   THR B   7       7.756  11.127  12.915  1.00 10.00           C
ATOM     28  O   THR B   7       7.000  12.070  13.145  1.00 10.00           O
ATOM     29  N   PRO B   8       7.291   9.907  12.609  1.00 10.00           N
ATOM     30  CA  PRO B   8       8.060   8.698  12.310  1.00 10.00           C
ATOM     31  C   PRO B   8       8.406   8.619  10.827  1.00 10.00           C
ATOM     32  O   PRO B   8       8.058   9.514  10.054  1.00 10.00           O
ATOM     33  N   LEU B   9       9.093   7.553  10.427  1.10 10.00           N
ATOM     34  CA  LEU B   9       9.459   7.385   9.023  1.10 10.00           C
ATOM     35  C   LEU B   9       8.232   7.000   8.213  1.10 10.00           C
ATOM     36  O   LEU B   9       8.026   7.506   7.113  1.10 10.00           O
TER
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines=good)
  hierarchy = pdb_inp.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = "wk1995")
  f_obs = abs(xrs.structure_factors(d_min=1,
    algorithm="direct").f_calc()).set_observation_type_xray_amplitude()
  sfp = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfp.algorithm="direct"
  #
  pdb_inp = iotbx.pdb.input(source_info=None, lines=bad)
  xrs = pdb_inp.xray_structure_simple()
  xrs.scattering_type_registry(table = "wk1995")
  #
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    xray_structure               = xrs,
    target_name                  = "ml",
    sf_and_grads_accuracy_params = sfp)
  alpha_beta = fmodel.alpha_beta()
  tg = mmtbx.utils.experimental_data_target_and_gradients(
    fmodel = fmodel,
    alpha_beta=alpha_beta)
  result = tg.group_occupancy_grads(
    pdb_hierarchy       = hierarchy,
    residues_per_window = 1)
  for r in result:
    print "chainID_resseqs: %s occupancy_grad: %-15.6f"%tuple(r)
  # get group gradients using custom selections
  selections = [flex.size_t([0,1,2,3]), flex.size_t([4,5,6,7,8])]
  result = tg.group_occupancy_grads(selections = selections)
  for i, r in enumerate(result):
    print "selection#: %s occupancy_grad: %-15.6f"%(str(i), r[1])

def exercise_get_atom_selections (verbose=False) :
  pdb_in = """\
CRYST1   15.000   15.000   15.000  90.00  90.00  90.00 P 212121
HETATM  115  O   HOH A  18       3.000   5.000   5.000  1.00 10.00           O
HETATM  115  O   HOH A  19       5.000   5.000   8.000  1.00 10.00           O
HETATM  115  O   HOH A  20       5.000   5.000   8.000  1.00 10.00           O
END"""
  log = null_out()
  if (verbose) :
    log = sys.stdout
  processed_pdb_files_srv = utils.process_pdb_file_srv(log=log)
  processed_pdb_file, pdb_inp = processed_pdb_files_srv.process_pdb_files(
    raw_records=pdb_in.splitlines())
  selections1 = utils.get_atom_selections(
    all_chain_proxies=processed_pdb_file.all_chain_proxies,
    xray_structure=processed_pdb_file.xray_structure(),
    selection_strings=["resseq 18", "resseq 19", "resseq 20"],
    parameter_name="refine.occupancy")
  try :
    selections2 = utils.get_atom_selections(
      all_chain_proxies=processed_pdb_file.all_chain_proxies,
      xray_structure=processed_pdb_file.xray_structure(),
      selection_strings=["resseq 18:19", "resseq 19:20"],
      parameter_name="refine.occupancy")
  except Sorry, s :
    assert (str(s) == """\
One or more overlapping selections for refine.occupancy:
resseq 18:19
resseq 19:20""")
  else :
    raise Exception_expected

def exercise_f_000():
  pdb_str="""
CRYST1    5.000    5.000    5.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       0.000   0.000   0.000  1.00  5.00           C
END
"""
  xrs = iotbx.pdb.input(source_info=None, lines=pdb_str).xray_structure_simple()
  #
  f = utils.f_000(xray_structure=xrs, mean_solvent_density=0)
  assert approx_equal(f.f_000, 6, 1.e-3)
  #
  f = utils.f_000(mean_solvent_density=1, unit_cell_volume=125)
  assert approx_equal(f.f_000, 125)
  #
  f = utils.f_000(mean_solvent_density=0.25, unit_cell_volume=125,
                  solvent_fraction=0.3, xray_structure = xrs)
  assert approx_equal(f.f_000, 0.25*125*0.3+6, 1.e-3)
  assert approx_equal(f.solvent_fraction, 0.3)
  #
  f = utils.f_000(mean_solvent_density=0.25, unit_cell_volume=125,
                  xray_structure = xrs)
  assert approx_equal(f.f_000, 0.25*125*0.687355324074+6, 1.e-3)

def exercise_cmdline_load_pdb_and_data() :
  pdb_str="""
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
END
"""
  pdb_in = iotbx.pdb.input(source_info=None,lines=pdb_str)
  hierarchy = pdb_in.construct_hierarchy()
  xrs = pdb_in.xray_structure_simple()
  xrs.scattering_type_registry(
    d_min=1.5,
    table="n_gaussian")
  file_base = "tmp_mmtbx_utils"
  open(file_base+".pdb", "w").write(
    hierarchy.as_pdb_string(crystal_symmetry=xrs))
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = fc.generate_r_free_flags()
  mtz = fc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(file_base+".mtz")
  open(file_base+".fa", "w").write(">Tyr\nY\n")
  cmdline = utils.cmdline_load_pdb_and_data(
    update_f_part1_for=None,
    args=[ file_base + ext for ext in [".pdb",".mtz",".fa"] ],
    master_phil=utils.cmdline_input_phil(),
    out=null_out())
  assert (cmdline.params.input.xray_data.file_name is not None)
  assert (cmdline.sequence is not None)
  r_factor = cmdline.fmodel.r_work()
  assert (r_factor < 0.001)

def exercise_detect_link_problems () :
  open("tmp_mmtbx_utils_asn_nag.pdb", "w").write("""\
CRYST1  124.702  124.702   71.573  90.00  90.00  90.00 P 4 21 2
ATOM   3196  N   ASN A 284      36.622 -19.654  35.782  1.00 19.63           N
ATOM   3197  CA  ASN A 284      36.491 -18.279  35.327  1.00 19.79           C
ATOM   3198  C   ASN A 284      35.037 -17.835  35.322  1.00 20.05           C
ATOM   3199  O   ASN A 284      34.751 -16.634  35.341  1.00 28.11           O
ATOM   3200  CB  ASN A 284      37.063 -18.130  33.915  1.00 21.07           C
ATOM   3201  CG  ASN A 284      38.549 -17.829  33.904  1.00 20.78           C
ATOM   3202  OD1 ASN A 284      39.040 -17.053  34.702  1.00 19.74           O
ATOM   3203  ND2 ASN A 284      39.263 -18.449  32.968  1.00 21.82           N
ATOM   3204  H   ASN A 284      36.875 -20.208  35.174  1.00 23.56           H
ATOM   3205 HD21 ASN A 284      38.893 -18.679  32.227  1.00 26.19           H
ATOM   3206  HA  ASN A 284      36.987 -17.719  35.944  1.00 19.79           H
ATOM   3207  HB2 ASN A 284      36.900 -18.947  33.418  1.00 21.07           H
ATOM   3208  HB3 ASN A 284      36.591 -17.419  33.454  1.00 21.07           H
ATOM   3209 HD22 ASN A 284      38.878 -18.991  32.422  1.00 21.82           H
HETATM 5988  C1  NAG A 467      40.601 -17.959  32.799  1.00 27.22           C
HETATM 5989  C2  NAG A 467      41.289 -19.314  32.714  1.00 22.16           C
HETATM 5990  C3  NAG A 467      42.783 -19.123  32.507  1.00 54.68           C
HETATM 5991  C4  NAG A 467      43.034 -18.265  31.278  1.00 23.55           C
HETATM 5992  C5  NAG A 467      42.261 -16.957  31.391  1.00 33.78           C
HETATM 5993  C6  NAG A 467      42.388 -16.125  30.141  1.00 24.49           C
HETATM 5994  C7  NAG A 467      41.114 -21.444  33.906  1.00 21.47           C
HETATM 5995  C8  NAG A 467      40.844 -22.107  35.214  1.00 20.34           C
HETATM 5996  N2  NAG A 467      41.041 -20.110  33.902  1.00 24.73           N
HETATM 5997  O3  NAG A 467      43.399 -20.391  32.338  1.00 54.77           O
HETATM 5998  O4  NAG A 467      44.417 -17.960  31.155  1.00 53.51           O
HETATM 5999  O5  NAG A 467      40.861 -17.215  31.600  1.00 31.39           O
HETATM 6000  O6  NAG A 467      41.470 -15.043  30.154  1.00 46.51           O
HETATM 6001  O7  NAG A 467      41.392 -22.081  32.897  1.00 22.76           O
HETATM 6002  H1  NAG A 467      40.952 -17.467  33.566  1.00 32.66           H
HETATM 6003  H2  NAG A 467      40.934 -19.790  31.940  1.00 26.59           H
HETATM 6004  H3  NAG A 467      43.163 -18.682  33.290  1.00 65.62           H
HETATM 6005  H4  NAG A 467      42.738 -18.746  30.482  1.00 28.26           H
HETATM 6006  H5  NAG A 467      42.608 -16.449  32.148  1.00 40.54           H
HETATM 6007  H61 NAG A 467      42.210 -16.687  29.363  1.00 29.39           H
HETATM 6008  H62 NAG A 467      43.296 -15.773  30.082  1.00 29.39           H
HETATM 6009  H81 NAG A 467      40.882 -23.076  35.101  1.00 24.41           H
HETATM 6010  H82 NAG A 467      39.958 -21.851  35.532  1.00 24.41           H
HETATM 6011  H83 NAG A 467      41.516 -21.829  35.865  1.00 24.41           H
HETATM 6012  HN2 NAG A 467      40.836 -19.681  34.680  1.00 29.67           H
HETATM 6013  HO3 NAG A 467      42.779 -20.998  32.145  1.00 65.72           H
HETATM 6014  HO4 NAG A 467      44.884 -18.471  31.711  1.00 64.21           H
HETATM 6015  HO6 NAG A 467      40.829 -15.206  30.746  1.00 55.81           H
""")
  result = utils.detect_hydrogen_nomenclature_problem(
    pdb_file="tmp_mmtbx_utils_asn_nag.pdb")
  assert (result.n_asn_hd22 == 1)

def run():
  verbose = "--verbose" in sys.argv[1:]
  exercise_00(verbose=verbose)
  exercise_01(verbose=verbose)
  exercise_02(verbose=verbose)
  exercise_03(verbose=verbose)
  exercise_05(verbose=verbose)
  exercise_06(verbose=verbose)
  exercise_07(verbose=verbose)
  exercise_08(verbose=verbose)
  exercise_09(verbose=verbose)
  exercise_10(verbose=verbose)
  exercise_11(verbose=verbose)
  exercise_12(verbose=verbose)
  exercise_13(verbose=verbose)
  exercise_14(verbose=verbose)
  exercise_15(verbose=verbose)
  exercise_16(verbose=verbose)
  exercise_17(verbose=verbose)
  exercise_18(verbose=verbose)
  exercise_19(verbose=verbose)
  exercise_20(verbose=verbose)
  exercise_21(verbose=verbose)
  exercise_22(verbose=verbose)
  exercise_23(verbose=verbose)
  exercise_24(verbose=verbose)
  exercise_25(verbose=verbose)
  exercise_26(verbose=verbose)
  exercise_27(verbose=verbose)
  exercise_28(verbose=verbose)
  exercise_29(verbose=verbose)
  exercise_30(verbose=verbose)
  exercise_d_data_target_d_atomic_params()
  exercise_d_data_target_d_atomic_params2()
  exercise_get_atom_selections(verbose=verbose)
  exercise_f_000()
  exercise_cmdline_load_pdb_and_data()
  exercise_detect_link_problems()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
