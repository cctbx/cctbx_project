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
  processed_pdb_files_srv = utils.process_pdb_file_srv(log = log)
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
    answer = [ [[0,1,2,3,4,5,6,7,8,9], 
                [10,11,12,13,14,15,16,17,18,19], 
                [20,21,22,23,24,25,26,27,28,29,67], 
                [30,31,32,33,34,35,36,37,38,39,68]], 
               [[45,46,47,48,49,50,51,52,53,54,55], 
                [56,57,58,59,60,61,62,63,64,65,66]]]
    
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
  n_bad = utils.detect_asparagine_link_problem("tmp_mmtbx_utils_asn_nag.pdb")
  assert (n_bad == 1)

def run():
  verbose = "--verbose" in sys.argv[1:]
  #exercise_00(verbose=verbose)
  #exercise_01(verbose=verbose)
  #exercise_02(verbose=verbose)
  #exercise_03(verbose=verbose)
  #exercise_05(verbose=verbose)
  #exercise_06(verbose=verbose)
  #exercise_07(verbose=verbose)
  #exercise_08(verbose=verbose)
  #exercise_09(verbose=verbose)
  #exercise_10(verbose=verbose)
  #exercise_11(verbose=verbose)
  #exercise_12(verbose=verbose)
  #exercise_13(verbose=verbose)
  #exercise_14(verbose=verbose)
  #exercise_15(verbose=verbose)
  #exercise_16(verbose=verbose)
  #exercise_17(verbose=verbose)
  #exercise_18(verbose=verbose)
  #exercise_19(verbose=verbose)
  #exercise_20(verbose=verbose)
  #exercise_21(verbose=verbose)
  #exercise_22(verbose=verbose)
  exercise_23(verbose=verbose)
  exercise_24(verbose=verbose)
  exercise_25(verbose=verbose)
  #exercise_d_data_target_d_atomic_params()
  #exercise_get_atom_selections(verbose=verbose)
  #exercise_f_000()
  #exercise_cmdline_load_pdb_and_data()
  #exercise_detect_link_problems()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
