import sys, os
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import libtbx.load_env
from cStringIO import StringIO
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
    fmodel.update_all_scales()
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
  exercise_d_data_target_d_atomic_params()
  print format_cpu_times()

if (__name__ == "__main__"):
  run()
