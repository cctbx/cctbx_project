from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.refinement.occupancies import occupancy_selections
from iotbx.cli_parser import run_program
from mmtbx.programs import fmodel
import mmtbx.model
import iotbx.pdb
import iotbx.phil
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times, null_out, Sorry
import libtbx.load_env
from six.moves import cStringIO as StringIO
import os
import sys
from six.moves import zip

def extract_serials(atoms, occ_groups):
  r = []
  # for atom in atoms:
  #   assert atom.serial == atom.i_seq, "%s %d" % (atom.serial, atom.i_seq)
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

def get_model(file_name, log):
  pdb_interpretation_params = iotbx.phil.parse(
          input_string=pdb_interpretation.grand_master_phil_str, process_includes=True).extract()
  pdb_interpretation_params.pdb_interpretation.sort_atoms=False
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      stop_for_unknowns = False,
      log=log)
  model.process(pdb_interpretation_params=pdb_interpretation_params)
  return model

def get_model_str(strings, log):
  pdb_interpretation_params = iotbx.phil.parse(
          input_string=pdb_interpretation.grand_master_phil_str, process_includes=True).extract()
  pdb_interpretation_params.pdb_interpretation.sort_atoms=False
  pdb_inp = iotbx.pdb.input(lines=strings, source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      stop_for_unknowns = False,
      log=log)
  model.process(pdb_interpretation_params=pdb_interpretation_params)
  return model

def exercise_00(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  base = [ [[2],[3]], [[6,7,8,9,10],[11,12,13,14,15]], [[16],[17]], [[24,25,26,27],[28,29,30,31]] ]
  # default
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  target = base[:]
  target.insert(3, [[21]])
  target.insert(4, [[23]])
  assert approx_equal(res, target)
  # default + add water
  res = occupancy_selections(
    model = model,
    add_water         = True,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  base_21_23 = target[:]
  target.extend([[[18]], [[19]], [[20]]])
  assert approx_equal(res, target)
  # 1
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  res = extract_serials(model.get_atoms(), res)
  target = base_21_23[:]
  target.extend([[[0]], [[1]], [[4]], [[5]]])
  assert approx_equal(res, target)
  res = occupancy_selections(
    model = model,
    add_water         = True,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and not (altloc A or altloc B)'])
  res = extract_serials(model.get_atoms(), res)
  target.extend([[[18]], [[19]], [[20]]])
  assert approx_equal(res, target)
  # 2
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['resseq 0 and (name S or name O1)'], ['resseq 0 and (name O3 or name O4)'] ])
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False,
    other_constrained_groups = other_constrained_groups)
  res = extract_serials(model.get_atoms(), res)
  target = base_21_23[:]
  target.extend([[[0, 1]], [[4, 5]]])
  assert approx_equal(res, target)
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['resseq 0 and (name S or name O1)'], ['resseq 0 and (name O3 or name O4)'] ])
  res = occupancy_selections(
    model = model,
    add_water         = True,
    as_flex_arrays    = False,
    other_constrained_groups = other_constrained_groups)
  res = extract_serials(model.get_atoms(), res)
  target.extend([[[18]], [[19]], [[20]]])
  assert approx_equal(res, target)
  # 3
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['resseq 0 and (name O3 or name O4)'] ])
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False,
    other_individual_selection_strings = ['resseq 0 and (name S or name O1)'],
    other_constrained_groups = other_constrained_groups)
  res = extract_serials(model.get_atoms(), res)
  target = base_21_23[:]
  target.extend([[[0]], [[1]], [[4, 5]]])
  assert approx_equal(res, target)

def exercise_01(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_h.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [ [[0,1,2,3,4,10,12,14,16,18,20,22], [5,6,7,8,9,11,13,15,17,19,21,23]] ]
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_02(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/occ_mix1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [ [[0,1,2,3,4,5,6,7,8,9,10,11,12,13], [14,15,16,17,18,19,20,21,22,23,24,25,26,27]] ]
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_03(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_hd.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [ [[7]], [[8]], [[9],[12]], [[10],[13]], [[11],[14]] ]
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_05(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ala_lys_arg_ser_tyr_neutron_hd.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [ [[9],[12]],  [[10],[13]], [[11],[14]], [[33],[37]], [[34],[38]],
           [[35],[39]], [[36],[40]], [[59],[65]], [[60],[66]], [[61],[67]],
           [[62],[68]], [[63],[69]], [[64],[70]], [[80],[82]], [[81],[83]],
           [[103],[105]], [[104],[106]]]
  res = occupancy_selections(
    model = model,
    as_flex_arrays      = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_06(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/NAD_594_HD.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [ [[62]], [[113]], [[65],[77]],  [[66],[78]],  [[67],[79]], [[68],[80]],
                            [[69],[81]],  [[70],[82]],  [[71],[83]], [[72],[84]],
                            [[73],[85]],  [[74],[86]],  [[75],[87]], [[76],[88]],
                            [[124],[127]],[[125],[128]],[[126],[129]]]
  res = occupancy_selections(
    model = model,
    as_flex_arrays      = False)
  assert approx_equal(res, base)

def exercise_07(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[0, 1, 2, 3, 4]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['resseq 0'] ])
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_08(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
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
    result = occupancy_selections(
      model = model,
      other_constrained_groups = other_constrained_groups,
      as_flex_arrays    = False)
    assert approx_equal(result, answer)

def exercise_09(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
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
    result = occupancy_selections(
      model = model,
      other_individual_selection_strings = [individual_selection],
      as_flex_arrays    = False)
    assert approx_equal(result, answer)

def exercise_10(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  e = None
  try:
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [ ['resseq 0'] ])
    result = occupancy_selections(
      model = model,
      other_constrained_groups = other_constrained_groups,
      other_individual_selection_strings = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception as e:
    assert str(e) == "Duplicate selection: same atoms selected for individual and group occupancy refinement."

def exercise_11(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  e = None
  try:
    result = occupancy_selections(
      model = model,
      remove_selection = ['resseq 0'],
      other_individual_selection_strings = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception as e:
    assert str(e) == "Duplicate selection: occupancies of same atoms selected to be fixed and to be refined."
  e = None
  try:
    other_constrained_groups = make_up_other_constrained_groups_obj(
      selections = [ ['resseq 0'] ])
    result = occupancy_selections(
      model = model,
      other_constrained_groups = other_constrained_groups,
      remove_selection = ['resseq 0'],
      as_flex_arrays    = False)
  except Exception as e:
    assert str(e) == "Duplicate selection: occupancies of same atoms selected to be fixed and to be refined."

def exercise_12(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[4],[5]], [[16],[17]], [[21]], [[23,24,25,26,27,28,29,30,31]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['resseq 6'] ])
  result = occupancy_selections(
    model = model,
    remove_selection = ['resseq 1'],
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)
  #
  answer = [ [[4],[5]], [[16],[17]], [[21]], [[23]], [[24]], [[25]], [[26]], [[27]], [[28]], [[29]], [[30]], [[31]] ]
  result = occupancy_selections(
    model = model,
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
  model = get_model(pdb_file, log)
  answer = [ [[8],[9]], [[10]], [[0],[1]], [[2],[3]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and name N','chain A and resseq 1 and name CA'],
                   ['chain A and resseq 1 and name C','chain A and resseq 1 and name O'] ]
    )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_14(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8],[9]], [[10]], [[0,1,2],[3,4]], [[5],[6]], [[7]] ]

  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and (name N or name CA or name C)', 'chain A and resseq 1 and (name O or name CB)'],
                   ['chain A and resseq 1 and name CG','chain A and resseq 1 and name CD'],
                   ['chain A and resseq 1 and name CE'] ]
    )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_15(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8],[9]], [[0,1,2],[10]], [[5,7]] ]

  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [ ['chain A and resseq 1 and (name N or name CA or name C)', 'chain S and resseq 1'],
                   ['chain A and resseq 1 and name CG or chain A and resseq 1 and name CE'] ]
    )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_16(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8],[9],[10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [
      ['chain A and resseq 1 and name NZ and altloc A', 'chain A and resseq 1 and name NZ and altloc B', 'chain S and resseq 1'] ]
    )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_17(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8,9,10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [
      ['chain A and resseq 1 and name NZ and altloc A or chain A and resseq 1 and name NZ and altloc B or chain S and resseq 1'] ]
    )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_18(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_2.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8],[9],[10]] ]
  other_constrained_groups = make_up_other_constrained_groups_obj(
   selections = [
    ['chain A and resseq 1 and name NZ and altloc A','chain A and resseq 1 and name NZ and altloc B','chain S and resseq 1 and altloc C']]
  )
  result = occupancy_selections(
    model = model,
    other_constrained_groups = other_constrained_groups,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_19(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/lys_1.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[8],[9],[10]] ]
  tmp = "chain A and resseq 1 and name XX and altloc A"
  other_constrained_groups = make_up_other_constrained_groups_obj(
    selections = [[
      tmp,
      'chain A and resseq 1 and name NZ and altloc B',
      'chain S and resseq 1']])
  try:
    result = occupancy_selections(
      model = model,
      other_constrained_groups = other_constrained_groups,
      as_flex_arrays    = False)
  except Exception as e:
    assert str(e) == \
      'Selection string results in empty selection (selects no atoms): "%s"' \
      % tmp

def exercise_20(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/ile_2conf_h.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  answer = [ [[4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], [19,20,21,22,23,24,25,26,27,28,29,30,31,32,33]] ]
  result = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  assert approx_equal(result, answer)

def exercise_21(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_3.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
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
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_22(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_4.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [[[0, 1, 2, 3, 8, 9, 10, 11, 12], [4, 5, 6, 7, 13, 14, 15, 16, 17]]]
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
  assert approx_equal(res, base)

def exercise_23(verbose):
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/gocr_5.pdb",
    test=os.path.isfile)
  if (verbose): log = sys.stdout
  else: log = StringIO()
  model = get_model(pdb_file, log)
  #
  base = [[[1, 2, 3, 4, 5, 6]], [[7, 8, 9, 10, 11], [12, 13, 14, 15, 16]]]
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  res = extract_serials(model.get_atoms(), res)
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
    model = get_model_str(pdb_str, log)
    res = occupancy_selections(
      model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
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
  model = get_model_str(pdb_str, log)
  res = occupancy_selections(
    model = model,
    as_flex_arrays    = False)
  answer = [ [[0]],
             [[4, 5, 6, 7, 8, 21],
              [9, 10, 11, 12, 13, 31]],
             [[14, 15, 16, 17, 18, 19, 20, 22, 23],
              [24, 25, 26, 27, 28, 29, 30, 32, 33]] ]
  assert approx_equal(res, answer)

def prepare_correlated_occupancy_inputs(
    prefix="tst_group_correlated_occupancy",
    create_mtz=False,
    d_min=1.0):
  pdb_raw = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.056   4.638   6.050  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.058   4.194   4.668  1.00 16.57           C
ATOM      3  C   GLY A   1      -7.993   3.144   4.430  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.521   2.511   5.374  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.616   2.953   3.169  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.526   2.044   2.840  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.216   2.527   3.434  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.943   3.727   3.466  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.382   1.888   1.330  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.632   1.344   0.685  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.042   0.216   0.957  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.247   2.142  -0.178  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.405   1.583   3.898  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.172   1.915   4.595  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.922   1.362   3.915  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.816   0.158   3.672  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.243   1.409   6.039  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.000   1.749   6.841  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.705   2.920   7.082  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.272   0.724   7.270  1.00 13.48           N
ATOM     21  N   GLN A   4      -0.987   2.256   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.361   1.860   3.201  1.00 10.53           C
ATOM     23  C   GLN A   4       1.398   2.605   4.031  1.00 10.24           C
ATOM     24  O   GLN A   4       1.454   3.834   4.025  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.626   2.117   1.712  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.924   1.459   1.221  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.465   2.050  -0.073  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.674   3.260  -0.178  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.708   1.192  -1.059  1.00  9.05           N
ATOM     30  N  AGLN A   5       2.202   1.848   4.775  0.62 10.38           N
ATOM     31  CA AGLN A   5       3.288   2.419   5.569  0.62 11.39           C
ATOM     32  C  AGLN A   5       4.638   1.844   5.123  0.62 11.52           C
ATOM     33  O  AGLN A   5       4.824   0.625   5.095  0.62 12.05           O
ATOM     34  CB AGLN A   5       3.046   2.170   7.063  0.62 11.96           C
ATOM     35  CG AGLN A   5       1.854   2.946   7.622  0.62 10.81           C
ATOM     36  CD AGLN A   5       1.361   2.406   8.951  0.62 13.10           C
ATOM     37  OE1AGLN A   5       0.800   1.312   9.019  0.62 10.65           O
ATOM     38  NE2AGLN A   5       1.562   3.175  10.016  0.62 12.30           N
ATOM     39  N  BGLN A   5       2.239   1.858   4.725  0.38 10.38           N
ATOM     40  CA BGLN A   5       3.326   2.476   5.450  0.38 11.39           C
ATOM     41  C  BGLN A   5       4.639   1.850   5.057  0.38 11.52           C
ATOM     42  O  BGLN A   5       4.814   0.627   5.020  0.38 12.05           O
ATOM     43  CB BGLN A   5       3.110   2.331   6.919  0.38 11.96           C
ATOM     44  CG BGLN A   5       2.695   0.980   7.141  0.38 10.81           C
ATOM     45  CD BGLN A   5       2.882   0.618   8.479  0.38 13.10           C
ATOM     46  OE1BGLN A   5       2.538   1.369   9.406  0.38 10.65           O
ATOM     47  NE2BGLN A   5       3.380  -0.597   8.664  0.38 12.30           N
ATOM     48  N   ASN A   6       5.565   2.732   4.753  1.00 11.99           N
ATOM     49  CA  ASN A   6       6.868   2.339   4.280  1.00 12.30           C
ATOM     50  C   ASN A   6       7.881   2.785   5.302  1.00 13.40           C
ATOM     51  O   ASN A   6       8.262   3.954   5.351  1.00 13.92           O
ATOM     52  CB  ASN A   6       7.133   2.954   2.915  1.00 12.13           C
ATOM     53  CG  ASN A   6       5.988   2.721   1.955  1.00 12.77           C
ATOM     54  OD1 ASN A   6       5.795   1.608   1.466  1.00 14.27           O
ATOM     55  ND2 ASN A   6       5.211   3.764   1.690  1.00 10.07           N
ATOM     56  N  ATYR A   7       8.304   1.849   6.146  0.59 14.70           N
ATOM     57  CA ATYR A   7       9.167   2.166   7.280  0.59 15.18           C
ATOM     58  C  ATYR A   7      10.622   2.326   6.868  0.59 15.91           C
ATOM     59  O  ATYR A   7      11.054   1.799   5.844  0.59 15.76           O
ATOM     60  CB ATYR A   7       9.044   1.086   8.356  0.59 15.35           C
ATOM     61  CG ATYR A   7       7.640   0.946   8.887  0.59 14.45           C
ATOM     62  CD1ATYR A   7       6.759   0.027   8.335  0.59 15.68           C
ATOM     63  CD2ATYR A   7       7.187   1.750   9.924  0.59 14.80           C
ATOM     64  CE1ATYR A   7       5.469  -0.098   8.810  0.59 13.46           C
ATOM     65  CE2ATYR A   7       5.899   1.633  10.407  0.59 14.33           C
ATOM     66  CZ ATYR A   7       5.044   0.707   9.845  0.59 15.09           C
ATOM     67  OH ATYR A   7       3.759   0.583  10.319  0.59 14.39           O
ATOM     68  OXTATYR A   7      11.394   2.990   7.558  0.59 17.49           O
ATOM     70  N  BTYR A   7       8.323   1.843   6.116  0.41 14.70           N
ATOM     71  CA BTYR A   7       9.149   2.183   7.247  0.41 15.18           C
ATOM     72  C  BTYR A   7      10.629   2.316   6.861  0.41 15.91           C
ATOM     73  O  BTYR A   7      11.084   1.756   5.864  0.41 15.76           O
ATOM     74  CB BTYR A   7       8.954   1.147   8.348  0.41 15.35           C
ATOM     75  CG BTYR A   7       9.942   1.356   9.417  0.41 14.45           C
ATOM     76  CD1BTYR A   7       9.807   2.381  10.320  0.41 15.68           C
ATOM     77  CD2BTYR A   7      11.054   0.580   9.473  0.41 14.80           C
ATOM     78  CE1BTYR A   7      10.746   2.569  11.248  0.41 13.46           C
ATOM     79  CE2BTYR A   7      11.968   0.749  10.405  0.41 14.33           C
ATOM     80  CZ BTYR A   7      11.858   1.724  11.252  0.41 15.09           C
ATOM     81  OH BTYR A   7      12.921   1.747  12.113  0.41 14.39           O
ATOM     82  OXTBTYR A   7      11.408   3.001   7.529  0.41 17.49           O
TER
HETATM   83  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   84  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   85  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   86  O  AHOH A  11      11.808   4.179   9.970  0.60 23.99           O
HETATM   87  O  AHOH A  12      13.605   1.327   9.198  0.60 26.17           O
HETATM   88  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   89  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
TER
"""
  pdb_in = "%s_in.pdb" % prefix
  with open(pdb_in, "w") as f:
    f.write(pdb_raw)
  if (create_mtz):
    args = [
      pdb_in,
      "high_resolution=%g" % d_min,
      "type=real",
      "label=F",
      "add_sigmas=True",
      "r_free_flags_fraction=0.1",
      "random_seed=12345",
      "output.file_name=%s.mtz" % prefix,
    ]
    run_program(program_class=fmodel.Program, args=args, logger=null_out())
  pdb_file = iotbx.pdb.input(pdb_in)
  hierarchy = pdb_file.construct_hierarchy()
  xrs = pdb_file.xray_structure_simple()
  for atom in hierarchy.atoms():
    atom.b = 5
    if (atom.occ < 1.0):
      atom.occ = 0.5
  with open("%s_start.pdb" % prefix, "w") as f:
    f.write(hierarchy.as_pdb_string(crystal_symmetry=xrs))

def exercise_regroup_3d(verbose):
  if (verbose): log = sys.stdout
  else: log = StringIO()
  prepare_correlated_occupancy_inputs()
  # File #1 (with homogenized occupancies) should work
  # File #2 should fail due to inconsistent occupancies
  pdb_files = [
    "tst_group_correlated_occupancy_start.pdb",
    "tst_group_correlated_occupancy_in.pdb",
  ]
  #
  constraint_groups = occupancy_selections(
        model = get_model("tst_group_correlated_occupancy_start.pdb", log),
        constrain_correlated_3d_groups=True,
        log=null_out())
  for g in constraint_groups:
     assert [ [ggg for ggg in gg] for gg in g ] == \
       [[29, 30, 31, 32, 33, 34, 35, 36, 37, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 84, 85],
        [38, 39, 40, 41, 42, 43, 44, 45, 46, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80]]
  #
  sorry_found = False
  try:
    constraint_groups = occupancy_selections(
      model = get_model("tst_group_correlated_occupancy_in.pdb", log),
      constrain_correlated_3d_groups=True,
      log=null_out())
  except Sorry as s:
    assert ("Inconsistent occupancies" in str(s)), str(s)
    sorry_found = True
  assert sorry_found

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
  exercise_regroup_3d(verbose=verbose)
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()
