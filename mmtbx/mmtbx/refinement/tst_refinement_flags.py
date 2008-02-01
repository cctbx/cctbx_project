from copy import deepcopy
from cctbx.array_family import flex
import sys, os, math, string
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args
from libtbx.utils import Sorry
from mmtbx.refinement import refinement_flags
from cStringIO import StringIO

expected_result_all = \
  """Refinement flags and selection counts:
  individual_sites       =  True (6 atoms)
  rigid_body             =  True (6 atoms in 3 groups)
  individual_adp         =  True (iso = 6 aniso = 6)
  group_adp              =  True (6 atoms in 3 groups)
  tls                    =  True (6 atoms in 3 groups)
  individual_occupancies =  True (6 atoms)
  group_occupancies      =  True (6 atoms in 3 groups)
  group_anomalous        =  True
"""
expected_result_none_false = \
  """Refinement flags and selection counts:
  individual_sites       = False (0 atoms)
  rigid_body             = False (0 atoms in 0 groups)
  individual_adp         = False (iso = 0 aniso = 0)
  group_adp              = False (0 atoms in 0 groups)
  tls                    = False (0 atoms in 0 groups)
  individual_occupancies = False (0 atoms)
  group_occupancies      = False (0 atoms in 0 groups)
  group_anomalous        = False
"""
expected_result_none_true = \
  """Refinement flags and selection counts:
  individual_sites       =  True (0 atoms)
  rigid_body             =  True (0 atoms in 0 groups)
  individual_adp         =  True (iso = 0 aniso = 0)
  group_adp              =  True (0 atoms in 0 groups)
  tls                    =  True (0 atoms in 0 groups)
  individual_occupancies =  True (0 atoms)
  group_occupancies      =  True (0 atoms in 0 groups)
  group_anomalous        =  True
"""
expected_result_mix = \
  """Refinement flags and selection counts:
  individual_sites       =  True (4 atoms)
  rigid_body             =  True (4 atoms in 2 groups)
  individual_adp         =  True (iso = 4 aniso = 4)
  group_adp              =  True (4 atoms in 2 groups)
  tls                    =  True (4 atoms in 2 groups)
  individual_occupancies =  True (4 atoms)
  group_occupancies      =  True (4 atoms in 2 groups)
  group_anomalous        =  True
"""

def all_defined():
  #                 0 1 2 3 4 5 6 7 8 9
  barr = flex.bool([0,1,1,0,1,0,1,1,1,0])
  iarr = [flex.size_t([1,2]), flex.size_t([4]), flex.size_t([6,7,8])]
  return refinement_flags.manager(
    individual_sites       = True,
    rigid_body             = True,
    individual_adp         = True,
    group_adp              = True,
    tls                    = True,
    individual_occupancies = True,
    group_occupancies      = True,
    group_anomalous        = True,
    sites_individual       = barr,
    sites_rigid_body       = iarr,
    adp_individual_iso     = barr,
    adp_individual_aniso   = barr,
    adp_group              = iarr,
    group_h                = iarr,
    adp_tls                = iarr,
    occupancies_individual = iarr,
    occupancies_group      = iarr)

def all_defined_1():
  #                 0 1 2 3 4 5 6 7 8 9
  barr = flex.bool([1,0,1,0,0,0,1,1,0,1])
  iarr = [flex.size_t([0]), flex.size_t([2]), flex.size_t([6,7]),
          flex.size_t([9])]
  return refinement_flags.manager(
    individual_sites       = True,
    rigid_body             = True,
    individual_adp         = True,
    group_adp              = True,
    tls                    = True,
    individual_occupancies = True,
    group_occupancies      = True,
    group_anomalous        = True,
    sites_individual       = barr,
    sites_rigid_body       = iarr,
    adp_individual_iso     = barr,
    adp_individual_aniso   = barr,
    adp_group              = iarr,
    group_h                = iarr,
    adp_tls                = iarr,
    occupancies_individual = iarr,
    occupancies_group      = iarr)

def compare(rm, expected_result, deep_copy=False, selection=None, show=False):
  assert [deep_copy, selection].count(True) != 2
  if(deep_copy):
    rm = rm.deep_copy()
  elif(selection is not None):
    rm = rm.select(selection = selection)
  out = StringIO()
  rm.show(log = out)
  if(show):
    print "-"*80
    print out.getvalue()
    print expected_result
  assert out.getvalue() == expected_result

def exercise_deepcopy_show_select():
  sel_all       = flex.bool([1,1,1,1,1,1,1,1,1,1])
  sel_none      = flex.bool([0,0,0,0,0,0,0,0,0,0])
  sel_all_true  = flex.bool([0,1,1,0,1,0,1,1,1,0])
  sel_all_false = flex.bool([1,0,0,1,0,1,0,0,0,1])
  sel_mix       = flex.bool([1,1,1,0,0,0,0,1,1,1])
  #
  compare(rm = all_defined(), expected_result = expected_result_all)
  compare(rm = all_defined(), expected_result = expected_result_all,       deep_copy = True)
  compare(rm = all_defined(), expected_result = expected_result_all,       selection = sel_all)
  compare(rm = all_defined(), expected_result = expected_result_none_true, selection = sel_none)
  compare(rm = all_defined(), expected_result = expected_result_all,       selection = sel_all_true)
  compare(rm = all_defined(), expected_result = expected_result_none_true, selection = sel_all_false)
  compare(rm = all_defined(), expected_result = expected_result_mix,       selection = sel_mix)
  #
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, deep_copy = True)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, selection = sel_all)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, selection = sel_none)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, selection = sel_all_true)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, selection = sel_all_false)
  compare(rm = refinement_flags.manager(), expected_result = expected_result_none_false, selection = sel_mix)
  #
  rm = all_defined()
  rm_dc = rm.deep_copy()
  rm_dc.individual_sites       = False
  rm_dc.rigid_body             = False
  rm_dc.individual_adp         = False
  rm_dc.group_adp              = False
  rm_dc.tls                    = False
  rm_dc.individual_occupancies = False
  rm_dc.group_occupancies      = False
  rm_dc.group_anomalous        = False
  rm_dc.sites_individual       = None
  rm_dc.sites_rigid_body       = None
  rm_dc.adp_individual_iso     = None
  rm_dc.adp_individual_aniso   = None
  rm_dc.adp_group              = None
  rm_dc.group_h                = None
  rm_dc.adp_tls                = None
  rm_dc.occupancies_individual = None
  rm_dc.occupancies_group      = None
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all
  out = StringIO()
  rm_dc.show(log = out)
  assert out.getvalue() == expected_result_none_false

def exercise_deepcopy_show_select_compare_arrays():
  rm = all_defined()
  sel_mix = flex.bool([1,1,1,0,0,0,0,1,1,1])
  rm_sel = rm.select(selection = sel_mix)
  assert rm_sel.individual_sites
  assert rm_sel.rigid_body
  assert rm_sel.individual_adp
  assert rm_sel.group_adp
  assert rm_sel.tls
  assert rm_sel.individual_occupancies
  assert rm_sel.group_occupancies
  assert rm_sel.group_anomalous
  barr_sel_mix  = flex.bool([0,1,1,1,1,0])
  iarr_sel_mix  = [flex.size_t([1,2]), flex.size_t([3,4])]
  assert approx_equal(rm_sel.sites_individual      , barr_sel_mix)
  assert approx_equal(rm_sel.sites_rigid_body      , iarr_sel_mix)
  assert approx_equal(rm_sel.adp_individual_iso    , barr_sel_mix)
  assert approx_equal(rm_sel.adp_individual_aniso  , barr_sel_mix)
  assert approx_equal(rm_sel.adp_group             , iarr_sel_mix)
  assert approx_equal(rm_sel.group_h               , iarr_sel_mix)
  assert approx_equal(rm_sel.adp_tls               , iarr_sel_mix)
  assert approx_equal(rm_sel.occupancies_individual, iarr_sel_mix)
  assert approx_equal(rm_sel.occupancies_group     , iarr_sel_mix)

def exercise_inflate():
  rm = all_defined()
  barr = flex.bool([1,0,1,1,0,1])
  iarr = [flex.size_t([10]), flex.size_t([12,13]), flex.size_t([15])]
  rm = rm.inflate(
    sites_individual       = barr,
    sites_rigid_body       = iarr,
    adp_individual_iso     = barr,
    adp_individual_aniso   = barr,
    adp_group              = iarr,
    group_h                = iarr,
    adp_tls                = iarr,
    occupancies_individual = iarr,
    occupancies_group      = iarr)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == \
  """Refinement flags and selection counts:
  individual_sites       =  True (10 atoms)
  rigid_body             =  True (10 atoms in 6 groups)
  individual_adp         =  True (iso = 10 aniso = 10)
  group_adp              =  True (10 atoms in 6 groups)
  tls                    =  True (10 atoms in 6 groups)
  individual_occupancies =  True (10 atoms)
  group_occupancies      =  True (10 atoms in 6 groups)
  group_anomalous        =  True
"""
  barr_result  = flex.bool([0,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1])
  iarr_result  = [flex.size_t([1,2]), flex.size_t([4]), flex.size_t([6,7,8]),
                  flex.size_t([10]), flex.size_t([12,13]), flex.size_t([15])]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.occupancies_individual, iarr_result)
  assert approx_equal(rm.occupancies_group     , iarr_result)

def exercise_add_1a():
  #
  rm = all_defined()
  rm = rm.add(next_to_i_seqs = flex.size_t([]))
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all
  #
  rm = rm.add(
    next_to_i_seqs = flex.size_t([]),
    sites_individual       = True,
    sites_rigid_body       = True,
    adp_individual_iso     = True,
    adp_individual_aniso   = True,
    adp_group              = True,
    group_h                = True,
    adp_tls                = True,
    occupancies_individual = True,
    occupancies_group      = True)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all

def exercise_add_1b():
  rm = all_defined()
  # [[1, 2], [4], [6, 7, 8]] - original
  # [[1, 2], [5], [8, 9, 10]] - after insertion for tuples
  #
  # [0,1,1,0,1,0,1,1,1,0] - original
  # [0,1,1,0,0,1,0,0,1,1,1,0,0] - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs         = flex.size_t([4,8,3]),
    sites_individual       = False,
    sites_rigid_body       = False,
    adp_individual_iso     = False,
    adp_individual_aniso   = False,
    adp_group              = False,
    group_h                = False,
    adp_tls                = False,
    occupancies_individual = False,
    occupancies_group      = False)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all
  barr_result = flex.bool([0,1,1,0,0,1,0,0,1,1,1,0,0])
  iarr_result = [flex.size_t([1,2]), flex.size_t([5]), flex.size_t([8,9,10])]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.occupancies_individual, iarr_result)
  assert approx_equal(rm.occupancies_group     , iarr_result)

def exercise_add_1c():
  rm = all_defined()
  # [[1, 2], [4], [6, 7, 8]] - original
  # [[1, 2], [4], [5], [6], [8, 9, 10, 11]] - after insertion for tuples
  #
  # [0,1,1,0,1,0,1,1,1,0]) - original
  # [0,1,1,0,1,1,1,0,1,1,1,1,0] - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs         = flex.size_t([4,8,3]),
    sites_individual       = True,
    sites_rigid_body       = True,
    adp_individual_iso     = True,
    adp_individual_aniso   = True,
    adp_group              = True,
    group_h                = True,
    adp_tls                = True,
    occupancies_individual = True,
    occupancies_group      = True)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == \
  """Refinement flags and selection counts:
  individual_sites       =  True (9 atoms)
  rigid_body             =  True (9 atoms in 5 groups)
  individual_adp         =  True (iso = 9 aniso = 9)
  group_adp              =  True (9 atoms in 5 groups)
  tls                    =  True (9 atoms in 5 groups)
  individual_occupancies =  True (9 atoms)
  group_occupancies      =  True (9 atoms in 5 groups)
  group_anomalous        =  True
"""
  barr_result = flex.bool([0,1,1,0,1,1,1,0,1,1,1,1,0])
  iarr_result = [flex.size_t([1,2]), flex.size_t([5]), flex.size_t([6]),
                 flex.size_t([8, 9, 10, 11]), flex.size_t([4])]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.occupancies_individual, iarr_result)
  assert approx_equal(rm.occupancies_group     , iarr_result)

def exercise_add_2b():
  rm = all_defined_1()
  # [[0], [2], [6,7], [9]] - original
  # [[0], [3], [8,10], [12]] - after insertion for tuples
  #
  # [1,0,1,0,0,0,1,1,0,1] - original
  # [1,0,0,1,0,0,0,0,1,0,1,0,1,0] - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs         = flex.size_t([0,4,6,9]),
    sites_individual       = False,
    sites_rigid_body       = False,
    adp_individual_iso     = False,
    adp_individual_aniso   = False,
    adp_group              = False,
    group_h                = False,
    adp_tls                = False,
    occupancies_individual = False,
    occupancies_group      = False)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == \
  """Refinement flags and selection counts:
  individual_sites       =  True (5 atoms)
  rigid_body             =  True (5 atoms in 4 groups)
  individual_adp         =  True (iso = 5 aniso = 5)
  group_adp              =  True (5 atoms in 4 groups)
  tls                    =  True (5 atoms in 4 groups)
  individual_occupancies =  True (5 atoms)
  group_occupancies      =  True (5 atoms in 4 groups)
  group_anomalous        =  True
"""
  barr_result = flex.bool([1,0,0,1,0,0,0,0,1,0,1,0,1,0])
  iarr_result = [flex.size_t([0]), flex.size_t([3]), flex.size_t([8,10]),
                 flex.size_t([12])]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.occupancies_individual, iarr_result)
  assert approx_equal(rm.occupancies_group     , iarr_result)

def exercise_add_2c():
  rm = all_defined_1()
  # [[0], [2], [6,7], [9]] - original
  # [[0], [1], [3], [6], [8,9,10], [12], [13]] - after insertion for tuples
  #
  # [1,0,1,0,0,0,1,1,0,1] - original
  # [1,1,0,1,0,0,1,0,1,1,1,0,1,1] - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs         = flex.size_t([0,4,6,9]),
    sites_individual       = True,
    sites_rigid_body       = True,
    adp_individual_iso     = True,
    adp_individual_aniso   = True,
    adp_group              = True,
    group_h                = True,
    adp_tls                = True,
    occupancies_individual = True,
    occupancies_group      = True)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == \
  """Refinement flags and selection counts:
  individual_sites       =  True (9 atoms)
  rigid_body             =  True (9 atoms in 7 groups)
  individual_adp         =  True (iso = 9 aniso = 9)
  group_adp              =  True (9 atoms in 7 groups)
  tls                    =  True (9 atoms in 7 groups)
  individual_occupancies =  True (9 atoms)
  group_occupancies      =  True (9 atoms in 7 groups)
  group_anomalous        =  True
"""
  barr_result = flex.bool([1,1,0,1,0,0,1,0,1,1,1,0,1,1])
  iarr_result = [flex.size_t([0]), flex.size_t([1]), flex.size_t([3]),
                 flex.size_t([8, 9, 10]), flex.size_t([12]), flex.size_t([13]),
                 flex.size_t([6])]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.occupancies_individual, iarr_result)
  assert approx_equal(rm.occupancies_group     , iarr_result)

if(__name__ == "__main__"):
  exercise_deepcopy_show_select()
  exercise_deepcopy_show_select_compare_arrays()
  exercise_inflate()
  exercise_add_1a()
  exercise_add_1b()
  exercise_add_1c()
  exercise_add_2b()
  exercise_add_2c()
