from cctbx.array_family import flex
from libtbx.test_utils import approx_equal, show_diff
from libtbx.utils import format_cpu_times
from mmtbx.refinement import refinement_flags
from cStringIO import StringIO

expected_result_all = \
  """Refinement flags and selection counts:
  individual_sites       = %s (9 atoms)
  torsion_angles         = %s (9 atoms)
  rigid_body             = %s (9 atoms in 6 groups)
  individual_adp         = %s (iso = 9 aniso = 9)
  group_adp              = %s (9 atoms in 6 groups)
  tls                    = %s (9 atoms in 6 groups)
  occupancies            = %s (9 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8)
expected_result_none_false = \
  """Refinement flags and selection counts:
  individual_sites       = %s (0 atoms)
  torsion_angles         = %s (0 atoms)
  rigid_body             = %s (0 atoms in 0 groups)
  individual_adp         = %s (iso = 0 aniso = 0)
  group_adp              = %s (0 atoms in 0 groups)
  tls                    = %s (0 atoms in 0 groups)
  occupancies            = %s (0 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(False)]*8)
expected_result_none_true = \
  """Refinement flags and selection counts:
  individual_sites       = %s (0 atoms)
  torsion_angles         = %s (0 atoms)
  rigid_body             = %s (0 atoms in 0 groups)
  individual_adp         = %s (iso = 0 aniso = 0)
  group_adp              = %s (0 atoms in 0 groups)
  tls                    = %s (0 atoms in 0 groups)
  occupancies            = %s (0 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8)
expected_result_mix = \
  """Refinement flags and selection counts:
  individual_sites       = %s (4 atoms)
  torsion_angles         = %s (4 atoms)
  rigid_body             = %s (4 atoms in 3 groups)
  individual_adp         = %s (iso = 4 aniso = 4)
  group_adp              = %s (4 atoms in 3 groups)
  tls                    = %s (4 atoms in 3 groups)
  occupancies            = %s (4 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8)

def all_defined():
  #                 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  barr = flex.int([0,1,1,0,1,1,0,1,1,0, 1, 1, 0, 1, 0, 0]).as_bool()
  iarr = [ flex.size_t([1]),
           flex.size_t([2]),
           flex.size_t([4,5]),
           flex.size_t([7,8]),
           flex.size_t([10]),
           flex.size_t([11,13]) ]
  oarr = [ [flex.size_t([1]),flex.size_t([2])],
           [flex.size_t([4,5]),flex.size_t([7,8])],
           [flex.size_t([10])],
           [flex.size_t([11,13])] ]
  return refinement_flags.manager(
    individual_sites     = True,
    torsion_angles       = True,
    rigid_body           = True,
    individual_adp       = True,
    group_adp            = True,
    tls                  = True,
    occupancies          = True,
    group_anomalous      = True,
    sites_individual     = barr,
    sites_torsion_angles = barr,
    sites_rigid_body     = iarr,
    adp_individual_iso   = barr,
    adp_individual_aniso = barr,
    adp_group            = iarr,
    group_h              = iarr,
    adp_tls              = iarr,
    s_occupancies        = oarr)

def all_defined_1():
  #                 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  barr = flex.int([1,1,1,0,1,1,0,0,1,0, 0, 0, 1, 1, 1, 1]).as_bool()
  iarr = [ flex.size_t([0]),
           flex.size_t([1,2]),
           flex.size_t([4,5]),
           flex.size_t([8]),
           flex.size_t([12,13,14]),
           flex.size_t([15]) ]
  oarr = [ [flex.size_t([0]),flex.size_t([1,2])],
           [flex.size_t([4,5])],
           [flex.size_t([8])],
           [flex.size_t([12,13,14]),flex.size_t([15])] ]
  return refinement_flags.manager(
    individual_sites     = True,
    torsion_angles       = True,
    rigid_body           = True,
    individual_adp       = True,
    group_adp            = True,
    tls                  = True,
    occupancies          = True,
    group_anomalous      = True,
    sites_individual     = barr,
    sites_torsion_angles = barr,
    sites_rigid_body     = iarr,
    adp_individual_iso   = barr,
    adp_individual_aniso = barr,
    adp_group            = iarr,
    group_h              = iarr,
    adp_tls              = iarr,
    s_occupancies        = oarr)

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
  assert not show_diff(out.getvalue(), expected_result)

def exercise_deepcopy_show_select():
  sel_all       = flex.int([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]).as_bool()
  sel_none      = flex.int([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]).as_bool()
  sel_all_true  = flex.int([0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,0]).as_bool()
  sel_all_false = flex.int([1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,1]).as_bool()
  sel_mix       = flex.int([1,1,0,0,1,1,0,0,0,1,1,0,0,0,1,1]).as_bool()
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
  rm_dc.individual_sites     = False
  rm_dc.torsion_angles       = False
  rm_dc.rigid_body           = False
  rm_dc.individual_adp       = False
  rm_dc.group_adp            = False
  rm_dc.tls                  = False
  rm_dc.occupancies          = False
  rm_dc.group_anomalous      = False
  rm_dc.sites_individual     = None
  rm_dc.sites_torsion_angles = None
  rm_dc.sites_rigid_body     = None
  rm_dc.adp_individual_iso   = None
  rm_dc.adp_individual_aniso = None
  rm_dc.adp_group            = None
  rm_dc.group_h              = None
  rm_dc.adp_tls              = None
  rm_dc.s_occupancies        = None
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all
  out = StringIO()
  rm_dc.show(log = out)
  assert out.getvalue() == expected_result_none_false

def exercise_deepcopy_show_select_compare_arrays():
  rm = all_defined()
  sel_mix = flex.int([1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1]).as_bool()
  rm_sel = rm.select(selection = sel_mix)
  assert rm_sel.individual_sites
  assert rm_sel.torsion_angles
  assert rm_sel.rigid_body
  assert rm_sel.individual_adp
  assert rm_sel.group_adp
  assert rm_sel.tls
  assert rm_sel.occupancies
  assert rm_sel.group_anomalous
  barr_sel_mix = flex.int([0,1,1,1,1,0,1,  0, 0]).as_bool()
  iarr_sel_mix = [ flex.size_t([1]),
                   flex.size_t([2,3]),
                   flex.size_t([4]),
                   flex.size_t([6])]
  oarr_sel_mix = [ [flex.size_t([1])],
                   [flex.size_t([2,3]),flex.size_t([4])],
                   [flex.size_t([6])]]
  assert approx_equal(rm_sel.sites_individual     , barr_sel_mix)
  assert approx_equal(rm_sel.sites_torsion_angles , barr_sel_mix)
  assert approx_equal(rm_sel.sites_rigid_body     , iarr_sel_mix)
  assert approx_equal(rm_sel.adp_individual_iso   , barr_sel_mix)
  assert approx_equal(rm_sel.adp_individual_aniso , barr_sel_mix)
  assert approx_equal(rm_sel.adp_group            , iarr_sel_mix)
  assert approx_equal(rm_sel.group_h              , iarr_sel_mix)
  assert approx_equal(rm_sel.adp_tls              , iarr_sel_mix)
  assert approx_equal(rm_sel.s_occupancies        , oarr_sel_mix)

def exercise_inflate():
  rm = all_defined()
  barr = flex.int([1,0,1,1,0,1]).as_bool()
  iarr = [flex.size_t([16]), flex.size_t([18,19]), flex.size_t([21])]
  oarr = [ [flex.size_t([16]), flex.size_t([18,19])], [flex.size_t([21])]]
  rm = rm.inflate(
    sites_individual     = barr,
    sites_torsion_angles = barr,
    sites_rigid_body     = iarr,
    adp_individual_iso   = barr,
    adp_individual_aniso = barr,
    adp_group            = iarr,
    group_h              = iarr,
    adp_tls              = iarr,
    s_occupancies        = oarr)
  out = StringIO()
  rm.show(log = out)
  assert not show_diff(out.getvalue(), """\
Refinement flags and selection counts:
  individual_sites       = %s (13 atoms)
  torsion_angles         = %s (13 atoms)
  rigid_body             = %s (13 atoms in 9 groups)
  individual_adp         = %s (iso = 13 aniso = 13)
  group_adp              = %s (13 atoms in 9 groups)
  tls                    = %s (13 atoms in 9 groups)
  occupancies            = %s (13 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8))
  #
  barr_result = flex.int([0,1,1,0,1,1,0,1,1,0,1,1,0,1,0,0,  1,0,1,1,0,1]) \
    .as_bool()
  iarr_result = [
    flex.size_t([1]),
    flex.size_t([2]),
    flex.size_t([4,5]),
    flex.size_t([7,8]),
    flex.size_t([10]),
    flex.size_t([11,13]),
       flex.size_t([16]),flex.size_t([18,19]),flex.size_t([21])]
  oarr_result = [
    [flex.size_t([1]),flex.size_t([2])],
    [flex.size_t([4,5]),flex.size_t([7,8])],
    [flex.size_t([10])],
    [flex.size_t([11,13])],
       [flex.size_t([16]), flex.size_t([18,19])], [flex.size_t([21])] ]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_torsion_angles  , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  assert approx_equal(rm.s_occupancies         , oarr_result)

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
    sites_individual     = True,
    sites_torsion_angles = True,
    sites_rigid_body     = True,
    adp_individual_iso   = True,
    adp_individual_aniso = True,
    adp_group            = True,
    group_h              = True,
    adp_tls              = True,
    s_occupancies        = True)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all

def exercise_add_1b():
  rm = all_defined()
  # [ [1], [2], [4,5], [7,8], [10], [11,13] ] - original
  # [ [1], [2], [5,7], [9,10], [13], [14,16] ] - after insertion for tuples
  #
  # [0,1,1,0,1,1,0,1,1,0, 1, 1, 0, 1, 0, 0] - original
  # [0,1,1,0,0,1,0,1,0,1,1,0,0, 1, 1, 0, 1, 0, 0] - after insertion for tuples
  #
  # [ [[1],[2]], [[4,5],[7,8]], [[10]], [[11,13]] ] - original
  # [ [[1],[2]], [[5,7],[9,10]], [[13]], [[14,16]] ] - after insertion for tuples
  rm = rm.add(
    next_to_i_seqs       = flex.size_t([4,8,3]),
    sites_individual     = False,
    sites_torsion_angles = False,
    sites_rigid_body     = False,
    adp_individual_iso   = False,
    adp_individual_aniso = False,
    adp_group            = False,
    group_h              = False,
    adp_tls              = False,
    s_occupancies        = False)
  out = StringIO()
  rm.show(log = out)
  assert out.getvalue() == expected_result_all
  barr_result = flex.int([0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,0,1,0,0]).as_bool()
  iarr_result = [flex.size_t([1]), flex.size_t([2]), flex.size_t([5,7]),
                 flex.size_t([9,10]), flex.size_t([13]), flex.size_t([14,16])]
  oarr_result = [[flex.size_t([1]), flex.size_t([2])], [flex.size_t([5,7]),
                 flex.size_t([9,10])], [flex.size_t([13])], [flex.size_t([14,16])]]
  assert approx_equal(rm.sites_individual     , barr_result)
  assert approx_equal(rm.sites_torsion_angles , barr_result)
  assert approx_equal(rm.sites_rigid_body     , iarr_result)
  assert approx_equal(rm.adp_individual_iso   , barr_result)
  assert approx_equal(rm.adp_individual_aniso , barr_result)
  assert approx_equal(rm.adp_group            , iarr_result)
  assert approx_equal(rm.group_h              , iarr_result)
  assert approx_equal(rm.adp_tls              , iarr_result)
  compare_selections(rm.s_occupancies, oarr_result)

def exercise_add_1c():
  rm = all_defined()
  # [ [1], [2], [4,5], [7,8], [10], [11,13] ] - original
  # [ [1], [2], [4], [5,6,7], [9,10,11], [13], [14,16] ] - after insertion for tuples
  #
  # [0,1,1,0,1,1,0,1,1,0, 1, 1, 0, 1, 0, 0] - original
  # [0,1,1,0,1,1,1,1,0,1,1,1,0, 1, 1, 0, 1, 0, 0] - after insertion for tuples
  #
  # [ [[1],[2]], [[4,5],[7,8]],  [[10]], [[11,13]] ] - original
  #
  # [ [[1],[2]], [[5,7],[9,10]], [[13]], [[14,16]], [[4]], [[6]], [[11]] ] - after insertion for tuples
  rm = rm.add(
    next_to_i_seqs       = flex.size_t([4,8,3]),
    sites_individual     = True,
    sites_torsion_angles = True,
    sites_rigid_body     = True,
    adp_individual_iso   = True,
    adp_individual_aniso = True,
    adp_group            = True,
    group_h              = True,
    adp_tls              = True,
    s_occupancies        = True)
  out = StringIO()
  rm.show(log = out)
  assert not show_diff(out.getvalue(), """\
Refinement flags and selection counts:
  individual_sites       = %s (12 atoms)
  torsion_angles         = %s (12 atoms)
  rigid_body             = %s (12 atoms in 7 groups)
  individual_adp         = %s (iso = 12 aniso = 12)
  group_adp              = %s (12 atoms in 7 groups)
  tls                    = %s (12 atoms in 7 groups)
  occupancies            = %s (12 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8))
  barr_result = flex.int([0,1,1,0,1,1,1,1,0,1,1,1,0, 1, 1, 0, 1, 0, 0]) \
    .as_bool()
  iarr_result = [ flex.size_t([1]), flex.size_t([2]),
    flex.size_t([5,6,7]), flex.size_t([9,10,11]), flex.size_t([13]),
    flex.size_t([14,16]), flex.size_t([4]) ]
  oarr_result = [ [flex.size_t([1]),flex.size_t([2])], [flex.size_t([5,7]),
    flex.size_t([9,10])], [flex.size_t([13])], [flex.size_t([14,16])],
    [flex.size_t([4])], [flex.size_t([6])], [flex.size_t([11])] ]
  assert approx_equal(rm.sites_individual    , barr_result)
  assert approx_equal(rm.sites_torsion_angles, barr_result)
  assert approx_equal(rm.sites_rigid_body    , iarr_result)
  assert approx_equal(rm.adp_individual_iso  , barr_result)
  assert approx_equal(rm.adp_individual_aniso, barr_result)
  assert approx_equal(rm.adp_group           , iarr_result)
  assert approx_equal(rm.group_h             , iarr_result)
  assert approx_equal(rm.adp_tls             , iarr_result)
  compare_selections(rm.s_occupancies, oarr_result)

def exercise_add_2b():
  rm = all_defined_1()
  # [ [0], [1,2], [4,5], [8], [12,13,14], [15] ] - original
  # [ [0], [2,3], [5,6], [9], [14,15,16], [18] ] - after insertion for tuples
  #
  # [ [[0],[1,2]], [[4,5]], [[8]], [[12,13,14],[15]] ] - original
  # [ [[0],[2,3]], [[5,6]], [[9]], [[14,15,16],[18]] ] - after insertion for tuples
  #
  #  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  # [1,1,1,0,1,1,0,0,1,0, 0, 0, 1, 1, 1, 1] - original
  # [1,0,1,1,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,0]  - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs       = flex.size_t([0,10,14,15]),
    sites_individual     = False,
    sites_torsion_angles = False,
    sites_rigid_body     = False,
    adp_individual_iso   = False,
    adp_individual_aniso = False,
    adp_group            = False,
    group_h              = False,
    adp_tls              = False,
    s_occupancies        = False)
  out = StringIO()
  rm.show(log = out)
  assert not show_diff(out.getvalue(), """\
Refinement flags and selection counts:
  individual_sites       = %s (10 atoms)
  torsion_angles         = %s (10 atoms)
  rigid_body             = %s (10 atoms in 6 groups)
  individual_adp         = %s (iso = 10 aniso = 10)
  group_adp              = %s (10 atoms in 6 groups)
  tls                    = %s (10 atoms in 6 groups)
  occupancies            = %s (10 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8))
  barr_result = flex.int([1,0,1,1,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,0]).as_bool()
  iarr_result = [ flex.size_t([0]), flex.size_t([2,3]), flex.size_t([5,6]),
    flex.size_t([9]), flex.size_t([14,15,16]), flex.size_t([18]) ]
  oarr_result = [ [flex.size_t([0]),flex.size_t([2,3])], [flex.size_t([5,6])],
    [flex.size_t([9])], [flex.size_t([14,15,16]),flex.size_t([18])] ]
  assert approx_equal(rm.sites_individual    , barr_result)
  assert approx_equal(rm.sites_torsion_angles, barr_result)
  assert approx_equal(rm.sites_rigid_body    , iarr_result)
  assert approx_equal(rm.adp_individual_iso  , barr_result)
  assert approx_equal(rm.adp_individual_aniso, barr_result)
  assert approx_equal(rm.adp_group           , iarr_result)
  assert approx_equal(rm.group_h             , iarr_result)
  assert approx_equal(rm.adp_tls             , iarr_result)
  compare_selections(rm.s_occupancies, oarr_result)

def exercise_add_2c():
  rm = all_defined_1()
  # [ [0], [1,2], [4,5], [8], [12,13,14], [15] ] - original
  # [ [0], [1], [2,3], [5,6], [9], [12], [14,15,16,17], [18], [19] ]  - after insertion for tuples
  #
  # [ [[0],[1,2]], [[4,5]], [[8]], [[12,13,14],[15]] ] - original
  # [ [[0],[2,3]], [[5,6]], [[9]], [[14,15,16],[18]], [[1]], [[12]], [[17]], [[19]] ] - after insertion for tuples
  #
  #  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
  # [1,1,1,0,1,1,0,0,1,0, 0, 0, 1, 1, 1, 1] - original
  # [1,1,1,1,0,1,1,0,0,1,0, 0,1, 0, 1, 1, 1,1, 1,1] - after insertion for bool single
  rm = rm.add(
    next_to_i_seqs       = flex.size_t([0,10,14,15]),
    sites_individual     = True,
    sites_torsion_angles = True,
    sites_rigid_body     = True,
    adp_individual_iso   = True,
    adp_individual_aniso = True,
    adp_group            = True,
    group_h              = True,
    adp_tls              = True,
    s_occupancies        = True)
  out = StringIO()
  rm.show(log = out)
  assert not show_diff(out.getvalue(), """\
Refinement flags and selection counts:
  individual_sites       = %s (14 atoms)
  torsion_angles         = %s (14 atoms)
  rigid_body             = %s (14 atoms in 9 groups)
  individual_adp         = %s (iso = 14 aniso = 14)
  group_adp              = %s (14 atoms in 9 groups)
  tls                    = %s (14 atoms in 9 groups)
  occupancies            = %s (14 atoms)
  group_anomalous        = %s
""" % tuple(["%5s" % str(True)]*8))
  barr_result = flex.int([1,1,1,1,0,1,1,0,0,1,0,0,1,0,1,1,1,1,1,1]).as_bool()
  iarr_result = [ flex.size_t([0]), flex.size_t([1]), flex.size_t([2,3]),
    flex.size_t([5,6]), flex.size_t([9]), flex.size_t([14,15,16,17]),
    flex.size_t([18]), flex.size_t([19]), flex.size_t([12])]
  oarr_result = [ [flex.size_t([0]),flex.size_t([2,3])], [flex.size_t([5,6])],
                  [flex.size_t([9])], [flex.size_t([14,15,16]),flex.size_t([18])],
                  [flex.size_t([1])], [flex.size_t([12])], [flex.size_t([17])],
                  [flex.size_t([19])] ]
  assert approx_equal(rm.sites_individual      , barr_result)
  assert approx_equal(rm.sites_torsion_angles  , barr_result)
  assert approx_equal(rm.sites_rigid_body      , iarr_result)
  assert approx_equal(rm.adp_individual_iso    , barr_result)
  assert approx_equal(rm.adp_individual_aniso  , barr_result)
  assert approx_equal(rm.adp_group             , iarr_result)
  assert approx_equal(rm.group_h               , iarr_result)
  assert approx_equal(rm.adp_tls               , iarr_result)
  compare_selections(rm.s_occupancies, oarr_result)

def compare_selections(x, y):
  x = [[list(k) for k in i] for i in x][:]
  y = [[list(k) for k in i] for i in y][:]
  x.sort()
  y.sort()
  assert approx_equal(x, y)


def exercise():
  exercise_deepcopy_show_select()
  exercise_deepcopy_show_select_compare_arrays()
  exercise_inflate()
  exercise_add_1a()
  exercise_add_1b()
  exercise_add_1c()
  exercise_add_2b()
  exercise_add_2c()
  print format_cpu_times()

if(__name__ == "__main__"):
   exercise()
