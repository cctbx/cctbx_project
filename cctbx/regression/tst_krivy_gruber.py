from cctbx.uctbx import reduction_base
from cctbx.uctbx import krivy_gruber_1976
from cctbx.uctbx import gruber_1973
from cctbx.uctbx import gruber_1973_table_1
from cctbx import uctbx
from cctbx import sgtbx
from scitbx import matrix
from scitbx.python_utils.misc import get_caller_name
from libtbx.utils import time_log
from libtbx.test_utils import Exception_expected, approx_equal
import math
import random
import sys

class check_is_niggli_cell(reduction_base.gruber_parameterization):

  def itva_is_niggli_cell(self):
    eq = self.eps_eq
    gt = self.eps_gt
    a,b,c,d,e,f = (self.a,self.b,self.c,self.d,self.e,self.f)
    if (not self.meets_main_conditions()): return False
    type = self.type()
    assert type in (1,2)
    if (type == 1):
      if (eq(a, b)):
        if (gt(d, e)): return False
      if (eq(b, c)):
        if (gt(e, f)): return False
      if (eq(d, b)):
        if (gt(f, e+e)): return False
      if (eq(e, a)):
        if (gt(f, d+d)): return False
      if (eq(f, a)):
        if (gt(e, d+d)): return False
    else:
      if (eq(a, b)):
        if (gt(abs(d), abs(e))): return False
      if (eq(b, c)):
        if (gt(abs(e), abs(f))): return False
      if (eq(abs(d), b)):
        if (not eq(f, 0)): return False
      if (eq(abs(e), a)):
        if (not eq(f, 0)): return False
      if (eq(abs(f), a)):
        if (not eq(e, 0)): return False
      if (eq(abs(d)+abs(e)+abs(f), a+b)):
        if (gt(a, abs(e) + abs(f)/2)): return False
    return True

class reduction_with_tracking(krivy_gruber_1976.reduction):

  def __init__(self, unit_cell, relative_epsilon=0, iteration_limit=100):
    self.action_log = []
    self.cell_log = []
    self.type_log = []
    self.meets_primary_conditions_log = []
    self.meets_main_conditions_log = []
    self.is_niggli_cell_log = []
    try:
      krivy_gruber_1976.reduction.__init__(self,
        unit_cell=unit_cell,
        relative_epsilon=relative_epsilon,
        iteration_limit=iteration_limit)
      self.iteration_limit_exceeded = False
    except StopIteration:
      self.iteration_limit_exceeded = True

  def cb_update(self, m_elems):
    caller = get_caller_name()
    if (caller == "n3_true_action"):
      id = 3
    elif (caller == "n3_false_action"):
      id = 4
    else:
      id = int(caller[1])
    self.action_log.append(id)
    self.cell_log.append(self.as_unit_cell())
    self.type_log.append(self.type())
    self.meets_primary_conditions_log.append(self.meets_primary_conditions())
    self.meets_main_conditions_log.append(self.meets_main_conditions())
    self.is_niggli_cell_log.append(self.is_niggli_cell())
    if (self._n_iterations == self._iteration_limit):
      raise StopIteration
    self._r_inv *= matrix.sqr(m_elems)
    self._n_iterations += 1

class reduction_with_tracking_and_eq_always_false(reduction_with_tracking):

  def __init__(self, unit_cell):
    reduction_with_tracking.__init__(self, unit_cell)

  def eps_eq(self, x, y):
    return False

relative_epsilon = None
track_infinite = False
eq_always_false = False
time_krivy_gruber_1976 = time_log("krivy_gruber_1976.reduction")
time_gruber_1973 = time_log("gruber_1973.reduction")
time_krivy_gruber_1976_minimum=time_log("krivy_gruber_1976.minimum_reduction")
time_gruber_1973_minimum = time_log("gruber_1973.minimum_reduction")
time_gruber_1973_fast_minimum = time_log("gruber_1973.fast_minimum_reduction")
time_uctbx_fast_minimum = time_log("uctbx.fast_minimum_reduction")
fast_minimum_reduction_max_n_iterations = 0

def do_reduce(inp):
  assert not inp.is_degenerate()
  time_krivy_gruber_1976.start()
  red = krivy_gruber_1976.reduction(
    inp, relative_epsilon=relative_epsilon)
  time_krivy_gruber_1976.stop()
  assert red.is_niggli_cell()
  red_cell = red.as_unit_cell()
  assert inp.change_basis(red.r_inv().elems).is_similar_to(red_cell)
  if (relative_epsilon is None):
    assert check_is_niggli_cell(red_cell).itva_is_niggli_cell()
  time_gruber_1973.start()
  gruber_reduction = gruber_1973.reduction(
    inp, relative_epsilon=relative_epsilon)
  time_gruber_1973.stop()
  assert gruber_reduction.is_buerger_cell()
  buerger_cell = gruber_reduction.as_unit_cell()
  assert inp.change_basis(gruber_reduction.r_inv().elems).is_similar_to(
    buerger_cell)
  time_krivy_gruber_1976_minimum.start()
  min_reduction = krivy_gruber_1976.minimum_reduction(inp)
  time_krivy_gruber_1976_minimum.stop()
  assert min_reduction.type() in (1,2)
  min_cell = min_reduction.as_unit_cell()
  assert approx_equal(min_cell.parameters()[:3], red_cell.parameters()[:3])
  assert inp.change_basis(min_reduction.r_inv().elems).is_similar_to(min_cell)
  time_gruber_1973_minimum.start()
  min_reduction = gruber_1973.minimum_reduction(inp)
  time_gruber_1973_minimum.stop()
  assert min_reduction.type() in (1,2)
  min_cell = min_reduction.as_unit_cell()
  assert approx_equal(min_cell.parameters()[:3], red_cell.parameters()[:3])
  assert inp.change_basis(min_reduction.r_inv().elems).is_similar_to(min_cell)
  time_gruber_1973_fast_minimum.start()
  min_reduction = gruber_1973.fast_minimum_reduction(inp)
  time_gruber_1973_fast_minimum.stop()
  assert min_reduction.type() in (1,2)
  min_cell = min_reduction.as_unit_cell()
  assert approx_equal(min_cell.parameters()[:3], red_cell.parameters()[:3])
  assert inp.change_basis(min_reduction.r_inv().elems).is_similar_to(min_cell)
  time_uctbx_fast_minimum.start()
  min_reduction = uctbx.fast_minimum_reduction(inp)
  time_uctbx_fast_minimum.stop()
  assert min_reduction.type() in (1,2)
  min_cell = min_reduction.as_unit_cell()
  assert approx_equal(min_cell.parameters()[:3], red_cell.parameters()[:3])
  assert inp.change_basis(min_reduction.r_inv()).is_similar_to(min_cell)
  global fast_minimum_reduction_max_n_iterations
  fast_minimum_reduction_max_n_iterations = max(
    fast_minimum_reduction_max_n_iterations, min_reduction.n_iterations())
  if (track_infinite):
    if (eq_always_false):
      red0 = reduction_with_tracking_and_eq_always_false(inp)
    else:
      red0 = reduction_with_tracking(inp)
    if (red0.iteration_limit_exceeded):
      n = 20
      print inp
      print "red0.action_log:", red0.action_log[-n:]
      print "red0.type_log:", red0.type_log[-n:]
      print "red0.meets_primary_conditions_log:", \
             red0.meets_primary_conditions_log[-n:]
      print "red0.meets_main_conditions_log:", \
             red0.meets_main_conditions_log[-n:]
      print "red0.is_niggli_cell_log:", red0.is_niggli_cell_log[-n:]
      if (1):
        for cell in red0.cell_log[-n:]:
          print cell
        print
      sys.stdout.flush()
  return red

def reduce(inp):
  try:
    return do_reduce(inp)
  except:
    print "Problem parameters:", inp.parameters()
    raise

def ucgmx((a,b,c,d,e,f)): # unit cell given Gruber matrix
  return uctbx.unit_cell(metrical_matrix=(a,b,c,f/2.,e/2.,d/2.))

def exercise_gruber_1973_example():
  start = ucgmx((4,136,76,-155,-31,44))
  assert start.is_similar_to(uctbx.unit_cell(
    (2, 11.66, 8.718, 139+40/60., 152+45/60., 19+24/60.)))
  buerger = ucgmx((4,16,16,-16,-1,-3))
  assert buerger.is_similar_to(uctbx.unit_cell(
    (2, 4, 4, 120, 93.5833, 100.807)))
  niggli = ucgmx((4,16,16,16,3,4))
  assert niggli.is_similar_to(uctbx.unit_cell(
    (2, 4, 4, 60, 79.1931, 75.5225)))
  red = reduction_base.gruber_parameterization(start)
  assert not red.is_buerger_cell()
  assert approx_equal(red.as_gruber_matrix(), (4,136,76,-155,-31,44))
  assert approx_equal(red.as_niggli_matrix(), (4,136,76,-155/2.,-31/2.,44/2.))
  assert approx_equal(red.as_sym_mat3(), (4,136,76,44/2.,-31/2.,-155/2.))
  assert red.as_unit_cell().is_similar_to(start)
  red = reduction_base.gruber_parameterization(buerger)
  assert red.is_buerger_cell()
  assert not red.is_niggli_cell()
  red = reduction_base.gruber_parameterization(niggli)
  assert red.is_niggli_cell()
  red = reduce(start)
  assert red.as_unit_cell().is_similar_to(niggli)
  assert red.r_inv().elems == (-1, 5, 9, 0, -1, -1, 0, 0, 1)
  assert red.n_iterations() == 29
  red = reduce(buerger)
  assert red.as_unit_cell().is_similar_to(niggli)
  assert red.r_inv().elems == (-1, 0, 0, 0, 1, 1, 0, 1, 0)
  assert red.n_iterations() == 4
  red = reduce(niggli)
  assert red.as_unit_cell().is_similar_to(niggli)
  assert red.r_inv().elems == (1, 0, 0, 0, 1, 0, 0, 0, 1)
  assert red.n_iterations() == 1
  try:
    red = krivy_gruber_1976.reduction(buerger, iteration_limit=1)
  except krivy_gruber_1976.iteration_limit_exceeded:
    pass
  else:
    raise Exception_expected
  assert not start.is_buerger_cell()
  assert not start.is_niggli_cell()
  assert buerger.is_buerger_cell()
  assert not buerger.is_niggli_cell()
  assert niggli.is_buerger_cell()
  assert niggli.is_niggli_cell()
  red = start.niggli_reduction()
  assert red.n_iterations() == 29
  assert start.niggli_cell().is_similar_to(niggli)

def exercise_krivy_gruber_1976_example():
  start = ucgmx((9,27,4,-5,-4,-22))
  assert start.is_similar_to(uctbx.unit_cell(
    (3, 5.196, 2, 103+55/60., 109+28/60., 134+53/60.)))
  for gmx in ((4,9,9,-8,-1,-4),
              (4,9,9,9,1,4)):
    red = reduction_base.gruber_parameterization(ucgmx(gmx))
    assert red.is_buerger_cell()
    assert not red.is_niggli_cell()
  niggli = ucgmx((4,9,9,9,3,4))
  assert niggli.is_similar_to(uctbx.unit_cell(
    (2, 3, 3, 60, 75+31/60., 70+32/60.)))
  red = reduction_base.gruber_parameterization(niggli)
  assert red.is_niggli_cell()
  red = reduce(start)
  assert red.as_unit_cell().is_similar_to(niggli)
  assert red.r_inv().elems == (0, 1, 2, 0, 0, 1, 1, 1, 2)
  assert red.n_iterations() == 11

def exercise_bravais_plus():
  for pg in ("1", "2", "2 2", "4", "3*", "6", "2 2 3"):
    for z in "PABCIRF":
      sgi = sgtbx.space_group_info("Hall: %s %s" % (z, pg))
      r_inv = sgi.group().z2p_op().c_inv().r()
      reduce(sgi.any_compatible_unit_cell(volume=100).change_basis(
        r_inv.num(),r_inv.den()))

def cos_deg(x):
  return math.cos(x*math.pi/180)

def exercise_grid(quick=False, verbose=0):
  if (quick):
    sample_lengths = (10,)
    sample_angles = (60,)
  else:
    sample_lengths = (10,20,30)
    sample_angles = (10,30,45,60,90,120,150,170)
  n_trials = 0
  for a in sample_lengths:
    for b in sample_lengths:
      for c in sample_lengths:
        for alpha in sample_angles:
          for beta in sample_angles:
            for gamma in sample_angles:
              a_b = a*b*cos_deg(gamma)
              a_c = a*c*cos_deg(beta)
              b_c = b*c*cos_deg(alpha)
              g = matrix.sqr((a*a,a_b,a_c,
                              a_b,b*b,b_c,
                              a_c,b_c,c*c))
              det_g = g.determinant()
              try: unit_cell = uctbx.unit_cell((a,b,c,alpha,beta,gamma))
              except:
                assert det_g <= 1.e-5
                continue
              assert abs(det_g-unit_cell.volume()**2) < 1.e-5
              if (unit_cell.volume() < a*b*c/1000.): continue
              n_trials += 1
              reduce(unit_cell)
  if (0 or verbose):
    print "exercise_grid n_trials:", n_trials

class random_unimodular_integer_matrix_generator(object):

  def __init__(self, reset_threshold=10):
    self.reset_threshold = reset_threshold
    self._m1 = matrix.sqr((0,0,1,1,0,0,0,1,0))
    self._m2 = matrix.sqr((1,-1,0,1,0,0,0,0,1))
    self._mi = self._m1 * self._m2

  def has_elements_which_are_to_large(self):
    e = self._mi.elems
    return max(abs(min(e)), abs(max(e))) >= self.reset_threshold

  def next(self):
    while 1:
      if (random.randrange(0,2)):
        self._mi = self._m2 * self._mi
      else:
        self._mi = self._m1 * self._mi
      if (not self.has_elements_which_are_to_large()):
        break
      self._mi = random.choice((self._m1, self._m2))
    return self._mi

class random_abcpq(object):

  def __init__(self, ck_type):
    rr = random.randrange
    self.a = rr(100,201)
    self.b = self.a
    self.p = rr(10,41)
    self.q = self.a
    if (ck_type[0] == "q"):
      if (ck_type[1] != "=" and rr(0,2)):
        self.q = rr(50,91)
      ck_type = ck_type[2:]
    if   (ck_type == "a=b<=c"):
      self.c = self.b
      if (rr(0,2)): self.c += rr(10,101)
    elif (ck_type == "a<b=c"):
      self.b += rr(10,101)
      self.c = self.b
    elif (ck_type == "a<=b<c"):
      if (rr(0,2)): self.b += rr(10,101)
      self.c = self.b + rr(10,101)
    elif (ck_type == "a<b<c"):
      self.b += rr(10,101)
      self.c = self.b + rr(10,101)
    elif (ck_type == "a=b<c"):
      self.c = self.b + rr(10,101)
    elif (ck_type == "a<b<=c"):
      self.b += rr(10,101)
      self.c = self.b
      if (rr(0,2)): self.c += rr(10,101)
    else:
      raise RuntimeError, "Unknown ck_type."

  def eval_defks(self, defks):
    a,b,c,p,q = tuple([float(x) for x in (self.a,self.b,self.c,self.p,self.q)])
    m = b/a
    n = (b-a)/a
    d,e,f = eval(defks)
    return a,b,c,d,e,f

def random_gruber_matrix(type_conditions):
  return random_abcpq(random.choice(
    type_conditions.ck_types)).eval_defks(type_conditions.defks)

def exercise_gruber_types(n_trials_per_type, dump=0, verbose=0):
  mk2_sets = gruber_1973_table_1.get_mk2_sets()
  type_conditions = gruber_1973_table_1.get_type_conditions()
  random_unimodular = random_unimodular_integer_matrix_generator()
  have_errors = False
  for k in xrange(1,29):
    set = mk2_sets[k]
    tc = type_conditions[k]
    if (0 or verbose):
      print " ", tc.ck_types, tc.defks
    n_errors = 0
    for i_trial in xrange(n_trials_per_type):
      gruber_matrix = random_gruber_matrix(tc)
      type_cell = ucgmx(gruber_matrix)
      if (0 or verbose):
        print " ", gruber_matrix
        print " ", type_cell
      red = reduction_base.gruber_parameterization(type_cell)
      assert red.is_niggli_cell()
      n_different_cells = 0
      for m in set:
        other_cell = type_cell.change_basis(m.inverse().transpose().elems, 1)
        if (0 or verbose):
          print " ", m.elems, m.determinant()
          print " ", other_cell
        red = reduction_base.gruber_parameterization(other_cell)
        if (not red.is_buerger_cell()):
          print "  Error: Transformed cell is not a Buerger cell."
          print "  gruber_matrix:", gruber_matrix
          n_errors += 1
        else:
          n_different_cells += 1
          if (red.is_niggli_cell()):
            if (not other_cell.is_similar_to(type_cell)):
              print "  Error: Transformed cell is a Niggli cell."
              print "  gruber_matrix:", gruber_matrix
              n_errors += 1
          else:
            krivy_cell = reduce(type_cell).as_unit_cell()
            assert krivy_cell.is_similar_to(type_cell)
            krivy_cell = reduce(other_cell).as_unit_cell()
            assert krivy_cell.is_similar_to(type_cell)
            r_inv = random_unimodular.next().elems
            random_cell = type_cell.change_basis(r_inv, 1)
            if (0 or verbose):
              print "  Random:", random_cell, r_inv
            red = reduce(random_cell)
            krivy_cell = red.as_unit_cell()
            if (dump):
              print "type_cell:", type_cell
              print "random_cell:", random_cell
              print "krivy_cell:", krivy_cell
              print
            if (not krivy_cell.is_similar_to(type_cell)):
              print "  Error: Random cell recovery."
              print "  gruber_matrix:", gruber_matrix
              print "  r_inv:", r_inv
              print "  red.as_gruber_matrix():", red.as_gruber_matrix()
              n_errors += 1
      if (n_different_cells == 0):
        print "  Error: Transformation does not yield different cells."
        n_errors += 1
        raise RuntimeError
    if ((0 or verbose) and n_errors != 0):
      print "Errors for type %d:" % k, n_errors
    if (n_errors != 0):
      have_errors = True
  assert not have_errors

def exercise_extreme():
  uc = uctbx.unit_cell((
    69.059014477286041, 48.674386086971339, 0.0048194797114296736,
    89.995145576185806, 89.999840576946085, 99.484656090034875))
  red = uc.niggli_reduction()
  assert red.as_unit_cell().is_similar_to(uctbx.unit_cell((
    0.00481948, 48.6744, 69.059,
    80.5153, 89.9962, 89.9951)))
  uc = uctbx.unit_cell((
    80.816186392181365, 81.021289502648813, 140.6784408482614,
    29.932540128999769, 89.92047105556459, 119.85301114570319))
  uc.niggli_reduction(iteration_limit=10000)

def exercise_real_world_examples():
  # SSZ-59, cell by Michael Treacy, infinite loop in GSAS rducll (Linux)
  uc = uctbx.unit_cell((
    12.7366, 29.2300, 5.0242,
    94.6570, 100.8630, 99.7561))
  nc = uc.niggli_cell()
  assert nc.is_similar_to(uctbx.unit_cell(
    (5.0242, 12.7366, 29.23, 99.7561, 94.657, 100.863)))
  # SSZ-59, Burton et al., Table 4
  uc = uctbx.unit_cell((
    12.7806, 12.7366, 29.457,
    103.42, 103.57, 22.71))
  red = uc.niggli_reduction()
  assert red.as_unit_cell().is_similar_to(nc)
  assert red.r_inv().elems == (-1, 0, 1, 1, -1, 0, 0, 0, 1)

def exercise_problem_parameters():
  problem_parameters = (
    (13.892443989449804, 13.892443989449804, 14.7648230602334,
     61.936000954634402, 61.936000954634402, 88.515487291567879),
    (10.0, 10.0, 20.0,
     90.0, 45.0, 120.0),
    (10.0, 20.0, 30.0,
     120.0, 60.0, 120.0),
    (10.816653826391969, 13.820274961085254, 13.820274961085254,
     60.0, 66.962544368849834, 66.962544368849834),
    (10.148891565092219, 13.379088160259652, 13.379088160259652,
     107.33461190548745, 107.94415159713115, 109.72759194290299),
    (19.798989873223331, 231.21851136965654, 14.352700094407323,
     133.37207519042573, 92.016673840743408, 134.55815348093702),
    (10.392304845413264, 13.19090595827292, 13.19090595827292,
     112.64730819498385, 104.36056979415913, 106.96527532101391),
    (16.046806535881213, 13.341664064126334, 197.64614845728718,
     153.28759931491018, 114.05435960569044, 92.543256980798247),
    (10.488088481701515, 13.820274961085254, 13.820274961085254,
     109.9226012907464, 104.00699650533103, 110.31922490992999),
    (10.04987562112089, 13.19090595827292, 13.19090595827292,
     118.05419482122835, 97.404049814230376, 106.92070123011929),
    (10.04987562112089, 13.45362404707371, 13.45362404707371,
     109.02416163919622, 105.88181549565937, 109.44017310001107),
    (11.357816691600547, 13.638181696985853, 13.638181696985855,
     115.81608733396159, 104.29612977641231, 104.29612977641233),
    (11.832159566199232, 13.784048752090222, 13.784048752090222,
     110.67521616123457, 104.95317005195066, 110.01926787579129))
  for parameters in problem_parameters:
    reduce(uctbx.unit_cell(parameters))

def exercise_one_pass(show_times=False):
  quick = "--Quick" in sys.argv[1:]
  verbose = "--Verbose" in sys.argv[1:]
  global relative_epsilon
  global track_infinite
  global eq_always_false
  if ("--zero_epsilon" in sys.argv[1:]):
    relative_epsilon = 0
  if ("--track_infinite" in sys.argv[1:]):
    track_infinite = True
  if ("--eq_always_false" in sys.argv[1:]):
    track_infinite = True
    eq_always_false = True
  exercise_problem_parameters()
  exercise_extreme()
  exercise_gruber_1973_example()
  exercise_krivy_gruber_1976_example()
  exercise_bravais_plus()
  exercise_grid(quick=quick, verbose=verbose)
  if (quick): n_trials_per_type=10
  else:       n_trials_per_type=100
  if ("--dump" in sys.argv[1:]):
    random.seed(0)
    exercise_gruber_types(n_trials_per_type, dump=True, verbose=verbose)
    return
  exercise_gruber_types(n_trials_per_type, verbose=verbose)
  exercise_real_world_examples()
  if (0 or verbose or show_times):
    print time_krivy_gruber_1976.report()
    print time_gruber_1973.report()
    print time_krivy_gruber_1976_minimum.report()
    print time_gruber_1973_minimum.report()
    print time_gruber_1973_fast_minimum.report()
    print time_uctbx_fast_minimum.report()
    if (time_uctbx_fast_minimum.accumulation != 0):
      print "fast_minimum Python/C++: %.3g" % (
          time_gruber_1973_fast_minimum.accumulation
        / time_uctbx_fast_minimum.accumulation)
    print "fast_minimum_reduction_max_n_iterations:", \
           fast_minimum_reduction_max_n_iterations
  print "OK"

def exercise():
  forever = "--Forever" in sys.argv[1:]
  while 1:
    exercise_one_pass(show_times=forever)
    sys.stdout.flush()
    if (not forever):
      break

if (__name__ == "__main__"):
  exercise()
