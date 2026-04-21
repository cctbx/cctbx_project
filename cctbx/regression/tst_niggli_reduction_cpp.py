"""Comprehensive comparison of C++ niggli_reduction vs Python krivy_gruber_1976.

Verifies that every public method on the C++ result (via
cpp_niggli_reduction_result) matches the Python reference implementation
for a wide range of unit cells.
"""
from __future__ import absolute_import, division, print_function
from cctbx.uctbx import krivy_gruber_1976
from cctbx.uctbx.reduction_base import iteration_limit_exceeded
from cctbx import uctbx
from libtbx.test_utils import Exception_expected

def ucgmx(gruber_matrix):
  (a,b,c,d,e,f) = gruber_matrix
  return uctbx.unit_cell(metrical_matrix=(a,b,c,f/2.,e/2.,d/2.))

def compare_reductions(unit_cell, relative_epsilon=None):
  """Compare all public methods between Python and C++ niggli reduction."""
  py_red = krivy_gruber_1976.reduction(
    unit_cell, relative_epsilon=relative_epsilon)
  cpp_red = unit_cell.niggli_reduction(relative_epsilon=relative_epsilon)

  # Reduced cell (exact: both construct unit_cell from identical G6 values
  # via the same C++ unit_cell(metrical_matrix=...) constructor)
  assert cpp_red.as_unit_cell().parameters() == \
    py_red.as_unit_cell().parameters(), \
    "as_unit_cell() mismatch: C++ %s vs Python %s" % (
      cpp_red.as_unit_cell().parameters(), py_red.as_unit_cell().parameters())

  # Change-of-basis matrix (exact integer comparison)
  assert cpp_red.r_inv().elems == py_red.r_inv().elems, \
    "r_inv() mismatch: C++ %s vs Python %s" % (
      cpp_red.r_inv().elems, py_red.r_inv().elems)

  # Iteration count
  assert cpp_red.n_iterations() == py_red.n_iterations(), \
    "n_iterations() mismatch: C++ %d vs Python %d" % (
      cpp_red.n_iterations(), py_red.n_iterations())

  # Analysis methods from gruber_parameterization
  assert cpp_red.is_niggli_cell() == py_red.is_niggli_cell(), \
    "is_niggli_cell() mismatch: C++ %s vs Python %s" % (
      cpp_red.is_niggli_cell(), py_red.is_niggli_cell())

  assert cpp_red.is_buerger_cell() == py_red.is_buerger_cell(), \
    "is_buerger_cell() mismatch: C++ %s vs Python %s" % (
      cpp_red.is_buerger_cell(), py_red.is_buerger_cell())

  assert cpp_red.type() == py_red.type(), \
    "type() mismatch: C++ %s vs Python %s" % (
      cpp_red.type(), py_red.type())

  assert cpp_red.meets_main_conditions() == py_red.meets_main_conditions(), \
    "meets_main_conditions() mismatch: C++ %s vs Python %s" % (
      cpp_red.meets_main_conditions(), py_red.meets_main_conditions())

  # change_of_basis_op
  cpp_cbop = cpp_red.change_of_basis_op()
  py_cbop = py_red.change_of_basis_op()
  assert cpp_cbop.c().r().as_double() == py_cbop.c().r().as_double(), \
    "change_of_basis_op() rotation mismatch: C++ %s vs Python %s" % (
      cpp_cbop.c().r().as_double(), py_cbop.c().r().as_double())
  assert cpp_cbop.c().t().as_double() == py_cbop.c().t().as_double(), \
    "change_of_basis_op() translation mismatch: C++ %s vs Python %s" % (
      cpp_cbop.c().t().as_double(), py_cbop.c().t().as_double())

  # Verify the change_of_basis_op actually works (tolerance-based because
  # change_basis reconstructs the cell through a different arithmetic path)
  assert unit_cell.change_basis(cpp_cbop).is_similar_to(
    cpp_red.as_unit_cell()), \
    "change_of_basis_op does not transform to reduced cell"

def exercise_textbook_examples():
  # Gruber 1973 example
  cells = [
    ucgmx((4,136,76,-155,-31,44)),   # start
    ucgmx((4,16,16,-16,-1,-3)),      # buerger
    ucgmx((4,16,16,16,3,4)),         # niggli (already reduced)
  ]
  for uc in cells:
    compare_reductions(uc)

  # Krivy-Gruber 1976 example
  cells = [
    ucgmx((9,27,4,-5,-4,-22)),       # start
    ucgmx((4,9,9,-8,-1,-4)),         # intermediate buerger
    ucgmx((4,9,9,9,1,4)),            # intermediate buerger
    ucgmx((4,9,9,9,3,4)),            # niggli (already reduced)
  ]
  for uc in cells:
    compare_reductions(uc)

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
    compare_reductions(uctbx.unit_cell(parameters))

def compare_reductions_with_limit(unit_cell, iteration_limit,
                                  relative_epsilon=None):
  """Like compare_reductions but with an explicit iteration_limit."""
  py_red = krivy_gruber_1976.reduction(
    unit_cell, relative_epsilon=relative_epsilon,
    iteration_limit=iteration_limit)
  cpp_red = unit_cell.niggli_reduction(
    relative_epsilon=relative_epsilon, iteration_limit=iteration_limit)

  assert cpp_red.as_unit_cell().parameters() == \
    py_red.as_unit_cell().parameters(), \
    "as_unit_cell() mismatch: C++ %s vs Python %s" % (
      cpp_red.as_unit_cell().parameters(), py_red.as_unit_cell().parameters())
  assert cpp_red.r_inv().elems == py_red.r_inv().elems, \
    "r_inv() mismatch: C++ %s vs Python %s" % (
      cpp_red.r_inv().elems, py_red.r_inv().elems)
  assert cpp_red.n_iterations() == py_red.n_iterations(), \
    "n_iterations() mismatch: C++ %d vs Python %d" % (
      cpp_red.n_iterations(), py_red.n_iterations())
  assert cpp_red.is_niggli_cell() == py_red.is_niggli_cell(), \
    "is_niggli_cell() mismatch: C++ %s vs Python %s" % (
      cpp_red.is_niggli_cell(), py_red.is_niggli_cell())
  assert cpp_red.type() == py_red.type(), \
    "type() mismatch: C++ %s vs Python %s" % (
      cpp_red.type(), py_red.type())

def exercise_extreme():
  uc = uctbx.unit_cell((
    69.059014477286041, 48.674386086971339, 0.0048194797114296736,
    89.995145576185806, 89.999840576946085, 99.484656090034875))
  compare_reductions(uc)
  # This cell requires a high iteration limit (>1000)
  uc = uctbx.unit_cell((
    80.816186392181365, 81.021289502648813, 140.6784408482614,
    29.932540128999769, 89.92047105556459, 119.85301114570319))
  compare_reductions_with_limit(uc, iteration_limit=10000)

def exercise_real_world_examples():
  uc = uctbx.unit_cell((
    12.7366, 29.2300, 5.0242,
    94.6570, 100.8630, 99.7561))
  compare_reductions(uc)
  uc = uctbx.unit_cell((
    12.7806, 12.7366, 29.457,
    103.42, 103.57, 22.71))
  compare_reductions(uc)

def exercise_grid():
  # Same angle set as tst_krivy_gruber.py (non-quick mode).
  n_tested = 0
  for a in (10, 20, 30):
    for b in (10, 20, 30):
      for c in (10, 20, 30):
        for alpha in (10, 30, 45, 60, 90, 120, 150, 170):
          for beta in (10, 30, 45, 60, 90, 120, 150, 170):
            for gamma in (10, 30, 45, 60, 90, 120, 150, 170):
              try:
                uc = uctbx.unit_cell((a,b,c,alpha,beta,gamma))
              except Exception:
                continue
              if uc.volume() < a*b*c/1000.:
                continue
              compare_reductions(uc)
              n_tested += 1
  assert n_tested > 1000, "Too few grid cells tested: %d" % n_tested

def exercise_bravais_plus():
  # Same space-group/cell combinations as tst_krivy_gruber.exercise_bravais_plus().
  from cctbx import sgtbx
  for pg in ("1", "2", "2 2", "4", "3*", "6", "2 2 3"):
    for z in "PABCIRF":
      sgi = sgtbx.space_group_info("Hall: %s %s" % (z, pg))
      r_inv = sgi.group().z2p_op().c_inv().r()
      uc = sgi.any_compatible_unit_cell(volume=100).change_basis(
        r_inv.num(), r_inv.den())
      compare_reductions(uc)

def exercise_gruber_types():
  # Same 28 Gruber types as tst_krivy_gruber.exercise_gruber_types(),
  # 10 random trials each.  Seed is fixed for reproducibility.
  from cctbx.uctbx import gruber_1973_table_1
  from cctbx.regression.tst_krivy_gruber import random_gruber_matrix
  import random
  random.seed(0)
  type_conditions = gruber_1973_table_1.get_type_conditions()
  for k in range(1, 29):
    tc = type_conditions[k]
    for _ in range(100):
      compare_reductions(ucgmx(random_gruber_matrix(tc)))

def exercise_iteration_limit():
  uc = ucgmx((4,16,16,-16,-1,-3))  # buerger cell, needs several iterations
  # Python raises iteration_limit_exceeded
  try:
    krivy_gruber_1976.reduction(uc, iteration_limit=1)
  except iteration_limit_exceeded:
    pass
  else:
    raise Exception_expected
  # C++ should also raise iteration_limit_exceeded
  try:
    uc.niggli_reduction(iteration_limit=1)
  except iteration_limit_exceeded:
    pass
  else:
    raise Exception_expected

def exercise_a3_a4_zero_component():
  """Targeted test for the zero-component edge case in do_a3_a4.

  A cell with alpha=beta=90 yields xi=eta=0 and zeta > 0 in G6.
  In sign_counts() this gives pos=1, nonneg=3 (pos != nonneg, pos%2 == 1),
  triggering the extra negate_column branch in the all-negative/mixed case.
  """
  uc = uctbx.unit_cell((10, 20, 30, 90, 90, 120))
  compare_reductions(uc)
  cpp_red = uc.niggli_reduction()
  assert cpp_red.is_niggli_cell(), \
    "Reduced cell is not Niggli: %s" % str(cpp_red.as_gruber_matrix())

def exercise_zero_epsilon():
  cells = [
    ucgmx((4,136,76,-155,-31,44)),
    ucgmx((4,16,16,16,3,4)),
    ucgmx((9,27,4,-5,-4,-22)),
    uctbx.unit_cell((10, 20, 30, 90, 90, 90)),
    uctbx.unit_cell((12.7366, 29.2300, 5.0242, 94.6570, 100.8630, 99.7561)),
    uctbx.unit_cell((13.892443989449804, 13.892443989449804, 14.7648230602334,
      61.936000954634402, 61.936000954634402, 88.515487291567879)),
  ]
  for uc in cells:
    compare_reductions(uc, relative_epsilon=0)

def exercise():
  exercise_textbook_examples()
  exercise_problem_parameters()
  exercise_extreme()
  exercise_real_world_examples()
  exercise_grid()
  exercise_bravais_plus()
  exercise_gruber_types()
  exercise_iteration_limit()
  exercise_a3_a4_zero_component()
  exercise_zero_epsilon()
  print("OK")

if (__name__ == "__main__"):
  exercise()
