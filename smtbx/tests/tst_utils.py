from cctbx import sgtbx
from libtbx.test_utils import approx_equal
from smtbx import development, utils

def exercise_connectivity_table():
  xs = development.sucrose()
  connectivity = utils.connectivity_table(xs)
  pair_counts = [
    2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 4, 1, 4, 1,
    4, 1, 4, 1, 1, 4, 1, 4, 1, 4, 4, 1, 1, 4, 1, 4, 1, 4, 1, 4, 1, 1]
  assert approx_equal(connectivity.pair_asu_table.pair_counts(), pair_counts)
  connectivity.add_bond(0, 1)
  assert approx_equal(
    connectivity.pair_asu_table.pair_counts(),
    [3, 3, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 4, 1, 4, 1,
     4, 1, 4, 1, 1, 4, 1, 4, 1, 4, 4, 1, 1, 4, 1, 4, 1, 4, 1, 4, 1, 1])
  connectivity.add_bond(5, 5, rt_mx_ji=sgtbx.rt_mx("x+1,y,z"))
  assert connectivity.pair_asu_table.pair_counts()[5] == 4
  connectivity.remove_bond(0, 1)
  connectivity.remove_bond(5,5, rt_mx_ji=sgtbx.rt_mx("x+1,y,z"))
  assert approx_equal(connectivity.pair_asu_table.pair_counts(), pair_counts)

def run():
  exercise_connectivity_table()
  print "OK"

if __name__ == '__main__':
  run()
