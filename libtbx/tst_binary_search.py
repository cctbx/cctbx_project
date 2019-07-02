from __future__ import absolute_import, division, print_function
from libtbx import binary_search
import random
import sys
from six.moves import range

class random_callback(object):

  def __init__(O, rng, low, high):
    O.low = low
    O.high = high
    O.critical_point = min(
      low + max(1, int((high-low+1) * rng.random())),
      high)
    n_bad = int((high-low-1) * 0.2 + 0.5)
    O.bad_points = set()
    while (len(O.bad_points) != n_bad):
      n_block = min(rng.randrange(1,21), n_bad - len(O.bad_points))
      begin = low + 1 + rng.randrange(high-low-1)
      assert begin < high
      for i in range(begin, min(begin+n_block,high)):
        if (i != O.critical_point):
          O.bad_points.add(i)
          if (len(O.bad_points) == n_bad):
            break

  def __call__(O, point):
    assert point > O.low
    assert point < O.high
    if (point in O.bad_points):
      return None
    return (point < O.critical_point)

def run(args):
  assert len(args) == 0
  #
  info = binary_search.true_false_bad_biased_up(
    low=17, high=18, callback=None)
  assert info.number_of_iterations == 0
  #
  points_tested = []
  def callback(point):
    points_tested.append(point)
    return None
  info = binary_search.true_false_bad_biased_up(
    low=0, high=9, callback=callback)
  assert points_tested == [4,6,7,8,5,2,3,1] # see example in docstring
  del points_tested[:]
  info = binary_search.true_false_bad_biased_up(
    low=30, high=39, callback=callback)
  assert points_tested == [34,36,37,38,35,32,33,31]
  #
  for i_trial in range(32):
    callback = random_callback(
      rng=random.Random(i_trial),
      low=i_trial*13,
      high=i_trial*14+567)
    info = binary_search.true_false_bad_biased_up(
      low=callback.low,
      high=callback.high,
      callback=callback)
    assert info.low_point_false == callback.critical_point
    assert info.high_point_true < info.low_point_false
    assert info.high_point_true not in callback.bad_points
    for i in range(info.high_point_true+1, info.low_point_false):
      assert i in callback.bad_points
    assert info.bad_points.issubset(callback.bad_points)
    assert info.number_of_iterations >= info.number_of_callbacks
    assert info.bad_gap_width >= 0
  #
  print("OK")

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
