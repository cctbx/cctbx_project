from __future__ import division

class true_false_bad_biased_up(object):
  """
Example application: binary search in svn history, with handling
of "bad" revisions that do not build.

Example callback:

  def callback(point):
    bad_points = set([3,8,10])
    critical_point = 7
    if (point in bad_points): return None
    return (point < critical_point)

The initial low point is assumed to be True, the initial high point
is assumed to be False; these assumptions are not verified, i.e.
callback is never called with point=low or point=high.

Illustration of binary search steps in region with bad points:
  T . . . . . . . F
  T . . . B . . . F
  T . . . B . B . F
  T . . . B . B B F
  T . . . B . B b F
  T . . . B . b b F
  T . . . B B b b F
  T . . . B b b b F
  T . . . b b b b F
  T . B . b b b b F
  T . B B b b b b F
  T . B b b b b b F
  T . b b b b b b F
  T B b b b b b b F
  T b b b b b b b F
Each row shows the state at a step of the search.
Each column is for a search point (e.g. svn revision).
T = point known to be True
F = point known to be False
. = point not tried yet
B = point known to be bad (e.g. the True/False decision cannot be made)
b = point known to be bad and contiguously preceding F or b

In the algorithm below, T and B are stored in the i_true_or_bad list.
This reflects that bad points are initially treated as if they were
True. If the binary search converges with a bad point preceding a False
one, the bad point is reinterpreted as False. The final result is the
last point certain to be True and the first point certain to be False.
All revisions in between, if any, are known to be bad.
"""

  __slots__ = [
    "low", "high", "bad_points", "high_point_true", "low_point_false",
    "number_of_iterations", "number_of_callbacks", "bad_gap_width"]

  def __init__(O, low, high, callback, known_bad_points=()):
    assert low < high
    O.low = low
    O.high = high
    O.bad_points = set(known_bad_points)
    i_true_or_bad = [low]
    O.low_point_false = i_false_or_bad = high
    O.number_of_iterations = 0
    O.number_of_callbacks = 0
    while (i_true_or_bad[0] + 1 != i_false_or_bad):
      if (i_true_or_bad[-1] + 1 == i_false_or_bad):
        i_false_or_bad = i_true_or_bad.pop()
      else:
        i_trial = i_true_or_bad[-1] + (i_false_or_bad - i_true_or_bad[-1]) // 2
        O.number_of_iterations += 1
        if (i_trial in O.bad_points):
          cb_result = None
        else:
          cb_result = callback(i_trial)
          O.number_of_callbacks += 1
          if (cb_result is None):
            O.bad_points.add(i_trial)
        if (cb_result is None):
          i_true_or_bad.append(i_trial)
        elif (cb_result):
          i_true_or_bad = [i_trial]
        else:
          del i_true_or_bad[1:]
          O.low_point_false = i_false_or_bad = i_trial
    assert len(i_true_or_bad) == 1
    O.high_point_true = i_true_or_bad[0]
    O.bad_gap_width = O.low_point_false - O.high_point_true - 1
