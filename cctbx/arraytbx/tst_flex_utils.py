import math
from cctbx_boost.arraytbx import flex
from cctbx_boost.arraytbx import flex_utils

def exercise_misc():
  x = flex.double((0, 1, 2, 3, 4))
  flex_utils.in_place_pow(x, 2)
  assert list(x) == [math.pow(i, 2) for i in xrange(5)]
  x = flex.double((2, 4, 3, 0, 1))
  flex_utils.set_if_less_than(x, 2.5, -1)
  assert tuple(x) == (-1, 4, 3, -1, -1)

def exercise_regression_and_statistics():
  x = flex.double((0, 1, 2, 3))
  y = flex.double((1, 3, 5, 7))
  r = flex_utils.linear_regression(x, y)
  assert r.is_well_defined()
  assert abs(r.b() - 1) < 1.e-6
  assert abs(r.m() - 2) < 1.e-6
  assert abs(r.cc() - 1) < 1.e-6
  for flex_type in (flex.double, flex.float):
    x = flex_type((0, 1, 2, 3))
    s = flex_utils.statistics(x)
    assert (s.min() - 0) < 1.e-6
    assert (s.max() - 3) < 1.e-6
    assert (s.mean() - 6./4.) < 1.e-6
    assert (s.mean2() - 14./4.) < 1.e-6
    assert (s.sigma() - math.sqrt(14./4. - 36./16.)) < 1.e-6

def exercise_map_utils():
  a = flex.double((1,2,0,3,4,0))
  flex_utils.inplace_unpad(a, (1,2,2), (1,2,3))
  assert tuple(a) == (1,2,3,4)
  flex_utils.inplace_unpad(a, (1,2,2), (1,2,2))
  assert tuple(a) == (1,2,3,4)

def exercise_export():
  a = flex.double(60)
  c = flex_utils.as_CObjectZYXfloat(a, (3,4,5), (0,0,0), (3,4,5), 0)
  a = flex.float(60)
  c = flex_utils.as_CObjectZYXfloat(a, (3,4,5), (0,0,0), (3,4,5), 1)

def run(iterations):
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_misc()
    exercise_regression_and_statistics()
    exercise_map_utils()
    exercise_export()
    i += 1

if (__name__ == "__main__"):
  import sys
  from cctbx.development import debug_utils
  Flags = debug_utils.command_line_options(sys.argv[1:], (
  ))
  n = 1
  if (len(sys.argv) > 1 + Flags.n):
    n = int(Flags.regular_args[0])
  run(n)
  print "OK"
