import math
from cctbx_boost.arraytbx import shared
from cctbx_boost.arraytbx import shared_map

def exercise_utils():
  a = shared.double((1,2,0,3,4,0))
  shared_map.inplace_unpad(a, (1,2,2), (1,2,3))
  assert tuple(a) == (1,2,3,4)
  shared_map.inplace_unpad(a, (1,2,2), (1,2,2))
  assert tuple(a) == (1,2,3,4)

def exercise_export():
  a = shared.double(60)
  c = shared_map.as_CObjectZYXfloat(a, (3,4,5), (0,0,0), (3,4,5), 0)
  a = shared.float(60)
  c = shared_map.as_CObjectZYXfloat(a, (3,4,5), (0,0,0), (3,4,5), 1)

def run(iterations):
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_utils()
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
