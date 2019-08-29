from __future__ import absolute_import, division, print_function
from scitbx import fftpack
from scitbx.array_family import flex
import omptbx # initializes OpenMP environment
import libtbx.utils
from libtbx.utils import time_log
import random
import sys
from six.moves import range

def one_series(num_threads, n_iterations, quick=False):
  if (num_threads > 0):
    omptbx.env.num_threads = num_threads
  print("working omptbx.env.num_threads:", omptbx.env.num_threads)
  use_wall_clock = (omptbx.env.num_threads > 1)
  print("use_wall_clock:", use_wall_clock)
  #
  if (quick):
    dims = (2*3, 3*4, 4*5)
  else:
    dims = (2*3*5*7,3*4*5*7,3*4*5*5)
  rfft = fftpack.real_to_complex_3d(dims)
  print("rfft.m_real():", rfft.m_real())
  #
  t_map = time_log(label="map", use_wall_clock=use_wall_clock)
  t_fill = time_log(label="fill", use_wall_clock=use_wall_clock)
  t_fft = time_log(label="fft", use_wall_clock=use_wall_clock)
  print(t_map.legend)
  sys.stdout.flush()
  for i_iteration in range(n_iterations):
    t_map.start()
    map = fftpack.zeros_parallel_double(flex_grid=flex.grid(rfft.m_real()))
    print(t_map.log())
    sys.stdout.flush()
    t_fill.start()
    for i in range(0, map.size(), 97):
      map[i] = random.random()
    print(t_fill.log())
    sys.stdout.flush()
    t_fft.start()
    map = rfft.forward(map)
    print(t_fft.log())
    sys.stdout.flush()

def run(args, num_threads=0, n_iterations=3):
  show_times = libtbx.utils.show_times()
  quick = False
  rest = []
  for arg in args:
    if (arg == "--quick"):
      quick = True
    else:
      rest.append(arg)
  args = rest
  assert len(args) <= 2, "[--quick] [num_threads, [n_iterations]]"
  if (len(args) > 0):
    num_threads = int(args[0])
    if (len(args) > 1):
      n_iterations = int(args[1])
  default_num_threads = omptbx.env.num_threads
  print("default omptbx.env.num_threads:", default_num_threads)
  if (num_threads >= 0):
    if (num_threads > 0):
      print("num_threads from command line:", num_threads)
    one_series(
      num_threads=num_threads, n_iterations=n_iterations, quick=quick)
  else:
    for num_threads in range(1, default_num_threads+1):
      one_series(
        num_threads=num_threads, n_iterations=n_iterations, quick=quick)
      print()
  show_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
