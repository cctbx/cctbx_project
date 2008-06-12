from scitbx import fftpack
from scitbx.array_family import flex
import omptbx # initializes OpenMP environment
import libtbx.utils
from libtbx.utils import time_log
import random
import sys

def run(args, num_threads=0, n_iterations=3):
  show_times = libtbx.utils.show_times()
  assert len(args) <= 2, "[num_threads, [n_iterations]]"
  if (len(args) > 0):
    num_threads = int(args[0])
    if (len(args) > 1):
      n_iterations = int(args[1])
  #
  print "default omptbx.env.num_threads:", omptbx.env.num_threads
  if (num_threads > 0):
    print "num_threads from command line:", num_threads
    omptbx.env.num_threads = num_threads
    print "working omptbx.env.num_threads:", omptbx.env.num_threads
  use_wall_clock = (omptbx.env.num_threads > 1)
  print "use_wall_clock:", use_wall_clock
  #
  rfft = fftpack.real_to_complex_3d((2*3*5*7,3*4*5*7,3*4*5*5))
  print "rfft.m_real():", rfft.m_real()
  #
  t_map = time_log(label="map", use_wall_clock=use_wall_clock)
  t_fill = time_log(label="fill", use_wall_clock=use_wall_clock)
  t_fft = time_log(label="fft", use_wall_clock=use_wall_clock)
  print t_map.legend
  for i_iteration in xrange(n_iterations):
    t_map.start()
    map = flex.double(flex.grid(rfft.m_real()))
    print t_map.log()
    t_fill.start()
    for i in xrange(0, map.size(), 97):
      map[i] = random.random()
    print t_fill.log()
    t_fft.start()
    map = rfft.forward(map)
    print t_fft.log()
  show_times()

if (__name__ == "__main__"):
  run(sys.argv[1:])
