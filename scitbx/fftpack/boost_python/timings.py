from scitbx import fftpack
from scitbx.array_family import flex
from scitbx.python_utils.misc import time_log
import random
import sys

def run(n_iterations=10):
  assert len(sys.argv) in (1,2)
  if (len(sys.argv) > 1):
    n_iterations = int(sys.argv[1])
  rfft = fftpack.real_to_complex_3d((2*3*5*7,4*5*7,3*5*5))
  print "rfft.m_real():", rfft.m_real()
  t_map = time_log("map")
  t_fill = time_log("fill")
  t_fft = time_log("fft")
  print t_map.legend()
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

if (__name__ == "__main__"):
  run()
