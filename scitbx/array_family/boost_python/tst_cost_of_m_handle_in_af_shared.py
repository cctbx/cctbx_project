from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
import time
from six.moves import range

def exercise(data_size, n_repeats):
  print("data_size=%i, n_repeats=%i" % (data_size, n_repeats))
  data = flex.random_double(data_size)
  cost = flex.cost_of_m_handle_in_af_shared(data)
  t0 = time.time()
  cost(n_repeats=0, test_id=0)
  t1 = time.time()
  t_overhead = t1 - t0
  for test_id in range(3):
    t0 = time.time()
    label = cost(n_repeats=n_repeats, test_id=test_id)
    t1 = time.time()
    print('%s: %.2f s' % (label, t1 - t0 - t_overhead))
  print()

def run():
  print("Loop over af::shared<double>:")
  for data_size in [1000, 10000, 100000, 1000000, 10000000]:
    n_repeats = 30000000 // data_size
    exercise(data_size, n_repeats)
  print('OK')

if __name__ == '__main__':
  run()
