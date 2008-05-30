from scitbx.array_family import flex
import time

def exercise(data_size, n_repeats):
  print "data_size=%i, n_repeats=%i" % (data_size, n_repeats)
  data = flex.random_double(data_size)
  cost = flex.cost_of_m_handle_in_af_shared(data)
  t0 = time.time()
  cost(n_repeats=0, use_af_shared_indexing=True)
  t1 = time.time()
  t_overhead = t1 - t0
  t0 = time.time()
  cost(n_repeats=n_repeats, use_af_shared_indexing=True)
  t1 = time.time()
  print 'a[i] = ... --> ',
  print '%.2f s' % (t1 - t0 - t_overhead)

  t0 = time.time()
  cost(n_repeats=n_repeats, use_af_shared_indexing=False)
  t1 = time.time()
  print 'r[i] = ... -->',
  print '%.2f s' % (t1 - t0 - t_overhead)
  print

def run():
  print ("Given af::shared<double> a and r = a.ref(), "
         "and using the following indexing in a loop:\n")
  for data_size in [1000, 10000, 100000, 1000000, 10000000]:
    n_repeats = 30000000 // data_size
    exercise(data_size, n_repeats)
  print 'OK'

if __name__ == '__main__':
  run()
