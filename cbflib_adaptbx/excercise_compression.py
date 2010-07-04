from scitbx.array_family import flex #mandatory import for cbflib_ext
from cbflib_ext import uncompress,compress,assert_equal
from libtbx.development.timers import Profiler

def create_random_data_with_gaussian_distribution(mu=0.0,sigma=1.0):
  slow=2048
  fast=2048
  total_sz = slow*fast
  random_data = flex.random_int_gaussian_distribution(total_sz,mu,sigma)
  random_data.reshape(flex.grid(slow,fast))
  return random_data

def basic_tests(verbose=True):
  initial_intdata = create_random_data_with_gaussian_distribution(0.0,100.0)

  #special deltas to test the compression algorithm
  addresses = [3,6,9,12,15,18]
  deltas = [-127,128,-32767,32768,-2147483647,2147483647]
  for x in xrange(6):
    initial_intdata[addresses[x]-1]=0
    initial_intdata[addresses[x]]=deltas[x]

  if verbose: P=Profiler("compress")
  array_shape = initial_intdata.focus()
  if verbose: print array_shape
  compressed = compress(initial_intdata)
  if verbose: print len(compressed)
  if verbose: P=Profiler("uncompress")
  decompressed_dat = uncompress(packed=compressed, fast=array_shape[1], slow=array_shape[0])

  if verbose: del P
  assert assert_equal(initial_intdata, decompressed_dat)

if __name__=="__main__":
  basic_tests(False)
  print "OK"
