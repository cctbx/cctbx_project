from __future__ import absolute_import, division, print_function
from six.moves import range
from scitbx.array_family import flex
from cbflib_adaptbx import uncompress,compress,assert_equal
from libtbx.development.timers import Profiler

def create_random_data_with_gaussian_distribution(mu=0.0,sigma=1.0):
  slow=2048
  fast=2048
  total_sz = slow*fast
  flex.set_random_seed(2048)
  random_data = flex.random_int_gaussian_distribution(total_sz,mu,sigma)
  random_data.reshape(flex.grid(slow,fast))
  return random_data

def basic_tests(verbose=True):
  initial_data = create_random_data_with_gaussian_distribution(0.0,100.0)

  #special deltas to test the compression algorithm
  addresses = [3,6,9,12,15,18]
  deltas = [-127,128,-32767,32768,-2147483647,2147483647]
  for x in range(6):
    initial_data[addresses[x]-1]=0
    initial_data[addresses[x]]=deltas[x]

  if verbose: P=Profiler("compress")
  array_shape = initial_data.focus()
  if verbose: print(array_shape)
  compressed_data = compress(initial_data)
  if verbose: print(len(compressed_data))
  if verbose: P=Profiler("uncompress")
  decompressed_data = uncompress(packed=compressed_data, fast=array_shape[1],
                                 slow=array_shape[0])

  if verbose: del P
  assert assert_equal(initial_data, decompressed_data)

if __name__=="__main__":
  basic_tests(False)
  print("OK")
